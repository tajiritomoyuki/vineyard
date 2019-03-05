#-*-coding:utf-8-*-
import os
import numpy as np
import pandas as pd
from scipy import stats, signal
import glob
import cv2
import h5py
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm
from astropy.io import fits
from astroquery.mast import Catalogs
import matplotlib.pyplot as plt

step1 = "/pike/pipeline/TIC1"
step2 = "/pike/pipeline/TIC2"

def load_h5():
    h5list = glob.glob(os.path.join(step1, "*.h5"))
    return h5list

def load_data(f):
    time = np.array(f["TPF"]["TIME"])
    flux = np.array(f["TPF"]["ROW_CNTS"])
    cx = f["header"]["cx"].value
    cy = f["header"]["cy"].value
    Tmag = f["header"]["Tmag"].value
    return time, flux, cx, cy, Tmag

def determine_area_thresh(Tmag):
    #Tmagによってapertureに使用するpixelに制限をかける
    area_len = 9 - np.fix(Tmag / 2)
    #最大値を7*7に制限
    area_len = min(area_len, 7)
    #最小値を3*3に制限
    area_len = max(area_len, 3)
    return area_len ** 2

def trim_aperture(img, sigma, mid_val, Q_std, area_thresh):
    #しきい値を決める
    thresh = mid_val + sigma * Q_std
    #しきい値以下のものを0、他を1にする
    thimg = np.where(img > thresh, 1, 0).astype(np.uint8)
    #特徴検出
    contours = cv2.findContours(thimg, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)[1]
    #領域が大きすぎるもの小さすぎるものは排除
    contours = [contour for contour in contours if 4 <= cv2.contourArea(contour) <= area_thresh]
    return contours

def make_aperture(img, center, area_thresh=9):
    mid_val = np.nanmedian(img)
    img = np.nan_to_num(img)
    #統計量を求める。
    flat_img = np.ravel(img)
    Q1 = stats.scoreatpercentile(flat_img, 25)
    Q3 = stats.scoreatpercentile(flat_img, 75)
    Q_std = Q3 - Q1
    #星中心を算出
    center_tuple = tuple(np.round(center).astype(np.uint8))
    #3sigma以上の切り出し領域を求める
    contours = trim_aperture(img, 3, mid_val, Q_std, area_thresh)
    #4sigma以上の切り出し領域を求める
    contours.extend(trim_aperture(img, 4, mid_val, Q_std, area_thresh))
    for contour in contours:
        #中心が含まれているか確認
        test = cv2.pointPolygonTest(contour, center_tuple, False)
        if test >= 0:
            #apertureを作成
            aperture = np.zeros_like(img).astype(np.uint8)
            cv2.fillConvexPoly(aperture, points=contour, color=1)
            break
    #決めかねてしまう場合
    else:
        #中心含む4pixをapertureにする
        offset = np.array([[0.5, 0.5], [0.5, -0.5], [-0.5, 0.5], [-0.5, -0.5]])
        aperture_contour = np.round(center + offset).astype(np.int32)
        aperture = np.zeros_like(img).astype(np.uint8)
        cv2.fillConvexPoly(aperture, points=aperture_contour, color=1)
    return aperture

def make_background(img, center, aperture, sigma=1, bg_thresh=20):
    mid_val = np.nanmedian(img)
    img = np.nan_to_num(img)
    #しきい値を決める
    flat_img = np.ravel(img)
    Q1 = stats.scoreatpercentile(flat_img, 25)
    Q3 = stats.scoreatpercentile(flat_img, 75)
    Q_std = Q3 - Q1
    thresh = mid_val + sigma * Q_std
    #apertrueより外側数pixelの部分をbackgroundとする
    #backgroundの内側領域は1pixelで固定
    background_inner = cv2.dilate(aperture, np.ones((3, 3), np.uint8))
    for i in range(4):
        #backgroundの外側領域を決定
        kernel = np.ones((i * 2 + 5, i * 2 + 5), np.uint8)
        background_outer = cv2.dilate(aperture, kernel)
        #backgroundの候補pixelを求める
        bg_aperture = background_outer - background_inner
        #sigmaが1以上の点を除去する
        bg_aperture[img > thresh] = 0
        #backgroundのpixel数が十分あるか確認
        if np.sum(bg_aperture) > bg_thresh:
            break
    #backgroundが定義できない場合
    else:
        bg_aperture = np.zeros_like(aperture)
    return bg_aperture

def save(dstpath, fr, aperture, bkg_aperture, calibrated_flux, time, sap_flux):
    with h5py.File(dstpath, "w") as fw:
        #グループを作成
        fw_header = fw.create_group("header")
        fw_TPF = fw.create_group("TPF")
        fw_LC = fw.create_group("LC")
        fw.create_group("APERTURE_MASK")
        #データコピー
        for item in fr["header"].keys():
            fr.copy("header/%s" % item, fw_header)
        for item in fr["TPF"].keys():
            fr.copy("TPF/%s" % item, fw_TPF)
        # fr.copy("TPF/QUALITY", fw_LC)
        #データを格納
        fw.create_dataset("TPF/FLUX", data=calibrated_flux)
        fw.create_dataset("LC/TIME", data=time)
        fw.create_dataset("LC/SAP_FLUX", data=sap_flux)
        fw.create_dataset("APERTURE_MASK/FLUX", data=aperture)
        fw.create_dataset("APERTURE_MASK/FLUX_BKG", data=bkg_aperture)

def main(h5path):
    h5name = os.path.basename(h5path)
    with h5py.File(h5path, "r") as f:
        #データをロード
        time, flux, cx, cy, Tmag = load_data(f)
        height, width = flux.shape[1:]
        if height == 11 and width == 11:
            #画像の時間積分の中央値を取得
            img = np.nanmedian(flux, axis=0)
            #apertrueに使用するpixel数の上限を計算
            area_thresh = determine_area_thresh(Tmag)
            #apertureを計算
            aperture = make_aperture(img, (cx, cy), area_thresh=area_thresh)
            #backgroundを計算
            bkg_aperture = make_background(img, (cx, cy), aperture)
            #バックグランウドを引く
            bkg_frame = np.where(bkg_aperture == 1, flux, np.nan)
            #MMMalgorithmを使用
            bkg_arr = 3. * np.nanmedian(bkg_frame, axis=(1, 2)) - 2. * np.nanmean(bkg_frame, axis=(1, 2))
            calibrated_flux = flux - bkg_arr.reshape((bkg_arr.shape[0], 1, 1))
            #light curveを作る
            aperture_frame = np.where(aperture == 1, calibrated_flux, 0)
            sap_flux = np.sum(aperture_frame, axis=(1, 2))
            #出力
            dstpath = os.path.join(step2, h5name)
            save(dstpath, f, aperture, bkg_aperture, calibrated_flux, time, sap_flux)

if __name__ == '__main__':
    #ファイルリストを取得
    h5list = load_h5()
    # for h5path in tqdm(h5list):
    #     #try:
    #     main(h5path)
    #     #except:
    #     #    pass
    Parallel(n_jobs=32)(delayed(main)(h5path) for h5path in tqdm(h5list))
