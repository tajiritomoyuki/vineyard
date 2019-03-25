#-*-coding:utf-8-*-
import os
import numpy as np
import glob
import h5py
import MySQLdb
from tqdm import tqdm
from itertools import product
from scipy.ndimage.morphology import binary_dilation
from joblib import Parallel, delayed

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

pattern = "CTL"

FFIdir = "/stingray/tess/data/FFI"

if pattern == "CTL":
    step2 = "/pike/pipeline/step2"
    step3 = "/pike/pipeline/step3"
elif pattern == "TIC":
    step2 = "/pike/pipeline/TIC2"
    step3 = "/pike/pipeline/TIC3"

def loadFFI(sector, camera, CCD):
    fitslist = glob.glob(os.path.join(FFIdir, "*%s-%s-%s-*ffic.fits" % (sector, camera, CCD)))
    fitslist.sort()
    return fitslist

def gather_hdf(sector):
    h5list = glob.glob(os.path.join(step2, "*_%s_1_?.h5" % sector))
    return h5list

def fits2pos(fitslist):
    x_center = np.array([])
    y_center = np.array([])
    for fitspath in tqdm(fitslist):
        #x, y座標を格納
        with fits.open(fitspath) as hdu:
            try:
                x = hdu[1].header["CRVAL1"]
                y = hdu[1].header["CRVAL2"]
            #何故かかけている時があるのでそういう場合は削除
            except:
                x = y = 0.
            x_center = np.hstack((x_center, x))
            y_center = np.hstack((y_center, y))
    return x_center, y_center

def pos2quality(x_center, y_center, sigma):
    #統計値を算出
    x_mean = np.nanmean(x_center)
    y_mean = np.nanmean(y_center)
    x_std = np.nanstd(x_center)
    y_std = np.nanstd(y_center)
    #sigma以上離れている点を1とする
    x_cond = np.logical_or(x_center > x_mean + sigma * x_std, x_center < x_mean - sigma * x_std)
    y_cond = np.logical_or(y_center > y_mean + sigma * y_std, y_center < y_mean - sigma * y_std)
    quality = np.logical_or(x_cond, y_cond)
    return quality

def make_quality_flag(sector, sigma=2.):
    #セクター中のすべてのqualityを統合
    quality_list = []
    print("caliculating quality...")
    for came, chi in product("1234", "1234"):
        print("camera%s chip%s" % (came, chi))
        #FFIを取得
        fitslist = loadFFI(sector, came, chi)
        #x, yの位置の時系列データを取得
        x_center, y_center = fits2pos(fitslist)
        #sector4のみ別処理
        if sector == 4:
            #差分を取得
            x_diff1 = np.diff(x_center)
            y_diff1 = np.diff(y_center)
            #差分でqualityを取得
            quality1 = pos2quality(x_diff1, y_diff1, sigma)
            #quality除去した差分を再取得
            x_diff2 = np.where(quality1, np.nan, x_diff1)
            y_diff2 = np.where(quality1, np.nan, y_diff1)
            #差分でqualityを再取得
            quality2 = pos2quality(x_diff2, y_diff2, 3.)
            #両者で少なくともひとつひっかかったqualityを集める
            quality = np.logical_or(quality1, quality2)
        else:
            #x, yの位置の時系列データからqualityを算出
            quality = pos2quality(x_center, y_center, sigma)
        quality_list.append(quality)
    #qualityが少なくとも1つ以上1であったらフラグを立てる
    quality_arr = np.sum(np.vstack(quality_list), axis=0)
    quality_arr = np.where(quality_arr > 0, 1, 0)
    #前後1点もqualityフラグを立てる
    quality_arr = binary_dilation(quality_arr)
    return quality_arr

def add_quality_flg(h5path, quality_arr):
    #path決める
    h5name = os.path.basename(h5path)
    dstpath = os.path.join(step3, h5name)
    #readhdf
    with h5py.File(h5path, "r") as fr:
        with h5py.File(dstpath, "w") as fw:
            #各アイテムごとにグループを作成してコピー
            group_list = ["header", "TPF", "LC", "APERTURE_MASK"]
            for group in group_list:
                #グループを作成
                fw_group = fw.create_group(group)
                #データコピー
                for item in fr[group].keys():
                    fr.copy("%s/%s" % (group, item), fw_group)
            #qualityデータを格納
            fw.create_dataset("TPF/QUALITY", data=quality_arr)
            fw.create_dataset("LC/QUALITY", data=quality_arr)

def main():
    #各セクターごとにクオリティフラグを作成
    # for sector in range(1, 6):
    sector = 7
    quality_arr = make_quality_flag(sector)
    #hdf集める
    h5list = gather_hdf(sector)
    #qualityを付与
    Parallel(n_jobs=32)(delayed(add_quality_flg)(h5path, quality_arr) for h5path in tqdm(h5list))

if __name__ == '__main__':
    main()
