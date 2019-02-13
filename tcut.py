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

rootdir = os.path.abspath(os.path.join(__file__))
#データの元ディレクトリ
datadir = "/pike/tess/data"
#FFIのディレクトリ
FFIdir = os.path.join(datadir, "FFI")
#出力先のディレクトリ
outputdir = os.path.join(datadir, "tesscut")

#SQLにログインするためのデータ
sql_data = {
      "user" : "root",
    "passwd" : "ocean",
      "host" : "localhost",
        "db" : "TESS",
}

def importTIC(sector, camera, CCD):
    with MySQLdb.connect(**sql_data) as cursor:
        query = "select ID, ra, `dec`, Tmag from CTLchip%s_%s_%s;" % (sector, camera, CCD)
        cursor.execute(query)
        result = cursor.fetchall()
    return result

def loadFFI(sector, camera, CCD):
    fitslist = glob.glob(os.path.join(FFIdir, "*%s-%s-%s-*ffic.fits" % (sector, camera, CCD)))
    fitslist.sort()
    return fitslist

def fits2data(fitslist):
    t_list = []
    f_list = []
    print("loading fits file...")
    for i, fitspath in enumerate(tqdm(fitslist)):
        with fits.open(fitspath) as hdu:
            #時間
            time_arr = hdu[1].header["TSTART"]
            #flux
            flux_arr = np.array(hdu[1].data)
            t_list.append(time_arr)
            f_list.append(flux_arr)
    time = np.array(t_list)
    flux = np.stack(tuple(f_list), axis=0)
    del t_list
    del f_list
    return time, flux

def get_wcs(fitsfile):
    hdu = fits.open(fitsfile)
    wcs = WCS(hdu[1].header)
    bounds = hdu[1].data.shape
    hdu.close()
    return wcs, bounds

def make_quality_flag_(x_center, y_center, sigma=2.):
    x_mean = np.nanmean(x)
    y_mean = np.nanmean(y)
    x_std = np.nanstd(x)
    y_std = np.nanstd(y)
    #sigma以上離れている点を1とする
    x_cond = np.logical_or(x > x_mean + sigma * x_std, x < x_mean - sigma * x_std)
    y_cond = np.logical_or(y > y_mean + sigma * y_std, y < y_mean - sigma * y_std)
    quality = np.logical_or(x_cond, y_cond)
    return quality

def make_quality_flag(sector, camera, CCD, sigma=2.):
    #カメラ中のすべてのqualityを統合
    quality_list = []
    print("caliculating quality...")
    for came, chi in product("1234", "1234"):
        fitslist = loadFFI(sector, came, chi)
        x_center = np.array([])
        y_center = np.array([])
        print("camera%s chip%s..." % (came, chi))
        for fitspath in tqdm(fitslist):
            #x, y座標を格納
            with fits.open(fitspath) as hdu:
                x = hdu[1].header["CRVAL1"]
                y = hdu[1].header["CRVAL2"]
                x_center = np.hstack((x_center, x))
                y_center = np.hstack((y_center, y))
        #統計値を算出
        x_mean = np.nanmean(x_center)
        y_mean = np.nanmean(y_center)
        x_std = np.nanstd(x_center)
        y_std = np.nanstd(y_center)
        #sigma以上離れている点を1とする
        x_cond = np.logical_or(x_center > x_mean + sigma * x_std, x_center < x_mean - sigma * x_std)
        y_cond = np.logical_or(y_center > y_mean + sigma * y_std, y_center < y_mean - sigma * y_std)
        quality = np.logical_or(x_cond, y_cond)
        quality_list.append(quality)
    quality_arr = np.sum(np.vstack(quality_list), axis=0)
    quality_arr = np.where(quality_arr > 0, 1, 0)
    quality_arr = binary_dilation(quality_arr)
    return quality_arr

def radec2pix(ra, dec, wcs):
    coord = SkyCoord(ra, dec, unit="deg")
    px, py = coord.to_pixel(wcs)
    return px, py

def cut(x, y, FFIflux, size=(11, 11)):
    x_int = np.round(x).astype(np.int32)
    y_int = np.round(y).astype(np.int32)
    height = (size[0] - 1) // 2
    width = (size[1] - 1) // 2
    flux = FFIflux[:, y_int-height:y_int+height+1, x_int-width:x_int+width+1]
    cx = width + x - x_int
    cy = height + y - y_int
    return flux, cx, cy

def save(TID, sector, camera, CCD, ra, dec, Tmag, x, y, cx, cy, time, flux, quality):
    tpfname = "tpf_%s_%s_%s_%s.h5" % (TID, sector, camera, CCD)
    tpfpath = os.path.join(outputdir, tpfname)
    with h5py.File(tpfpath, "w") as f:
        f.create_group("header")
        f.create_group("TPF")
        f.create_dataset("header/TID", data=TID)
        f.create_dataset("header/sector", data=sector)
        f.create_dataset("header/camera", data=camera)
        f.create_dataset("header/CCD", data=CCD)
        f.create_dataset("header/ra", data=ra)
        f.create_dataset("header/dec", data=dec)
        f.create_dataset("header/Tmag", data=Tmag)
        f.create_dataset("header/x", data=x)
        f.create_dataset("header/y", data=y)
        f.create_dataset("header/cx", data=cx)
        f.create_dataset("header/cy", data=cy)
        f.create_dataset("TPF/TIME", data=time)
        f.create_dataset("TPF/ROW_CNTS", data=flux)
        f.create_dataset("TPF/QUALITY", data=quality)

def main(sector, camera, CCD):
    #ID, ra, dec, Tmagデータを読み込み
    data = importTIC(sector, camera, CCD)
    #fitsファイルのリストを取得
    fitslist = loadFFI(sector, camera, CCD)
    #fitsファイルからtime, fluxを取得
    time, FFIflux = fits2data(fitslist)
    #wcsを取得
    wcs, bounds = get_wcs(fitslist[0])
    #quality_flagを作成
    quality = make_quality_flag(sector, camera, CCD)
    print("making h5file...")
    for TID, ra, dec, Tmag in tqdm(data):
        #ra, decからpixelを抽出
        x, y = radec2pix(ra, dec, wcs)
        #pixel情報からFFIを切り出し
        flux, cx, cy = cut(x, y, FFIflux)
        #出力
        save(TID, sector, camera, CCD, ra, dec, Tmag, x, y, cx, cy, time, flux, quality)
    del FFIflux

if __name__ == '__main__':
    #Parallel(n_jobs=2)([delayed(main)(sector, camera, CCD) for sector, camera, CCD in product("12", "1234", "1234")])
    for sector, camera, CCD in product("34", "1234", "1234"):
        main(sector, camera, CCD)
