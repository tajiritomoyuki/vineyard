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

#FFIのディレクトリ
FFIdir = "/manta/tess/data/FFI"
#出力先のディレクトリ
outputdir = "/pike/pipeline/step1"

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

def save(TID, sector, camera, CCD, ra, dec, Tmag, x, y, cx, cy, time, flux):
    tpfname = "tess_%s_%s_%s_%s.h5" % (TID, sector, camera, CCD)
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

def main(sector, camera, CCD):
    #ID, ra, dec, Tmagデータを読み込み
    data = importTIC(sector, camera, CCD)
    #fitsファイルのリストを取得
    fitslist = loadFFI(sector, camera, CCD)
    #fitsファイルからtime, fluxを取得
    time, FFIflux = fits2data(fitslist)
    #wcsを取得
    wcs, bounds = get_wcs(fitslist[0])
    print("making h5file...")
    for TID, ra, dec, Tmag in tqdm(data):
        #ra, decからpixelを抽出
        x, y = radec2pix(ra, dec, wcs)
        #pixel情報からFFIを切り出し
        flux, cx, cy = cut(x, y, FFIflux)
        #出力
        save(TID, sector, camera, CCD, ra, dec, Tmag, x, y, cx, cy, time, flux)
    del FFIflux

if __name__ == '__main__':
    #Parallel(n_jobs=2)([delayed(main)(sector, camera, CCD) for sector, camera, CCD in product("12", "1234", "1234")])
    for sector, camera, CCD in product("6", "1234", "1234"):
        main(sector, camera, CCD)
