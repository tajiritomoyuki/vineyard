#coding:utf-8
import os
import glob
from astropy.io import fits
import subprocess
import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed
from itertools import product

datadir = "/stingray/tess/data/FFI"

def check_lack(date, sector, camera, chip, yonmoji):
    tarname = "%s-s%04d-%s-%s-%s-s_ffic.fits" % (date, sector, camera, chip, yonmoji)
    tarpath = os.path.join(datadir, tarname)
    if not os.path.exists(tarpath):
        cmd = "curl -C - -L -o %s https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/%s" % (tarpath, tarname)
        subprocess.run(cmd, shell=True)

def main():
    #各セクターごとに足りないFFIがないか検索
    # for sector in range(1, 6):
    sector = 10
    fitslist = glob.glob(os.path.join(datadir, "*s%04d*ffic.fits" % sector))
    datelist = list(set([os.path.basename(fitspath).split("-")[0] for fitspath in fitslist]))
    yonmoji = fitslist[0].split("-")[4]
    Parallel(n_jobs=20)(delayed(check_lack)(date, sector, camera, chip, yonmoji) for date, camera, chip in product(datelist, "1234", "1234"))


if __name__ == '__main__':
    main()
