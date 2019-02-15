#coding:utf-8
import os
import glob
from astropy.io import fits
import subprocess
import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed

datadir = "/manta/tess/data/FFI"

def search(fitspath):
    #不正なfitsファイルがないか検索
    try:
        hdu = fits.open(fitspath)
        data = np.array(hdu[1].data)
        hdu.close()
    #もし不正なfitsファイルが有ったばあい、再ダウンロード
    except:
        fitsname = os.path.basename(fitspath)
        cmd = "curl -C - -L -o %s https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/%s" % (fitspath, fitsname)
        subprocess.run(cmd)

def main():
    fitslist = glob.glob(os.path.join(datadir, "*ffic.fits"))
    Parallel(n_jobs=16)(delayed(search)(fitspath) for fitspath in tqdm(fitslist))

if __name__ == '__main__':
    main()
