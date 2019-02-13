#coding:utf-8
import os
import glob
from astropy.io import fits
import csv
import numpy as np
from tqdm import tqdm

csvpath = "/home/tajiri/tess/pipeline/incorrect3.csv"
datadir = "/pike/tess/data/FFI"
fitslist = glob.glob(os.path.join(datadir, "*s0003*.fits"))
with open(csvpath, "w") as f:
    writer = csv.writer(f, lineterminator="\n")
    for fitspath in tqdm(fitslist):
        try:
            hdu = fits.open(fitspath)
            data = np.array(hdu[1].data)
            hdu.close()
        except:
            writer.writerow([fitspath])
