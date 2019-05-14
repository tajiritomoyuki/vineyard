#coding:utf-8
import os
import numpy as np
import glob
from tqdm import tqdm
from itertools import product
import MySQLdb
from joblib import Parallel, delayed

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS, NoConvergence
data = {
    "user" : "root",
    "passwd" : "ocean",
    "host" : "localhost",
    "db" : "TESS"
}
datadir_org = "/%s/tess/data/FFI"

def process(s, x, y, result):
    #sectorで条件分岐
    if int(s) < 7:
        datadir = datadir_org % "manta"
    else:
        datadir = datadir_org % "stingray"
    fitslist = glob.glob(os.path.join(datadir, "*%s-%s-%s*.fits" % (s, x, y)))
    hdu = fits.open(fitslist[0])
    wcs = WCS(hdu[1].header)
    bounds = hdu[1].data.shape
    hdu.close()
    with MySQLdb.connect(**data) as cursor:
        for ID, ra, dec in tqdm(result):
            coord = SkyCoord(ra, dec, unit="deg")
            try:
                px, py = coord.to_pixel(wcs)
            except NoConvergence:
                continue
            if (px < bounds[0]) and (px > 0) and (py < bounds[1]) and (py > 0):
                query = "INSERT INTO CTLchip%s_%s_%s SELECT * FROM CTLv7_has_key WHERE ID=%s" % (s, x, y, ID)
                cursor.execute(query)

def main():
    with MySQLdb.connect(**data) as cursor:
        query = "select ID, ra, `dec` from CTLv7 where `dec`<30;"
        cursor.execute(query)
        result = cursor.fetchall()
    Parallel(n_jobs=4)([delayed(process)(s, x, y, result) for s, x, y in product("9", "1234", "1234")])
    # for s, x, y in product("7", "1234", "1234"):
    #     process(s, x, y, result)


if __name__ == '__main__':
    main()
