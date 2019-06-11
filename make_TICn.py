#coding:utf-8
import os
import numpy as np
import pandas as pd
import pandas.io.sql as pdsql
import glob
from tqdm import tqdm
from itertools import product
import MySQLdb
import multiprocessing as mp

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS, NoConvergence

data = {
    "user" : "fisher",
    "passwd" : "atlantic",
    "host" : "133.11.231.118",
    "db" : "TESS"
}
datadir_org = "/%s/tess/data/FFI"

class Register():
    def __init__(self, sector, camera, chip):
        self.sector = sector
        self.camera = camera
        self.chip = chip
        self.wcs = None
        self.bounds = None

    def get_wcs(self):
        #sectorで条件分岐
        if int(self.sector) < 7:
            datadir = datadir_org % "manta"
        else:
            datadir = datadir_org % "stingray"
        fitslist = glob.glob(os.path.join(datadir, "*%s-%s-%s*.fits" % (self.sector, self.camera, self.chip)))
        hdu = fits.open(fitslist[0])
        self.wcs = WCS(hdu[1].header)
        self.bounds = hdu[1].data.shape
        hdu.close()

    def check_coord(self, row):
        ID, ra, dec = row
        coord = SkyCoord(ra, dec, unit="deg")
        try:
            px, py = coord.to_pixel(self.wcs)
            if (px < self.bounds[0]) and (px > 0) and (py < self.bounds[1]) and (py > 0):
                with MySQLdb.connect(**data) as cursor:
                    query = "INSERT INTO chip%s_%s_%s SELECT * FROM TICv7n_has_key WHERE ID=%s" % (self.sector, self.camera, self.chip, ID)
                    cursor.execute(query)
        except NoConvergence:
            pass

    def register(self, TICresult):
        ctx = mp.get_context("spawn")
        with ctx.Pool(10) as p:
            p.map(self.check_coord, tqdm(TICresult))

def get_TIC(Tmag_limit):
    with MySQLdb.connect(**data) as cursor:
        query = "select ID, ra, `dec` from TICv7n where Tmag < %s;" % Tmag_limit
        cursor.execute(query)
        result = cursor.fetchall()
    return result

def get_CTL(Tmag_limit):
    con = MySQLdb.connect(**data)
    query = "select ID from CTLv7 where Tmag < %s and `dec`<0;" % Tmag_limit
    CTLdf = pdsql.read_sql(query, con)
    con.close()
    return CTLdf

def omit_dupilication(TICdf, CTLdf):
    iterator = TICdf["ID"].__iter__()
    ctx = mp.get_context("spawn")
    manager = mp.Manager()
    CTLdf_shared = manager.list(CTLdf.values.tolist())
    with ctx.Pool(mp.cpu_count()) as p:
        v = p.map(functools.partial(omit, CTLdf=CTLdf_shared), iterator)
    for TID in v:
        print(TID)

def main():
    Tmag_limit = 13
    TICresult = get_TIC(Tmag_limit)
    for sector, camera, chip in product("8", "1", "1234"):
        regi = Register(sector, camera, chip)
        regi.get_wcs()
        regi.register(TICresult)

if __name__ == '__main__':
    main()
