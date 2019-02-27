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
    "user" : "root",
    "passwd" : "ocean",
    "host" : "localhost",
    "db" : "TESS"
}
datadir = "/manta/tess/data/FFI"

class Register():
    def __init__(self, sector, camera, chip):
        self.sector = sector
        self.camera = camera
        self.chip = chip
        self.wcs = None
        self.bounds = None

    def get_wcs(self):
        fitslist = glob.glob(os.path.join(datadir, "*%s-%s-%s*.fits" % (self.sector, self.camera, self.chip)))
        hdu = fits.open(fitslist[0])
        self.wcs = WCS(hdu[1].header)
        self.bounds = hdu[1].data.shape
        hdu.close()

    def check_coord(self, row):
        print(row, self.sector)
        coord = SkyCoord(ra, dec, unit="deg")
        try:
            px, py = coord.to_pixel(wcs)
            if (px < bounds[0]) and (px > 0) and (py < bounds[1]) and (py > 0):
                # with MySQLdb.connect(**data) as cursor:
                query = "INSERT INTO chip%s_%s_%s SELECT * FROM TICv7s_has_key WHERE ID=%s" % (self.sector, self.camera, self.chip, ID)
                print(query)
                    # cursor.execute(query)
        except NoConvergence:
            pass

    def register(self, TICdf):
        ctx = mp.get_context("spawn")
        with ctx.Pool(mp.cpu_count()) as p:
            print("now")
            p.map(self.check_coord, TICdf.iterrows())

def get_TIC(Tmag_limit):
    con = MySQLdb.connect(**data)
    query = "select ID, ra, `dec` from TICv7s where Tmag < %s limit 20000;" % Tmag_limit
    TICdf = pdsql.read_sql(query, con)
    con.close()
    return TICdf

def get_CTL(Tmag_limit):
    con = MySQLdb.connect(**data)
    query = "select ID from CTLv7 where Tmag < %s and `dec`<0;" % Tmag_limit
    CTLdf = pdsql.read_sql(query, con)
    con.close()
    return CTLdf

def omit_dupilication(TICdf, CTLdf):
    iterator = TICdf["ID"].__iter__()
    #https://nb4799.neu.edu/wordpress/?p=783
    ctx = mp.get_context("spawn")
    manager = mp.Manager()
    CTLdf_shared = manager.list(CTLdf.values.tolist())
    with ctx.Pool(mp.cpu_count()) as p:
        v = p.map(functools.partial(omit, CTLdf=CTLdf_shared), iterator)
    for TID in v:
        print(TID)

def main():
    Tmag_limit = 13
    TICdf = get_TIC(Tmag_limit)
    CTLdf = get_CTL(Tmag_limit)
    print("load")
    for sector, camera, chip in product("12345", "1234", "1234"):
        regi = Register(sector, camera, chip)
        regi.get_wcs()
        regi.register(TICdf)

if __name__ == '__main__':
    main()
