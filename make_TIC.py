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

def get_TIC(Tmag_limit):
    con = MySQLdb.connect(**data)
    query = "select ID, ra, `dec` from TICv7s where Tmag < %s limit 50;" % Tmag_limit
    TICdf = pdsql.read_sql(query, con)
    con.close()
    return TICdf

def get_CTL(Tmag_limit):
    con = MySQLdb.connect(**data)
    query = "select ID from CTLv7 where Tmag < %s and `dec`<0;" % Tmag_limit
    CTLdf = pdsql.read_sql(query, con)
    con.close()
    return CTLdf

def make_iter(TICdf, CTLdf):
    for TID in TICdf["ID"].__iter__():
        yield TID, CTLdf

def omit(TID, CTLdf):
    return TID

def omit_dupilication(TICdf, CTLdf):
    iterator = make_iter(TICdf, CTLdf)
    ctx = mp.get_context("spawn")
    with ctx.Pool(mp.cpu_count()) as p:
        v = p.imap_unordered(omit, iterator)
    for TID in v:
        print(TID)

def main():
    Tmag_limit = 13
    TICdf = get_TIC(Tmag_limit)
    CTLdf = get_CTL(Tmag_limit)
    omit_dupilication(TICdf, CTLdf)

if __name__ == '__main__':
    main()
