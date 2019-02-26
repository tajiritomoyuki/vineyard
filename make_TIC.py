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
import functools

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

def omit(TID, CTLdf):
    return TID

def omit_dupilication(TICdf, CTLdf):
    iterator = TICdf["ID"].__iter__()
    #https://nb4799.neu.edu/wordpress/?p=783
    ctx = mp.get_context("spawn")
    manager = mp.Manager()
    CTLdf_shared = manager.list(CTLdf)
    with ctx.Pool(mp.cpu_count()) as p:
        v = p.imap_unordered(functools.partial(omit, CTLdf=CTLdf_shared), iterator)
    for TID in v:
        print(TID)

def main():
    Tmag_limit = 13
    TICdf = get_TIC(Tmag_limit)
    CTLdf = get_CTL(Tmag_limit)
    omit_dupilication(TICdf, CTLdf)

if __name__ == '__main__':
    main()
