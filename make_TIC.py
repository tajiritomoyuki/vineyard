#coding:utf-8
import os
import numpy as np
import pandas as pd
import pandas.io.sql as pdsql
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
datadir = "/manta/tess/data/FFI"

def get_TIC(Tmag_limit):
    with MySQLdb.connect(**data) as cursor:
        query = "select ID, ra, `dec` from TICv7s where Tmag < %s;" % Tmag_limit
        TICdf = pdsql.read_sql(query, cursor)
    return TICdf

def get_CTL(Tmag_limit):
    with MySQLdb.connect(**data) as cursor:
        query = "select ID from CTLv7 where Tmag < %s and `dec`<0;" % Tmag_limit
        CTLdf = pdsql.read_sql(query, cursor)
    return CTLdf

def main():
    Tmag_limit = 13
    TIC = get_TIC(Tmag_limit)
    CTL = get_CTL(Tmag_limit)


if __name__ == '__main__':
    main()
