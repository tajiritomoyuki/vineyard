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

data = {
    "user" : "root",
    "passwd" : "ocean",
    "host" : "localhost",
    "db" : "TESS"
}

class extractCTL(CTL):
    def __init__(self):
        self.CTL = None

    def get_CTL(self):
        con = MySQLdb.connect(**data)
        query = "select ID from CTLv7;"
        self.CTL = pdsql.read_sql(query, con)
        con.close()

    def extract(self, sector, camera, chip):
        with MySQLdb.connect(**data) as cursor:
            query = "select ID from chip%s_%s_%s" % (sector, camera, chip)
            cursor.execute(query)
            TIC = cursor.fetchall()
            for TID_row in TIC:
                TID = TID_row[0]
                if TID in self.CTL:
                    query = "delete from chip%s_%s_%s where ID=%s" % (sector, camera, chip, TID)
                    # cursor.execute(query)
                    print(query)

    def extract_all(self):
        ctx = mp.get_context("spawn")
        with ctx.Pool(5) as p:
            p.map(self.extract, tqdm(product("12345", "1234", "1234")))


if __name__ == '__main__':
    exCTL = extractCTL()
    exCTL.get_CTL()
    exCTL.extract_all()
