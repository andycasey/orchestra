Make measurements on HARPS spectra
==================================

This page provides instructions on how to write code that will measure something on all the HARPS spectra, then upload the result of that measurement to the `harps.guru` database.
This guide is incomplete. For example, you should contact Andy to get him to actually *run* your code, and this example assumes that you have already set up some table schema in
`db/schema.sql`.

With those assumptions in mind, here is some example code:
````python


"""
A template script for measuring something on all HARPS data, and uploading the
result to the harps.guru database.
"""

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import os
import logging
import numpy as np
import multiprocessing as mp
import psycopg2 as pg
import yaml
from astropy.io import fits
from glob import glob

# Enable logging.
logger = logging.getLogger("orchestra")
logger.setLevel(logging.INFO) 

hdl = logging.StreamHandler()
hdl.setFormatter(logging.Formatter("%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(hdl)

# Set some configuration values.
THREADS = 10
cwd = os.path.dirname(__file__)
DATA_DIR = os.path.join(cwd, "../data/spectra/")

# Find all spectra.
filenames = glob(os.path.join(DATA_DIR, "data/reduced/*/*_s1d_A.fits"))
N = len(filenames)
logger.info("There are {} HARPS spectra matching '*_s1d_A.fits'".format(N))

# Database credentials
with open(os.path.join(cwd, "../db/credentials.yaml"), "r") as fp:
    credentials = yaml.load(fp)


# Create a function to process a single spectrum.
def _process_harps_spectrum(filename, connection):
    """
    Measure something from a spectrum.
    """

    with fits.open(filename) as image:
        object_name = image[0].header["OBJECT"]

    cursor = connection.cursor()
    cursor.execute(
        """ INSERT INTO my_objects (object_name)
            VALUES (%s)
            ON CONFLICT DO NOTHING;""",
        (object_name, ))
    cursor.close()
    return None


def process_harps_spectra(*filenames):
    connection = pg.connect(**credentials)
    for filename in filenames:
        _process_harps_spectrum(filename, connection)

    connection.commit()
    connection.close()


# Chunk it out.
pool = mp.Pool(THREADS)
s = int(np.ceil(float(N)/THREADS))

results  =  [pool.apply_async(process_harps_spectra, filenames[t * s:(t + 1) * s]) \
                for t in range(THREADS)]      

results = [result.get() for result in results]
pool.join()
pool.close()
````
