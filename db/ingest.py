
"""
Ingest header information from HARPS data products.
"""

__author__ = "Andrew R. Casey <arc@ast.cam.ac.uk>"

import os
import multiprocessing as mp
import psycopg2 as pg
import yaml
from glob import glob
from astropy.io import fits


# Database credentials
dirname = os.path.dirname(__file__)
with open(os.path.join(dirname, "credentials.yaml"), "r") as fp:
    credentials = yaml.load(fp)

# Find the relevant *_bis_G2_A.fits files.
filenames = glob("data/spectra/data/reduced/*/*_bis_G2_A.fits")

# Get translators for the column/header names.
with open(os.path.join(dirname, "obs-columns.yaml"), "r") as fp:
    columns = yaml.load(fp)

    column_names = ", ".join(columns.values())
    fillers = ", ".join(["%s"] * len(columns))


def ingest(filename):
    """
    Ingest the headers from a reduced HARPS spectrum (*_bis_G2_A.fits product).

    :param filename:
        The local path of the reduced FITS file.
    """

    print("Ingesting from {}".format(filename))

    connection = pg.connect(**credentials)
    cursor = connection.cursor()

    with fits.open(filename) as image:
        cursor.execute(
            """SELECT EXISTS(SELECT 1 FROM obs WHERE date_obs=%s)""",
            (image[0].header["DATE-OBS"], ))
        exists, = cursor.fetchone()
        
        if not exists:
            values = [image[0].header.get(k, None) for k in columns.keys()]
            cursor.execute(
                """INSERT INTO obs ({column_names}) VALUES ({fillers})"""\
                .format(column_names=column_names, fillers=fillers), values)

    cursor.close()
    connection.commit()

    return None


# Ingest the headers from each file.
pool = mp.Pool(20)
r = pool.map_async(ingest, filenames).get()
pool.close()
pool.join()
