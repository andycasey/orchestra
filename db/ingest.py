
"""
Ingest header information from HARPS data products.
"""

__author__ = "Andrew R. Casey <arc@ast.cam.ac.uk>"

import os
import logging
import multiprocessing as mp
import psycopg2 as pg
import yaml
from glob import glob
from astropy.io import fits
from astropy.table import Table


# Database credentials
dirname = os.path.dirname(__file__)
with open(os.path.join(dirname, "credentials.yaml"), "r") as fp:
    credentials = yaml.load(fp)

#----------------------------------------------------------
#   Results from Phase 3 product searches
#----------------------------------------------------------

# Find relevant ESO Phase 3 product catalogs
phase3_filenames = glob("data/eso-phase3-products/*.fits")

def ingest_phase3_catalog(filename):
    """
    Ingest records from a table containing records from the ESO Phase 3 archive.
    These tables were produced by the `scripts/eso_search_phase3.py` file.

    :param filename:
        The local path of the catalog file.
    """

    table = Table.read(filename)

    print("Ingesting {} Phase 3 records from {}".format(len(table), filename))
    connection = pg.connect(**credentials)

    for record in table:
        
        cursor = connection.cursor()
        cursor.execute(
            """SELECT EXISTS(SELECT 1 FROM phase3_products WHERE arcfile=%s)""",
            (record["ARCFILE"], ))
        exists, = cursor.fetchone()

        if not exists:
            try:
                cursor.execute(
                    """INSERT INTO phase3_products (arcfile, object, ra, dec,
                        wavelength, snr, resolution, instrument, date_obs,
                        exptime, program_id, origfile, dataset) VALUES (%s, %s, 
                        %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)""",
                    [record[k] for k in record.dtype.names])

            except pg.IntegrityError:
                logging.exception("IntegrityError on {}/{}:\n{}\n".format(
                    filename, record["ARCFILE"], record))
                connection.rollback()

            else:
                connection.commit()

        cursor.close()
    connection.close()

    return None


# Ingest the headers from each file.
pool = mp.Pool(20)
r = pool.map_async(ingest_phase3_catalog, phase3_filenames).get()
pool.close()
pool.join()

raise a


#----------------------------------------------------------
#   Headers from reduced data products
#----------------------------------------------------------

# Find the relevant *_bis_G2_A.fits files.
obs_filenames = glob("data/spectra/data/reduced/*/*_bis_G2_A.fits")
obs_filenames += glob("data/spectra/HARPS*_bis_G2_A.fits")

# Get translators for the column/header names.
with open(os.path.join(dirname, "obs-columns.yaml"), "r") as fp:
    columns = yaml.load(fp)

    column_names = ", ".join(columns.values())
    fillers = ", ".join(["%s"] * len(columns))


def ingest_obs_headers(filename):
    """
    Ingest the headers from a reduced HARPS spectrum (*_bis_G2_A.fits product).

    :param filename:
        The local path of the reduced FITS file.
    """

    print("Ingesting observed headers from {}".format(filename))

    connection = pg.connect(**credentials)
    cursor = connection.cursor()

    with fits.open(filename) as image:
        cursor.execute(
            """SELECT EXISTS(SELECT 1 FROM obs WHERE date_obs=%s)""",
            (image[0].header["DATE-OBS"], ))
        exists, = cursor.fetchone()
        
        if not exists:
            try:
                values = [image[0].header.get(k, None) for k in columns.keys()]

            except fits.VerifyError:
                image.verify("fix")
                values = [image[0].header.get(k, None) for k in columns.keys()]                
            
            cursor.execute(
                """INSERT INTO obs ({column_names}) VALUES ({fillers})"""\
                .format(column_names=column_names, fillers=fillers), values)

    cursor.close()
    connection.commit()

    return None


# Ingest the headers from each file.
pool = mp.Pool(20)
r = pool.map_async(ingest_obs_headers, obs_filenames).get()
pool.close()
pool.join()
