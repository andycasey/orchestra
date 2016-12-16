
"""
Measure the S_HK statistic as a proxy for stellar activity.
"""

import os
import logging
import numpy as np
import multiprocessing as mp
import psycopg2 as pg
import yaml
from astropy.io import fits
from glob import glob

from orchestra.stellar_activity import shk_index as measure_stellar_activity

# Enable logging.
logger = logging.getLogger("orchestra")
logger.setLevel(logging.INFO) 

hdl = logging.StreamHandler()
hdl.setFormatter(logging.Formatter("%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(hdl)

# Set some configuration values.
cwd = os.path.dirname(__file__)
OVERWRITE = True
THREADS = 10
DATA_DIR = os.path.join(cwd, "../data/spectra/")

# Find all spectra.
filenames = glob(os.path.join(DATA_DIR, "data/reduced/*/*_s1d_A.fits"))
N = len(filenames)
logger.info("There are {} HARPS spectra matching '*_s1d_A.fits'".format(N))

# Database credentials
with open(os.path.join(cwd, "../db/credentials.yaml"), "r") as fp:
    credentials = yaml.load(fp)

# Create a function to process all spectra in parallel.
def _process_stellar_activity(filename, connection):
    """
    Measure the S_HK index as a proxy of stellar activity from a reduced HARPS
    spectrum, and ingest that measurement into the database.

    :param filename:
        The local path of a HARPS reduced data product (of the name *_s1d_A.fits)
    """

    basename = os.path.basename(filename)

    with fits.open(filename) as image:
        naxis = image[0].header["NAXIS1"]
        crval, cdelt = image[0].header["CRVAL1"], image[0].header["CDELT1"]
    
        date_obs = image[0].header["DATE-OBS"]
        object_name = image[0].header["OBJECT"]

        wavelength = np.arange(naxis) * cdelt + crval
        flux = image[0].data

    # Get the measured radial velocity from the database.
    cursor = connection.cursor()
    cursor.execute(
        """ SELECT drs_ccf_rvc AS stellar_rv
            FROM obs WHERE date_obs = %s""", (date_obs, ))
    if 1 > cursor.rowcount:
        logger.warn(
            "No headers ingested for observation with date_obs = '{}'"\
            .format(date_obs))
        cursor.close()
        return None
    
    rv, = cursor.fetchone()
    rv = float(rv)

    if not np.isfinite(rv):
        s_hk, e_s_hk = (np.nan, np.nan)

        logger.warn(
            "RV for date_obs = '{}' (filename: {}) is not finite"\
            .format(date_obs, basename))       
    else:
        s_hk, e_s_hk = measure_stellar_activity(wavelength, flux, rv)

        logger.info(
            "Measured S_HK = {:.2f} ({:.2e}) in '{}' from {} (date_obs = {})"\
            .format(s_hk, e_s_hk, object_name, basename, date_obs))

    # Insert or update this measurement in the database.
    try:
        cursor.execute(
            """ INSERT INTO stellar_activity (date_obs, filename, s_hk, e_s_hk)
                VALUES (%s, %s, %s, %s)
                ON CONFLICT DO NOTHING;""",
            (date_obs, basename, s_hk, e_s_hk))

    except (pg.IntegrityError, pg.DataError):
        logger.warn("Exception when ingesting {}".format(filename))
        connection.rollback()
        
    else:
        connection.commit()

    cursor.close()
    return None


def process_stellar_activity_in_parallel(*filenames):
    connection = pg.connect(**credentials)
    for filename in filenames:
        _process_stellar_activity(filename, connection)

    connection.commit()
    connection.close()


# Chunk it out.
pool = mp.Pool(THREADS)
s = int(np.ceil(float(N)/THREADS))

results = []
for t in range(THREADS):
    results.append(
        pool.apply_async(
            process_stellar_activity_in_parallel, filenames[t * s:(t + 1) * s]))

results = [each.get() for each in results]
pool.join()
pool.close()

