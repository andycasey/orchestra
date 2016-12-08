
"""
Measure the S_HK statistic as a proxy for stellar activity.
"""

import os
import logging
import numpy as np
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
THREADS = 8
DATA_DIR = os.path.join(cwd, "../data/spectra/")

# Find all spectra.
filenames = glob(os.path.join(DATA_DIR, "data/reduced/*/*_s1d_A.fits"))
N = len(filenames)
logger.info("There are {} HARPS spectra matching '*_s1d_A.fits'".format(N))

# Database credentials
with open(os.path.join(cwd, "../db/credentials.yaml"), "r") as fp:
    credentials = yaml.load(fp)

# Create a function to process all spectra in parallel.
def _process_stellar_activity(filename):
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

    connection = pg.connect(**credentials)
    cursor = connection.cursor()

    def _cleanup():
        cursor.close()
        connection.close()
        return None

    # Check to see if we have already measured the activity from this spectrum.
    cursor.execute(
        """SELECT EXISTS(SELECT 1 FROM stellar_activity WHERE date_obs=%s)""",
        (date_obs, ))
    exists, = cursor.fetchone()

    if exists and not OVERWRITE:
        logger.info(
            "Activity already measured for date_obs = '{}'; not overwriting"\
            .format(date_obs))
        return _cleanup()

    # Get the measured radial velocity from the database.
    cursor.execute(
        """ SELECT drs_ccf_rvc AS stellar_rv
            FROM obs WHERE date_obs = %s""", (date_obs, ))
    if 1 > cursor.rowcount:
        logger.warn(
            "No headers ingested for observation with date_obs = '{}'"\
            .format(date_obs))
        return _cleanup()
    
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
    if exists:
        cursor.execute(
            """ UPDATE stellar_activity
                SET s_hk = %s,
                    e_s_hk = %s,
                    filename = %s
                WHERE date_obs = %s """,
            (s_hk, e_s_hk, basename, date_obs))

    else:
        cursor.execute(
            """ INSERT INTO stellar_activity (date_obs, filename, s_hk, e_s_hk)
                VALUES (%s, %s, %s, %s)""",
            (date_obs, basename, s_hk, e_s_hk))

    connection.commit()

    return _cleanup()
    

for filename in filenames:
    _process_stellar_activity(filename)


raise a

pool = mp.Pool(THREADS)
result = pool.map_async(_process_stellar_activity, filenames).get()
pool.close()
pool.join()
