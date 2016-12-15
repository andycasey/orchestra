
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

cwd = os.path.dirname(__file__)

DEBUG = False
THREADS = 10
DATA_DIR = os.path.realpath(os.path.join(cwd, "../data/spectra/"))

# Database credentials
with open(os.path.join(cwd, "../db/credentials.yaml"), "r") as fp:
    credentials = yaml.load(fp)

# Find the relevant *_bis_G2_A.fits files.
obs_filenames = glob(os.path.join(DATA_DIR, "data/reduced/*/*_bis_G2_A.fits"))
N = len(obs_filenames)

# Get translators for the column/header names.
with open(os.path.join(cwd, "../db/obs-columns.yaml"), "r") as fp:
    columns = yaml.load(fp)


def _ingest_obs_headers(filename):
    """
    Ingest the headers from a reduced HARPS spectrum (*_bis_G2_A.fits product).

    :param filename:
        The local path of the reduced FITS file.
    """
    connection = pg.connect(**credentials)
    cursor = connection.cursor()

    required_keys = ("drs_ccf_rvc", "mjd_obs", "drs_dvrms")
    print("Ingesting {}".format(filename))

    try:
        with fits.open(filename) as image:
            cursor.execute(
                """SELECT EXISTS(SELECT 1 FROM obs WHERE date_obs=%s)""",
                (image[0].header["DATE-OBS"], ))
            exists, = cursor.fetchone()
            
            if not exists:

                keys = []
                values = []
                for header_key, sql_key in columns.items():
                    try:
                        value = image[0].header.get(header_key, None)

                    except fits.VerifyError:
                        print("Couldn't parse key '{}' from {}".format(
                            header_key, filename))

                        if sql_key == "date_obs":
                            # If this is the case, we are just fucked.
                            raise

                    else:
                        keys.append(sql_key)

                        # Some 2004 data fail because they record NaNs like this:
                        if value == "NaN." \
                        or (sql_key in required_keys and value is None):
                            value = "NaN"

                        values.append(value)

                # Sometimes the DATE-OBS and filename can be a millisecond different
                # from each other, so we are storing the filename stub here too.
                keys += ["filename"]
                values += ["/".join(filename.split("/")[-2:])]
                
                for required_key in required_keys:
                    if required_key not in keys:
                        keys.append(required_key)
                        values.append(np.nan)
                        print("ADDED FOR ", required_key)

                # Prepare inputs for the request.
                fillers = ", ".join(["%s"] * len(keys))
                column_names = ", ".join(keys)
                
                try:
                    cursor.execute(
                        """INSERT INTO obs ({column_names}) VALUES ({fillers})"""\
                        .format(column_names=column_names, fillers=fillers), 
                        values)

                except (pg.IntegrityError, pg.DataError):
                    logging.exception("Ingestion error on {}:".format(filename))
                    connection.rollback()
                    if DEBUG: raise
                    
                else:
                    connection.commit()

    except IOError:
        logging.exception("IOError while opening {}".format(filename))

    cursor.close()
    connection.close()

    return None


# Failure on at least:
#/media/My Book/research/orchestra/data/spectra/data/reduced/2004-10-05/HARPS.2004-10-05T03:16:30.268_bis_G2_A.fits:


#print("There are {} files to check (and potentially ingest)".format(N))
#for filename in obs_filenames:
#    _ingest_obs_headers(filename)


# Ingest the headers from each file.
pool = mp.Pool(THREADS)
result = pool.map_async(_ingest_obs_headers, obs_filenames).get()
pool.close()
pool.join()
