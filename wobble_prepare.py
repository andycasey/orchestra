
""" 
Prepare HARPS data for Wobble. 
"""

import logging
import os
import argparse
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astroquery.gaia import Gaia

from harps.client import Harps

parser = argparse.ArgumentParser(description="""
Search for, retrieve, and prepare HARPS data for use with Wobble.

Given a provided source position, this script will perform a cone search against
Gaia DR2 (Brown et al. 2018) and use Gaia observables to select likely spectra
of this source from the HARPS archive. Quality constraints as to what should be
considered a good spectrum can be provided by the user. This script will then
retrieve all selected HARPS spectra, and process them into a HDF5 file that can
be used with Wobble (Bedell et al. 2019).
""")

def ValidFloatActionFactory(lower, upper):
    """ A factory for actions that will validate float inputs. """
    class ValidateAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            if not (upper >= values >= lower):
                raise ValueError(f"{self.dest} must be between [{lower}, {upper}]")
            setattr(namespace, self.dest, values)
    return ValidateAction

parser.add_argument("ra", metavar="ra",
                    type=float, action=ValidFloatActionFactory(0, 360),
                    help="right ascension of the source [degrees]")
parser.add_argument("dec", metavar="dec",
                    type=float, action=ValidFloatActionFactory(-90, 90),
                    help="declination of the source [degrees]")
parser.add_argument("--name", metavar="source_name", type=str,
                    help="the name of the source (used for local naming only)")
parser.add_argument("--radius", metavar="radius", type=float,
                    action=ValidFloatActionFactory(0, np.inf),
                    default=5 * u.arcsecond,
                    help="cone search radius [arcseconds] to use for external "\
                         "services (e.g., Gaia, Simbad)")
parser.add_argument("--eso-credentials-path", metavar="eso_credentials_path",
                    default="eso_credentials.yml",
                    help="local path containing ESO credentials")
parser.add_argument("--working-directory", metavar="working_directory",
                    help="local working directory for files")
parser.add_argument("--verbose", "-v", action="store_true",
                    help="verbose logging")
args = parser.parse_args()

# Prepare logging.
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG if args.verbose else logging.INFO)
ch.setFormatter(logging.Formatter("%(message)s"))
logger.addHandler(ch)

# Step 1: Gaia search.
input_coord = SkyCoord(ra=args.ra, dec=args.dec, unit=(u.degree, u.degree),
                       frame="icrs") # todo: allow user-specified frame

logger.info(f"Querying Gaia at {input_coord}\nwith radius {args.radius}")

j = Gaia.cone_search(input_coord, args.radius)
gaia_results = j.get_results()
G = len(gaia_results)

if G < 1:
    raise NotImplementedError

if G > 1:
    raise NotImplementedError

#python wobble_prepare.py 344.36658333333327 20.768833333333333
gaia_result = gaia_results[0]
if args.name is None:
    args.name = gaia_result["designation"].decode()

logger.info(f"Source identified as {gaia_result['designation'].decode()}")

gaia_ra, gaia_dec = (gaia_result["ra"], gaia_result["dec"])

# Step 3: Search HARPS.
logger.info("Logging in to ESO archive..")
harps = Harps(args.eso_credentials_path)

# todo: check that 'box' should actually be in degrees
harps_results = harps.query_position(gaia_ra, gaia_dec,
                                     box=f"{args.radius.to(u.deg)}")

logger.info(f"HARPS archive returned {len(harps_results)} rows:")
logger.info(harps_results)

if len(harps_results) < 1:
    raise NotImplementedError

# Step 4: Download files.
# todo: don't limit just to first 10 -- this is for testing only
rid, paths = harps.get_dataset_identifiers(harps_results["dataset_identifier"][:10])
N = len(paths)

if args.working_directory is None:
    args.working_directory = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        args.name)

os.makedirs(args.working_directory, exist_ok=True)
logger.info(f"Downloading {N} files to {args.working_directory}")

headers = dict()
header_keys = (
    "HIERARCH ESO DRS CAL TH FILE", 
    "HIERARCH ESO DRS DRS CCF RVC", 
    "HIERARCH ESO DRS BERV", 
    "HIERARCH ESO DRS BJD",
    "HIERARCH ESO DRS BERVMX"
)

for i, remote_path in enumerate(paths, start=1):
    local_path = os.path.join(args.working_directory, os.path.basename(remote_path))
    logger.info(f"  ({i}/{N}): {remote_path} -> {local_path}")
    if os.path.exists(local_path):
        logger.info("  -> Skipping")

    else:
        harps.get_remote_path(remote_path, local_path)

    # If it's a FITS file, open it and get the wavelength calibration file.
    if local_path.lower().endswith(".fits"):
        with fits.open(local_path) as image:
            headers[remote_path] = \
                dict([(k, image[0].header.get(k, None)) for k in header_keys])

# Download the unique calibration files.
calibration_remote_paths = []
for remote_path, h in headers.items():
    for k, v in h.items():
        if k == "HIERARCH ESO DRS CAL TH FILE":
            calibration_remote_paths.append(f"SAF+{v[:-12]}.tar")
calibration_remote_paths = list(set(calibration_remote_paths))

cal_rid, cal_paths = harps.get_dataset_identifiers(calibration_remote_paths)

logger.info(f"Downloading {len(cal_paths)} calibration files")

for i, remote_path in enumerate(cal_paths, start=1):
    local_path = os.path.join(args.working_directory, os.path.basename(remote_path))
    logger.info(f"  ({i}/{len(cal_paths)}): {remote_path} -> {local_path}")
    if os.path.exists(local_path):
        logger.info("  -> Skipping")
        continue

    harps.get_remote_path(remote_path, local_path)

# Step 5: Extract all downloaded TAR files.


# Step 6: Run Bedell script.


print("Fin")