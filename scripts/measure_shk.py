
"""
Measure the S_HK statistic as a proxy for stellar activity.
"""

import os
from glob import glob

from code.stellar_activity import shk_index as measure_stellar_activity

cwd = os.path.dirname(__file__)

OVERWRITE = False
DATA_DIR = os.path.join(cwd, "../data/spectra/")
 
# GET THE GIT HASH.

# Find all spectra.
filenames = glob(os.path.join(DATA_DIR, "data/reduced/*/*_s1d_A.fits"))

# Fit each one in parallel.
raise a
# Ingest in serial.

# Ingest date_obs, filename, s_hk, e_shk, shk_version, git_hash

