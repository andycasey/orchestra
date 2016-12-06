
"""
Fix the folder structure for extracted products from the HARPS pipeline.
"""

__author__ = "Andrew R. Casey <arc@ast.cam.ac.uk>"


import os
from glob import glob

cwd = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.realpath(os.path.join(cwd, "../data/spectra/"))

filenames = glob(os.path.join(DATA_DIR, "HARPS.*.fits")) \
          + glob(os.path.join(DATA_DIR, "HARPS.*.tbl"))
          
N = len(filenames)

for i, from_filename in enumerate(filenames):

    print("At {}/{}: {}".format(i + 1, N, from_filename))

    basename = os.path.basename(from_filename)
    date = basename.split("T")[0].split(".")[1]
    folder = os.path.join(DATA_DIR, "data/reduced/{}".format(date))
    if not os.path.exists(folder):
        os.makedirs(folder)

    to_filename = os.path.join(folder, basename)
    os.system('mv -v "{}" "{}"'.format(from_filename, to_filename))

# Fix the hierarchy so we just have data/spectra/reduced/YYYY-MM-DD/?
# TODO: This is an unusually big decision....

print("OK, start the ingest!")