
"""
Find HARPS data from the ESO archive and produce a script to download the data.
"""

__author__ = "Andrew R. Casey <arc@ast.cam.ac.uk>"


# CRITICAL NOTE:  
#    You will need to authenticate with ESO in a Python terminal before running
#    this script. Here's how:

#    >> from astroquery.eso import Eso as ESO

#    >> eso = ESO()
#    >> eso.login("MY_USERNAME", store_password=True)

#    This will store your password locally so you don't need to provide it in
#    future sessions.

# This script also requires 'keyring.alt' package (available through pip)

import os
import re
import time
from astropy.table import Table, vstack
from astroquery.eso import Eso as ESO 
from glob import glob

# Login to ESO.
eso = ESO()
eso.login("andycasey")

cwd = os.path.dirname(os.path.realpath(__file__))
data_dir = "{}/../data/harps-spectra/".format(cwd)
eso.cache_location = "{}/cache/".format(cwd)
# See https://github.com/astropy/astroquery/issues/800
eso._survey_list = ["HARPS"]
eso.ROW_LIMIT = 10000

no_data = []
no_records = []
local_catalog = Table.read("data/HARPS_all.csv")
targets = list(set(local_catalog["Name"]))[1:]

N = len(targets)
M = 0

for i, target in enumerate(targets):

    eso_catalog_path = os.path.join(data_dir, "eso-{}.fits".format(target))
    
    if os.path.exists(eso_catalog_path):
        print("Skipping {}".format(target))
        continue

    # Search the ESO archive for HARPS data.
    r = eso.query_survey("HARPS", target=target)
    if r is None:
        print("No results found for star name {}".format(target))
        no_records.append(target)
        continue

    # For some reason we get a lot of shit that we don't want.
    # Require:
    keep = (r["Instrument"] == "HARPS") \
         * (r["R (&lambda;/&delta;&lambda;)"] == 115000)
    eso_catalog = r[keep]

    # Search again but return in HTML so that we can get the PHASE3 identifiers.
    survey_form = eso._request(
        "GET", "http://archive.eso.org/wdb/wdb/adp/phase3_main/form", 
        cache=False)
    query_dict = dict(wdbo="html/display", phase3_program="HARPS",
        max_rows_returned=eso.ROW_LIMIT, target=target)

    survey_response = eso._activate_form(survey_form, form_index=0, 
        inputs=query_dict, cache=False)

    r_html = survey_response.content

    # Find all the PHASE3 data sets. These identifiers will give us back the
    # intermediate pipeline products, not just the final spectrum.
    
    datasets = []
    for arcfile in eso_catalog["ARCFILE"]:
        match = re.findall(
            "PHASE3\+[0-9]+\+{}".format(arcfile.replace(".", "\.")), r_html)
        assert len(match) == 1
        datasets.extend(match)

    eso_catalog["dataset"] = datasets
    
    K = len(eso_catalog)
    M += K

    print("Found {} datasets ({} total so far) for {} ({}/{})".format(
        K, M, target, i, N))
    
    # Save the catalog.
    eso_catalog.write(eso_catalog_path, overwrite=True)
    print("Saved dataset records to {}".format(eso_catalog_path))
