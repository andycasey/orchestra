
"""
Download HARPS data from the ESO archive.
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
eso.cache_location = "{}/cache/".format(cwd)
# See https://github.com/astropy/astroquery/issues/800
eso._survey_list = ["HARPS"]
eso.ROW_LIMIT = 10000

eso_catalog = []
local_catalog = Table.read("data/HARPS_all.csv")
targets = set(local_catalog["Name"])
N = len(targets)
M = 0

for i, target in enumerate(targets):

    # Search the ESO archive for HARPS data.
    r = eso.query_survey("HARPS", target=target)
    if r is None:
        print("No results found for star name {}".format(target))
        continue

    # For some reason we get a lot of shit that we don't want.
    # Require:
    keep = (r["Instrument"] == "HARPS") \
         * (r["R (&lambda;/&delta;&lambda;)"] == 115000)
    r = r[keep]

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
    for arcfile in r["ARCFILE"]:

        match = re.findall(
            "PHASE3\+[0-9]+\+{}".format(arcfile.replace(".", "\.")), r_html)
        assert len(match) == 1

        datasets.extend(match)

    r["dataset"] = datasets

    # Stringify the Filter columns because these can differ in format.
    stringify_columns = ("Filter", "REFERENC")
    for column_name in stringify_columns:
        _ = map(str, r[column_name])
        del r[column_name]
        r[column_name] = _

    _ = len(r)
    M += _

    print("Found {} datasets ({} total so far) for {} ({}/{})".format(
        _, M, target, i, N))
    
    eso_catalog.append(r)


# Stack the records from ESO.
eso_catalog = vstack(eso_catalog)
N = len(eso_catalog)
print("{} datasets found (~{} MB)".format(N, N * 11))

eso_catalog_path = os.path.join(cwd, "eso-datasets.fits")
eso_catalog.write(eso_catalog_path, overwrite=True)
print("Saved dataset records to {}".format(eso_catalog_path))


# Request the data!
data = [("dataset", dataset) for dataset in eso_catalog["dataset"]]

prepare_response = eso._session.request("POST",
    "http://dataportal.eso.org/rh/confirmation", data=data)

# Additional payload items required for confirmation.
data += [
    ("requestDescription", ""),
    ("deliveryMediaType", "WEB"), # OR USB_DISK --> Holy shit what the fuck!
    ("requestCommand", "SELECTIVE_HOTFLY"),
    ("submit", "Submit")
]

confirmation_response = eso._session.request("POST", 
    "http://dataportal.eso.org/rh/requests/{}/submission".format(eso.USERNAME),
    data=data)

# Get the request number and get the download script from ESO.
_ = re.findall("Request #[0-9]+\w", confirmation_response.text)[0].split()[-1]
request_number = int(_.lstrip("#"))

remote_path = "https://dataportal.eso.org/rh/requests/{username}/{request}/script"\
    .format(username=eso.USERNAME, request=request_number)
local_path = os.path.join(cwd, "download_spectra.sh")

print("Giving ESO 10 seconds to prepare the download script at {}..."\
    .format(remote_path))

time.sleep(10)
completed_response = eso._download_file(remote_path, local_path)
print("Download script from {} saved to {}".format(remote_path, local_path))

# Remove things from the cache.
for path in glob(eso.cache_location + "/*"):
    os.remove(path)
