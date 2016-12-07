
"""
Retrieve reduced and intermediate HARPS data products from the ESO archive.


This script will download catalog files for each object which contain the 
Phase 3 identifier required to download the reduced and intermediate data
products.
"""

__author__ = "Andrew R. Casey <arc@ast.cam.ac.uk>"


import cPickle as pickle
import os
import psycopg2 as pg
import re
import time
import yaml
from astropy.table import Table
from glob import glob
from bs4 import BeautifulSoup

from astroquery.eso import Eso as ESO 

cwd = os.path.dirname(os.path.realpath(__file__))

SKIP = 0 # Skip how many batches at the start? (for if you are re-running this..)
BATCH = 2000 # How many datasets should we get per ESO request?
WAIT_TIME = 60 # Seconds between asking ESO if they have prepared our request
DATA_DIR = "{}/../data/spectra/".format(cwd) # Where's the spectra at?

# Connect to the PostgreSQL database.
with open(os.path.join(cwd, "../db/credentials.yaml"), "r") as fp:
    credentials = yaml.load(fp)
connection = pg.connect(**credentials)

remote_paths_path = os.path.join(cwd, "eso_retrieve_paths.pkl")

# Cross-match tables to find Phase 3 data products that we don't have yet.
cursor = connection.cursor()
cursor.execute(
    """ SELECT dataset FROM phase3_products
        WHERE NOT EXISTS(
            SELECT 1 FROM obs
            WHERE obs.date_obs = phase3_products.date_obs)
        ORDER BY phase3_products.date_obs ASC
    """)
records = cursor.fetchall()
cursor.close()
connection.close()

N = len(records)
I = N / BATCH + (1 if N % BATCH else 0)

assert not os.path.exists(remote_paths_path), "This will overwrite old paths!"

remote_paths = []

print("In total there are {} records we will request".format(N))

for i in range(I):

    if i < SKIP:
        print("Skipping {}".format(i + 1))
        continue

    print("Starting with batch number {}/{}".format(i + 1, I))

    data = [("dataset", dataset) for dataset in records[i*BATCH:(i + 1)*BATCH]]

    # Login to ESO.
    eso = ESO()
    eso.login("andycasey")

    prepare_response = eso._session.request("POST",
        "http://dataportal.eso.org/rh/confirmation", data=data)
    assert prepare_response.ok

    # Additional payload items required for confirmation.
    data += [
        ("requestDescription", ""),
        ("deliveryMediaType", "WEB"), # OR USB_DISK --> Holy shit what the fuck!
        ("requestCommand", "SELECTIVE_HOTFLY"),
        ("submit", "Submit")
    ]

    confirmation_response = eso._session.request(
        "POST", 
        "http://dataportal.eso.org/rh/requests/{}/submission".format(eso.USERNAME),
        data=data)
    assert confirmation_response.ok

    # Parse the request number so that we can get a download script from ESO later
    _ = re.findall("Request #[0-9]+\w", confirmation_response.text)[0].split()[-1]
    request_number = int(_.lstrip("#"))

    # Get the download scripts for our request number.
    url = "https://dataportal.eso.org/rh/requests/{}".format(eso.USERNAME)
    print("Retrieving remote paths for request number {}/{}: {}".format(
        i + 1, I, request_number))

    # Check if ESO is ready for us.
    while True:    

        # Login to ESO.
        eso = ESO()
        eso.login("andycasey")

        check_state = eso._request("GET", "{}/recentRequests".format(url))
        root = BeautifulSoup(check_state.text, "html5lib")

        link = root.find(href="/rh/requests/{}/{}".format(
            eso.USERNAME, request_number))
        image = link.find_next("img")
        state = image.attrs["alt"]

        print("Current state {} on request {} ({}/{})".format(
            state, request_number, i + 1, I))

        if state != "COMPLETE":

            # Remove anything from the astroquery cache.
            for cached_file in glob(os.path.join(eso.cache_location, "*")):
                os.remove(cached_file)

            print("Sleeping for {} seconds..".format(WAIT_TIME))
            time.sleep(WAIT_TIME)

        else:
            break

    response = eso._request(
        "GET", "{}/{}/script".format(url, request_number))
    
    paths = response.text.split("__EOF__")[-2].split("\n")[1:-2]
    print("Found {} remote paths for request_number {}".format(
        len(paths), request_number))
    remote_paths.extend(paths)

    # Save the remote paths in case of a catastrophic crash!
    with open(remote_paths_path, "wb") as fp:
        pickle.dump(remote_paths, fp, -1)

    # Remove anything from the astroquery cache.
    for cached_file in glob(os.path.join(eso.cache_location, "*")):
        os.remove(cached_file)
    
    # We have all the remote paths for this request. At ESO's advice, let's
    # wait another 60 seconds before starting our new request.
    time.sleep(60)

# Prepare the script for downloading.
template_path = os.path.join(cwd, "download_template.sh")
with open(template_path, "r") as fp:
    contents = fp.read()

script_path = os.path.join(DATA_DIR, "download_spectra.sh")
with open(script_path, "w") as fp:
    fp.write(contents.replace("$$REMOTE_PATHS$$", "\n".join(remote_paths)))

print("""
OK now do:
    cd "{}"
    sh download_spectra.sh

Then untar everything new with:
    grep tar download_spectra.sh | awk '{n=split($1,a,"/"); print "tar -xvf \"" a[n] " --keep-old-files --force-local "}' > untar.sh
    sh untar.sh

Then ingest everything by running:
    cd "{}"
    python scripts/correct_folder_structure.py
    python scripts/db_ingest_headers.py
""".format(
    os.path.realpath(DATA_DIR),
    os.path.realpath(os.path.join(cwd, ".."))))