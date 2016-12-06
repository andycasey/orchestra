
"""
Find HARPS data from the ESO archive for a set of target positions.

This script will download catalog files for each object which contain the 
Phase 3 identifier required to download the reduced and intermediate data
products.
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
import pg
import re
import time
from astropy.extern.six import BytesIO, cPickle as pickle
from astropy.table import Table
from astroquery.eso import Eso as ESO 
from astroquery.eso.core import _check_response


# Load local catalog of positions.
local_catalog = Table.read("data/HARPS_all.csv")

# Login to ESO.
eso = ESO()
eso.login("andycasey")
eso.ROW_LIMIT = 100000 # Maximum possible number of observations per star

cwd = os.path.dirname(os.path.realpath(__file__))
catalog_dir = "{}/../data/eso-phase3-products/".format(cwd)

# Connect to the PostgreSQL database.
with open("db/credentials.yaml", "r") as fp:
    credentials = yaml.load(fp)
connection = pg.connect(**credentials)

def query_harps_phase3_by_position(ra, dec, **kwargs):
    """
    Query the ESO Phase 3 science archive by position.

    :param ra:
        Right ascension [degrees].

    :param dec:
        Declination [degrees].
    """

    payload = [
        ("wdbo", ("", "html/display")),
        ("max_rows_returned", ("", "{:.0f}".format(eso.ROW_LIMIT))),
        ("target", ("", "")),
        ("resolver", ("", "simbad")),
        ("wdb_input_file", ("", "", "application/octet-stream")),
        ("coord_sys", ("", "eq")),
        ("coord1", ("", str(ra))),
        ("coord2", ("", str(dec))),
        ("box", ("", "02 09 00")),
        ("tab_ra", ("", "on")),
        ("tab_dec", ("", "on")),
        ("tab_filter", ("", "on")),
        ("filter", ("", "Any")),
        ("tab_wavelength", ("", "on")),
        ("wavelength", ("", "Any")),
        ("tab_dataproduct_type", ("", "on")),
        ("dataproduct_type", ("", "Any")),
        ("tel_id", ("", "Any")),
        ("tab_ins_id", ("", "on")),
        ("ins_id", ("", "HARPS")),
        ("obstech", ("", "Any")),
        ("tab_date_obs", ("", "on")),
        ("date_obs", ("", "")),
        ("mjd_obs", ("", "")),
        ("tab_exptime", ("", "on")),
        ("exptime", ("", "")),
        ("multi_ob", ("", "%")),
        ("tab_collection_name", ("", "on")),
        ("tab_prog_id", ("", "on")),
        ("prog_id", ("", "")),
        ("username", ("", "")),
        ("p3orig", ("", "%")),
        ("tab_origfile", ("", "on")),
        ("origfile", ("", "")),
        ("tab_dp_id", ("", "on")),
        ("dp_id", ("", "")),
        ("rel_date", ("", "")),
        ("tab_referenc", ("", "on")),
        ("referenc", ("", "")),
        ("batch_id", ("", "")),
        ("publication_date", ("", "")),
        ("wdb_input_file_raw", ("", "", "application/octet-stream")),
        ("order_main", ("", "dummy"))
    ]
    
    url = "http://archive.eso.org/wdb/wdb/adp/phase3_main/query"

    survey_response = eso._request("POST", url, cache=False, files=payload)
    content = survey_response.content

    if not _check_response(content):
        return None

    rows = "\n".join([r for r in content.split("\n") if "PHASE3+" in r])
    rows = rows.replace("[doc&nbsp;id:", "[doc:")
    html_content = "<table>{}</table>".format(rows)
    table = Table.read(BytesIO(html_content), format="ascii.html",
        names=("Mark", "More", "ARCFILE", "HDR", "Object", "RA", "DEC", "Filter",
            "ABMAGLIM", "Wavelength", "SNR", "Resolution", "Product category",
            "Instrument", "Date Obs", "Exptime", "Collection", "Product version",
            "Release Description", "Run/Program ID", "ORIGFILE", "REFERENCE Catalog",
            "Interface"))

    # Delete unnecessary columns.
    for column_name in ("Mark", "More", "HDR", "Filter", "ABMAGLIM",
        "Product category", "Collection", "Product version", "Release Description",
        "Interface", "REFERENCE Catalog"):
        del table[column_name]

    # Parse the PHASE3 identifiers.
    table["dataset"] = re.findall(
        "PHASE3\+[0-9]+\+ADP\.[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}\.[0-9]{3}", 
        content)
    return table


warnings = {}
failures = {}
M, N = (0, len(local_catalog))

for i, target in enumerate(local_catalog):

    eso_catalog_path = os.path.join(
        catalog_dir, "eso-{}.fits".format(target["Name"]))

    if os.path.exists(eso_catalog_path):
        print("Skipping {} because {} exists".format(
            target["Name"], eso_catalog_path))

        t = Table.read(eso_catalog_path)
        K = len(t)
        M += K

        if target["N_exp"] > K:
            warnings[target["Name"]] = "Expected {}; found {}".format(
                target["N_exp"], K)
            print("Warning: Expected {} and found {}".format(target["N_exp"], K))
        continue

    # Search the ESO archive for HARPS data.
    try:
        response = query_harps_phase3_by_position(target["RA"], target["Dec"])
    
    except ValueError:
        failures[target["Name"]] = "ValueError: only one result?"
        print("ValueError: Only one result for {}?".format(target["Name"]))
        continue

    if response is None:
        print("No results found for star name {}".format(target["Name"]))
        failures[target["Name"]] = "No results found in Phase 3 search"
        continue

    # Let's be spoilt as fuck
    keep = response["Resolution"] == 115000
    response = response[keep]

    if len(response) == 0:
        print("No R ~ 115,000 spectra found for star {}".format(target["Name"]))
        failures[target["Name"]] = "Only R ~ 80,000 spectra found"
        continue

    K = len(response)
    M += K

    print("({}/{}) Found {} datasets ({} expected; {} total so far) for {} ({:.3f} / {:.3f})"\
        .format(i, N, K, target["N_exp"], M, target["Name"], target["RA"], target["Dec"]))

    if target["N_exp"] > K:
        warnings[target["Name"]] = "Expected {}; found {}".format(
            target["N_exp"], K)
        print("Warning: Expected {} and found {}".format(target["N_exp"], K))

    # Ingest the catalog.
    print("Ingesting {} Phase 3 records from {}".format(len(response), filename))
    
    for record in response:
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

# Close the PostgreSQL database.
connection.close()

# Save any warnings or failures.
with open(os.path.join(cwd, "eso-search-phase3-output.pkl"), "wb") as fp:
    pickle.dump((warnings, failures), fp, -1)