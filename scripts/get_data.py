
import os
import re
from astropy.table import Table
from glob import glob

from astroquery.eso import Eso as ESO 


files = glob("data/harps-spectra/eso-*.fits")[:10]


# Login to ESO.
eso = ESO()
eso.login("andycasey")

cwd = os.path.dirname(os.path.realpath(__file__))
data_dir = "{}/../data/harps-spectra/".format(cwd)
eso.cache_location = "{}/cache/".format(cwd)
# See https://github.com/astropy/astroquery/issues/800
eso._survey_list = ["HARPS"]
eso.ROW_LIMIT = 10000


data = []
for filename in files:
    t = Table.read(filename)
    data.extend([("dataset", dataset) for dataset in t["dataset"]])


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

confirmation_response = eso._session.request("POST", 
    "http://dataportal.eso.org/rh/requests/{}/submission".format(eso.USERNAME),
    data=data)
assert confirmation_response.ok

# Get the request number and get the download script from ESO.
_ = re.findall("Request #[0-9]+\w", confirmation_response.text)[0].split()[-1]
request_number = int(_.lstrip("#"))

remote_path = "https://dataportal.eso.org/rh/requests/{username}/{request}/script"\
    .format(username=eso.USERNAME, request=request_number)
local_path = os.path.join(data_dir, "eso-{}-{}.sh".format(request_number, target))

raise a


completed_response = eso._download_file(remote_path, local_path)
print("Download script from {} saved to {}".format(remote_path, local_path))


# Save the failures.
with open(os.path.join(cwd, "no_records-catalog.txt"), "w") as fp:
    fp.write("\n".join(no_records))

# Remove things from the cache.
for path in glob(eso.cache_location + "/*"):
    os.remove(path)
