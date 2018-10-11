
import re
import asyncio # so you know we are serious
import requests
import tempfile
import yaml

from astropy.table import Table
from astropy.extern.six import BytesIO
from collections import OrderedDict

from bs4 import BeautifulSoup


class Harps(object):

    row_limit = 100000

    def __init__(self, credentials_path):
        self._credentials_path = credentials_path

        with open(credentials_path, "r") as fp:
            credentials = yaml.load(fp)["eso"]

        self.login(**credentials)

        return None



    def login(self, username, password):

        session = requests.Session()
        prepare = session.get("https://www.eso.org/sso")

        root = BeautifulSoup(prepare.content, 'html5lib')
        login_input = root.find(name='input', attrs={'name': 'execution'})

        if login_input is None:
            raise ValueError("ESO page did not have the correct attribute")

        execution = login_input.get('value')

        login = session.post(prepare.url, 
                             data=dict(
                                username=username,
                                password=password,
                                execution=login_input.get("value"),
                                _eventId="submit",
                                geolocation=""
                            ))
        login.raise_for_status()

        root = BeautifulSoup(login.content, "html5lib")
        authenticated = (root.find("h4").text == "Login successful")

        assert authenticated

        self.session = session    



    def query_position(self, ra, dec, **kwargs):
        return self._query(coord1=("", f"{ra}"), coord2=("", f"{dec}"))

    def query_target(self, target_name, **kwargs):
        return self._query(target=("", f"{target_name}"))
        

    def _query(self, **params):
        payload = OrderedDict([
            ("wdbo", ("", "html/display")),
            ("max_rows_returned", ("", "{:.0f}".format(self.row_limit))),
            ("target", ("", "")),
            ("resolver", ("", "simbad")),
            ("wdb_input_file", ("", "", "application/octet-stream")),
            ("coord_sys", ("", "eq")),
            ("coord1", ("", "")),
            ("coord2", ("", "")),
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
        ])
        payload.update(params)


        response = self.session.post(
            "http://archive.eso.org/wdb/wdb/adp/phase3_main/query", 
            files=payload)

        content = str(response.content)

        rows = "\n".join([r for r in content.split("\n") if "PHASE3+" in r])

        names = ("Mark", "More", "ARCFILE", "HDR", "Object", "RA", "DEC", "Filter",
                "ABMAGLIM", "Wavelength", "SNR", "Resolution", "Product category",
                "Instrument", "Date Obs", "Exptime", "Collection", "Product version",
                "Release Description", "Run/Program ID", "ORIGFILE", "REFERENCE Catalog",
                "Interface")
        delete_names = ("Mark", "More", "HDR", "Filter", "ABMAGLIM",
            "Product category", "Collection", "Product version", "Release Description",
            "Interface", "REFERENCE Catalog")

        _ = "<TR"
        rows = _ + _.join(rows.replace("[doc&nbsp;id:", "[doc:").split(_)[1:])
        if rows == _:
            return Table(names=[n for n in names if n not in delete_names])

        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(bytes(f"<table>{rows}</table>", encoding="utf-8"))
        f.close()


        table = Table.read(f.name, format="ascii.html", names=names)

        # Delete unnecessary columns.
        for column_name in delete_names:
            del table[column_name]

        # Parse the PHASE3 identifiers.
        table["phsae3_dataset_identifier"] = re.findall(
            "PHASE3\+[0-9]+\+ADP\.[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}\.[0-9]{3}", 
            content)

        return table




