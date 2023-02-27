import h5py as h5
from glob import glob
from tqdm import tqdm
import os
from astropy.io import fits
from astropy.table import Table

import untangle



data_dir = "51Peg/"

bis_header_keys = (
    "ESO DRS DVRMS", "ESO DRS CCF RVC", "ESO DRS BERV", "ESO DRS BJD", "ESO DRS BERVMX", "EXPTIME",
    "ESO DRS BIS RV", "ESO DRS BIS SPAN", "ESO DRS CCF CONTRAST", "ESO DRS CCF FWHM",
    "ESO DRS CCF LINES",
    "ESO DRS CCF MASK", 
    "ESO DRS CCF NOISE",
    "RA", "DEC", "MJD-OBS", "DATE-OBS", "UTC", "LST",   
)

compression = "gzip"
include_associated_categories = ("THAR_THAR", "FLAT", "THAR_FP")
require_associated_categories = ("THAR_THAR", "FLAT")

slices = [
    (0, 256),
    (240, 496),
    (480, 736),
    (720, 976),
    (960, 1216),
    (1200, 1456),
    (1440, 1696),
    (1680, 1936),
    (1920, 2176),
    (2160, 2416),
    (2400, 2656),
    (2640, 2896),
    (2880, 3136),
    (3120, 3376),
    (3360, 3616),
    (3600, 3856),
    (3840, 4096),
]


fits_paths = glob(f"{data_dir}*.fits.Z")
xml_paths = glob(f"{data_dir}*.xml")


def store_fits_image_as_h5_group(parent_group, name, x_slice=None, y_slice=None):
    group = parent_group.create_group(name)

    image_name = os.path.join(data_dir, f"{name}.fits")

    assert os.path.exists(image_name)
    
    with fits.open(image_name) as image:
        assert len(image) == 3

        for j, hdu in enumerate(image):
            if hdu.data is None:
                data_shape = (0, )
            else:
                if x_slice is not None or y_slice is not None:
                    data_shape = hdu.data[x_slice, y_slice].shape
                else:
                    data_shape = hdu.data.shape

            dataset = group.create_dataset(
                f"hdu_{j}", 
                data_shape, 
                dtype='i', 
                compression=compression,
                track_order=True
            )

            if hdu.data is not None:
                if x_slice is not None or y_slice is not None:
                    dataset[:] = hdu.data[x_slice, y_slice]
                else:
                    dataset[:] = hdu.data

            for key, value in hdu.header.items():
                dataset.attrs[key] = value
    return group



for xs, xe in slices:
    for ys, ye in slices:
        output_path = f"51Peg_{xs}-{xe}_{ys}-{ye}.h5"

        if os.path.exists(output_path):
            print(f"Skipping because {output_path} exists")
            continue

        print(f"Doing {output_path}")

        x_slice = slice(xs, xe)
        y_slice = slice(ys, ye)

        output = h5.File(output_path, "w")

        visits = output.create_group("visits")
        images = output.create_group("images")

        for i, xml_path in enumerate(tqdm(xml_paths)):
            struct = untangle.parse(xml_path)

            main_file = struct.association.mainFiles.file
            name = main_file["name"]
            category = main_file["category"]

            print(f"Main file: {name} ({category})")

            try:
                associations = struct.association.associatedFiles.association
            except AttributeError:
                print(f"\tNo associations. Continuing")
                continue
            else:
                calib_categories = [ea["category"] for ea in associations]
                if set(require_associated_categories).difference(calib_categories):
                    print(f"\tMissing categories")
                    continue
                else:
                    print(f"\tHas all categories")
                    # It has the calibrations we require
                    # Get the associated radial velocity first.
                    bis_paths = glob(f"{data_dir}/data/reduced/*/{name}_bis*.fits")
                    if len(bis_paths) == 0:
                        # No radial velocity.
                        print(f"\tNo radial velocity file from the DRS")
                        continue

                    elif len(bis_paths) > 1:
                        print(f"\tMultiple BIS files for {name}: {bis_paths}")
                        print(f"\tTaking the first one")

                    bis_path = bis_paths[0]
                    bis_values = {}
                    with fits.open(bis_path) as image:
                        for key in bis_header_keys:
                            bis_values[key] = image[0].header[f"HIERARCH {key}"]

                    image_group = store_fits_image_as_h5_group(images, name, x_slice=x_slice, y_slice=y_slice)
                    visit = visits.create_group(name)

                    # Add BIS headers 
                    for key, value in bis_values.items():
                        visit.attrs[key] = value

                    for association in associations:
                        if association["category"] not in include_associated_categories:
                            continue

                        for associated_main_file in association.mainFiles.file:
                            associated_name = associated_main_file['name']
                            associated_category = associated_main_file['category']

                            if associated_name not in images:
                                store_fits_image_as_h5_group(images, associated_name, x_slice=x_slice, y_slice=y_slice)

                            visit.attrs[associated_name] = associated_category

        output.close()
        print(f"Created {output_path}")