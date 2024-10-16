import pytest
import os
from pathlib import Path

import numpy as np
from astropy.nddata import CCDData
from photutils.datasets import make_4gaussians_image
from astropy.io import fits
import ccdproc

MAINDIR = os.getcwd()

class Helper:
    def show_files(path):
        print(f"Scanning Directory: {path}")
        with os.scandir(path) as iterable_items:
            print("Name\t\tType")
            for entry in iterable_items:
                if entry.is_file():
                    entry_type = "file"
                    print(f"{entry.name}\t{entry_type}")
                
                if entry.is_dir():
                    entry_type = "dir"
                    print(f"{entry.name}\t\t{entry_type}")
                    with os.scandir(path / entry.name ) as iterable_dir:
                        for subfile in iterable_dir:
                            print(f" |-> {subfile.name}")
    
@pytest.fixture(scope="session")
def helper():
    return Helper

@pytest.fixture(scope="session")
def obj_dir(tmp_path_factory):
    # Make directories for fit files
    #  needs data, wcs, align, flat, dark, raw, cal

    fit_files_path = tmp_path_factory.mktemp("data")
    
    wcs_path = fit_files_path / "wcs"
    wcs_path.mkdir()
    
    align_path = fit_files_path / "align"
    align_path.mkdir()

    flat_path = fit_files_path / "flat"
    flat_path.mkdir()

    dark_path = fit_files_path / "dark"
    dark_path.mkdir()
    
    raw_path = fit_files_path / "raw"
    raw_path.mkdir()

    cal_path = fit_files_path / "cal"
    cal_path.mkdir()

    # Make fit file data
    hdu = example_hdu()
    flat_data = example_flat_data(data=hdu.data.mean(dtype="float"), shape=hdu.data.shape)
    flat = example_hdu(data=flat_data)
    dark_data = example_flat_data(data=1.0, shape=hdu.data.shape)
    dark = example_hdu(data=dark_data)

    print(f"Making file path: {fit_files_path}")
    _ = [hdu.writeto(fit_files_path / f"ccd-{i}.fit") for i in range(3)]
    _ = [hdu.writeto(raw_path / f"raw-{i}.fit") for i in range(3)]
    _ = [flat.writeto(flat_path / f"flat-{i}.fit") for i in range(3)]
    _ = [hdu.writeto(wcs_path / f"wcs-{i}.fit") for i in range(3)]
    _ = [dark.writeto(dark_path/ f"dark-{i}.fit") for i in range(3)]

    Helper.show_files(fit_files_path)

    return fit_files_path

@pytest.fixture(scope="session")
def flat_dir(obj_dir):
    return obj_dir / 'flat'

@pytest.fixture(scope="session")
def dark_dir(obj_dir):
    return obj_dir / 'dark'

@pytest.fixture(scope="session")
def object_frame(obj_dir):
    object_frame_path = obj_dir / 'ccd-0.fit'
    return object_frame_path

def example_ccd_data():
    #return CCDData(np.random.rand(2,2), unit="adu")
    return make_4gaussians_image()

def example_flat_data(data = 1.0, shape=(100,200)):
    data = np.ones(shape)*data
    return CCDData(data, unit="adu")

def example_header_with_wcs():
    path = Path(MAINDIR)
    path = path / 'tests'/ 'example_header.txt'
    with open(path, 'r') as file:
        data = file.read()
    header = fits.Header()
    header = header.fromstring(data)
    return header

def example_hdu(data = example_ccd_data()):
    hdu = fits.PrimaryHDU(data, example_header_with_wcs())
    return hdu
