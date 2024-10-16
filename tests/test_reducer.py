import pytest
import os
import astropy
import numpy as np

from pathlib import Path
from libraries.reducer import Reducer

def test_align_frames(obj_dir):
   assert isinstance(Reducer.align_frames(obj_dir), list)
   aligned_frames = os.listdir(obj_dir / 'align')
   print(f"aligned frames dir: {aligned_frames}")
   assert len(os.listdir(obj_dir / 'align')) > 0

def test_align_frames_bad_path():
    bad_path = 'qweruiopasdfjkl;'
    with pytest.raises(FileNotFoundError):
        Reducer.align_frames(bad_path)

def test_align_frames_frames_exist(capsys, obj_dir):
    Reducer.align_frames(obj_dir)
    Reducer.align_frames(obj_dir)
    out, err = capsys.readouterr()
    skip_statement = "Skipping alignment on frame"
    assert skip_statement.lower() in out.lower()

def test_make_dark(dark_dir, helper):
    assert isinstance(Reducer.make_dark(dark_dir), astropy.nddata.CCDData)

def test_make_dark_exist(capsys, dark_dir):
    Reducer.make_dark(dark_dir)
    Reducer.make_dark(dark_dir)
    out, err = capsys.readouterr()
    extant_dark_statement = "Reading extant master dark"
    assert extant_dark_statement.lower() in out.lower()

def test_make_flat(flat_dir, dark_dir, helper):
    assert isinstance(Reducer.make_flat(flat_dir, dark_dir), astropy.nddata.CCDData)

def test_make_flat_exist(capsys, flat_dir, dark_dir):
    Reducer.make_flat(flat_dir, dark_dir)
    Reducer.make_flat(flat_dir, dark_dir)
    out, err = capsys.readouterr()
    extant_flat_statement = "Reading extant flat"
    assert extant_flat_statement.lower() in out.lower()

def test_make_mask(object_frame):
    mask, boxes = Reducer.make_mask(object_frame)
    assert isinstance(mask, np.ndarray)
    assert isinstance(boxes, dict)

def test_make_stack(obj_dir):
    assert isinstance(Reducer.make_stack(obj_dir), astropy.nddata.CCDData)

def test_make_stack_exist(capsys, obj_dir):
    Reducer.make_stack(obj_dir)
    Reducer.make_stack(obj_dir)
    out, err = capsys.readouterr()
    extant_stack_statement = "Reading extant stack"
    assert extant_stack_statement.lower() in out.lower()

def test_reduce_objects(obj_dir, flat_dir, dark_dir, bkg_method='flat'):
    assert isinstance(Reducer.reduce_objects(obj_dir, flat_dir, dark_dir,  bkg_method), list)

def test_reduce_objects_skip(capsys, obj_dir, flat_dir, dark_dir, bkg_method='flat'):
    Reducer.reduce_objects(obj_dir, flat_dir, dark_dir, bkg_method)
    Reducer.reduce_objects(obj_dir, flat_dir, dark_dir, bkg_method)
    out, err = capsys.readouterr()
    skip_statement = "Skipping reduction on frame"
    assert skip_statement.lower() in out.lower()

def test_solve_plates(capsys, obj_dir):
    assert isinstance(Reducer.solve_plates(obj_dir), list)
    out, err = capsys.readouterr()
    not_solved_statement = "did not solve"
    assert not_solved_statement.lower() in out.lower()

#def test_solve_plates_skip(capsys, helper, obj_dir):
#    Reducer.solve_plates(obj_dir)
#    Reducer.solve_plates(obj_dir)
#    out, err = capsys.readouterr()
#    skip_statement = "Skipping plate solve on frame"
#    assert skip_statement.lower() in out.lower()

def test_scandir(obj_dir, helper):
    helper.show_files(obj_dir)
