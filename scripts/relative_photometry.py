"""
relative_photometry

"""

from lib.calculator import *
from lib.photometer import *
from lib.querier import *
from lib.reducer import *

import configparser
import os
import time

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter("ignore", category=AstropyWarning)

os.system("clear")
print("Running script [relative_photometry]")
print()
start_time = time.time()

# --- Create class instances
calculator = Calculator()
photometer = Photometer()
querier = Querier()
reducer = Reducer()

# --- Read configuration file
config = configparser.ConfigParser()
config.read("config.ini")

main_dir = config["Directory"]["main_dir"]

camera = config["Observatory"]["camera"]
latitude = config["Observatory"]["latitude"]
longitude = config["Observatory"]["longitude"]
height = config["Observatory"]["height"]

airmass_method = config["Photometry"]["airmass_method"]
catalog = config["Photometry"]["catalog"]
catalog_f1 = config["Photometry"]["catalog_f1"]
catalog_f2 = config["Photometry"]["catalog_f2"]
field_ra = config["Photometry"]["field_ra"]
field_dec = config["Photometry"]["field_dec"]
phot_f1 = config["Photometry"]["phot_f1"]
phot_f2 = config["Photometry"]["phot_f2"]
radius_an_in = config["Photometry"]["radius_an_in"]
radius_an_out = config["Photometry"]["radius_an_out"]
radius_aper = config["Photometry"]["radius_aper"]
radius_solve = config["Photometry"]["radius_solve"]
radius_query = config["Photometry"]["radius_query"]

combine_method = config["Preprocessing"]["combine_method"]
dilate_size = config["Preprocessing"]["dilate_size"]
dtype = config["Preprocessing"]["dtype"]
mem_limit = config["Preprocessing"]["mem_limit"]
npixels = config["Preprocessing"]["npixels"]
sigma = config["Preprocessing"]["sigma"]

# --- Build directory tree
dark_flat_f1_dir = main_dir + "/dark/flat-" + phot_f1
dark_flat_f2_dir = main_dir + "/dark/flat-" + phot_f2
dark_obj_dir = main_dir + "/dark/obj"

flat_f1_dir = main_dir + "/flat/" + phot_f1
flat_f2_dir = main_dir + "/flat/" + phot_f2

obj_early_f1_dir = main_dir + "/obj/early/" + phot_f1
obj_early_f2_dir = main_dir + "/obj/early/" + phot_f2
obj_main_dir = main_dir + "/obj/main"
obj_late_f1_dir = main_dir + "/obj/late/" + phot_f1
obj_late_f2_dir = main_dir + "/obj/late/" + phot_f2

# --- Make master darks
reducer.make_dark(dark_flat_f1_dir)
reducer.make_dark(dark_flat_f2_dir)
reducer.make_dark(dark_obj_dir)

# --- Make flatfields
reducer.make_flat(flat_f1_dir, dark_flat_f1_dir)
reducer.make_flat(flat_f2_dir, dark_flat_f2_dir)

# --- Reduce object frames
reducer.reduce_objects(obj_early_f1_dir, flat_f1_dir, dark_obj_dir)
reducer.reduce_objects(obj_early_f2_dir, flat_f2_dir, dark_obj_dir)
reducer.reduce_objects(obj_main_dir, flat_f1_dir, dark_obj_dir)
reducer.reduce_objects(obj_late_f1_dir, flat_f1_dir, dark_obj_dir)
reducer.reduce_objects(obj_late_f2_dir, flat_f2_dir, dark_obj_dir)

# --- Plate solve object frames
reducer.solve_plate(obj_early_f1_dir, field_ra, field_dec, radius_solve)
reducer.solve_plate(obj_early_f2_dir, field_ra, field_dec, radius_solve)
reducer.solve_plate(obj_main_dir, field_ra, field_dec, radius_solve)
reducer.solve_plate(obj_late_f1_dir, field_ra, field_dec, radius_solve)
reducer.solve_plate(obj_late_f2_dir, field_ra, field_dec, radius_solve)

# --- Align object frames
obj_early_f1_list = reducer.align_frames(obj_early_f1_dir)
obj_early_f2_list = reducer.align_frames(obj_early_f2_dir)
obj_main_list = reducer.align_frames(obj_main_dir)
obj_late_f1_list = reducer.align_frames(obj_late_f1_dir)
obj_late_f2_list = reducer.align_frames(obj_late_f2_dir)

# --- Stack object frames
reducer.make_stack(obj_early_f1_dir)
reducer.make_stack(obj_early_f2_dir)
reducer.make_stack(obj_main_dir)
reducer.make_stack(obj_late_f1_dir)
reducer.make_stack(obj_late_f2_dir)

# --- Submit PS1 query
os.chdir(obj_main_dir)

query_table = querier.query_ps1_cone(field_ra, field_dec, radius_query)
query_table.write("table-query.cat", format="ascii.fixed_width", overwrite=True)

stack_frame = "stack.fit"

source_table = photometer.extract_sky_sources(stack_frame)
source_table.write("table-source.cat", format="ascii.fixed_width", overwrite=True)

match_table = photometer.match_catalogs(source_table, query_table)
match_table.write("table-match.cat", format="ascii.fixed_width", overwrite=True)

master_table = photometer.photometry_field(stack_frame, match_table, catalog_f1, catalog_f2, phot_f1, phot_f2)
master_table.write("table-master.cat", format="ascii.fixed_width", overwrite=True)

os.chdir(obj_early_f1_dir)
obj_early_f1_table = photometer.extract_sky_sources(stack_frame)

os.chdir(obj_early_f2_dir)
obj_early_f2_table = photometer.extract_sky_sources(stack_frame)

os.chdir(obj_late_f1_dir)
obj_late_f1_table = photometer.extract_sky_sources(stack_frame)

os.chdir(obj_late_f2_dir)
obj_late_f2_table = photometer.extract_sky_sources(stack_frame)

end_time = time.time()
total_time = end_time - start_time
print()
print("Script ended in", "%.1f" % total_time, "seconds")