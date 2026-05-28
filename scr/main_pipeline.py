from config import Configuration
from lib.photometry import Photometry
from lib.query import Query
from lib.survey import Survey
import glob
import numpy as np
import os
import time

from astropy import log
log.setLevel('WARNING')

os.system('clear')
#print('[\033[1m' + 'Πλούτων' + '\033[0m]')
print('[\033[1m' + 'ᾍδης ζῇ' + '\033[0m]')
print('\nRunning main pipeline\n')
st = time.time()

# configure input directories
input_dir = '/media/epimetheus/692e5e1e-6b16-4928-a25f-46dbabe207e6/toros_data/'
dark_dir = input_dir + 'dark/'
flat_dir = input_dir + 'flat/'
raw_dir = input_dir + 'raw/'

# configure output directories
output_directory = '/media/epimetheus/ExtremeSSD/'
dark_directory = output_directory + 'dark/'
flat_directory = output_directory + 'flat/'
clean_frame_directory = output_directory + 'clean/'
output_dir_list = [dark_directory, flat_directory, clean_frame_directory]

for dr in output_dir_list:
	if not os.path.exists(dr):
		os.mkdir(dr)

# grab the dates
field = Configuration.FIELD
field_id = field.split('_')[1]
ra = float(Survey.get_field(field_id)[0])
dec = float(Survey.get_field(field_id)[1])

date_list = Configuration.DATES
num_dates = len(date_list)

# run the main pipeline
for dte in range(num_dates):
	date = date_list[dte]

	output_name = 'stk_' + str(field) + '_' + str(date)
	output_directory = Configuration.OUTPUT_DATA_DIRECTORY + 'clean/' + date + '/' + field + '/'

	query_table_path = output_directory + 'table_query' + Configuration.TABLE_EXTENSION
	source_table_path = output_directory + 'table_source' + Configuration.TABLE_EXTENSION
	match_table_path = output_directory + 'table_match' + Configuration.TABLE_EXTENSION
	master_table_path = output_directory + 'table_master' + Configuration.TABLE_EXTENSION

	# make a master (flat) dark
	Photometry.make_dark(date, 'flat')

	# make a master (light) dark
	Photometry.make_dark(date, 'light')

	# make a normalized flatfield
	Photometry.make_flat(date)

	# reduce the frames
	frame_table = Photometry.clean_raw_frames(date, field)

	# make a stack
	stack_data, stack_header = Photometry.make_stack(date, field, frame_table)

	# create a queried source catalog
	query_table = Query.gaia_cone(ra, dec, Configuration.QUERY_RADIUS, query_table_path)

	# create an extracted source catalog
	source_table = Photometry.extract_sources(stack_data, stack_header, source_table_path)

	# match the catalogs
	match_table = Photometry.match_catalogs(source_table, query_table, match_table_path)

	# perform photometry on the stack
	master_table = Photometry.frame_aperture_photometry(date, field, stack_data, stack_header, match_table, master_table_path, output_name=output_name)

	# perform timeseries on clean frames
	test_ra = 148.2712725
	test_dec = -6.5137330
	timeseries_table = Photometry.timeseries(field, date, test_ra, test_dec)

	#jd, airmass, flux, flux_error, inst_mag, inst_mag_error = Photometry.point_aperture_photometry(date, field, stack_data, stack_header, test_ra, test_dec)

	#print(jd, airmass, flux, flux_error, inst_mag, inst_mag_error)

	# difference frames
	Photometry.difference_frames(field, date)

fn = time.time()
dt = np.around(fn - st, decimals=2)
print('\nMain pipeline finished in', dt, 's')