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
print('[\033[1m' + 'ᾍδης ζῇ' + '\033[0m] - Running main pipeline\n')
st = time.time()

# configure input directories
input_data_directory = Configuration.INPUT_DATA_DIRECTORY
dark_directory = input_data_directory + 'dark/'
flat_directory = input_data_directory + 'flat/'
raw_frame_directory = input_data_directory + 'raw/'

# configure output directories
output_data_directory = Configuration.OUTPUT_DATA_DIRECTORY
dark_directory = output_data_directory + 'dark/'
flat_directory = output_data_directory + 'flat/'
clean_frame_directory = output_data_directory + 'clean/'
output_data_directory_list = [dark_directory, flat_directory, clean_frame_directory]
for dr in output_data_directory_list:
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

	out_name = str(field) + '_' + str(date)
	out_directory = clean_frame_directory + date + '/' + field + '/'

	bkg_directory = out_directory + 'bkg/'
	hst_directory = out_directory + 'hst/'
	img_directory = out_directory + 'img/'
	plt_directory = out_directory + 'plt/'
	res_directory = out_directory + 'res/'
	tbl_directory = out_directory + 'tbl/'
	out_directory_list = [bkg_directory, hst_directory, img_directory, plt_directory, res_directory, tbl_directory]
	for dr in out_directory_list:
		if not os.path.exists(dr):
			os.mkdir(dr)

	tbl_query_aavso_path = tbl_directory + 'tbl_query_aavso' + Configuration.TABLE_EXTENSION
	tbl_query_gaia_path = tbl_directory + 'tbl_query_gaia' + Configuration.TABLE_EXTENSION
	tbl_query_glade_path = tbl_directory + 'tbl_query_glade' + Configuration.TABLE_EXTENSION
	tbl_source_path = tbl_directory + 'tbl_source' + Configuration.TABLE_EXTENSION
	tbl_match_path = tbl_directory + 'tbl_match' + Configuration.TABLE_EXTENSION
	tbl_master_path = tbl_directory + 'tbl_master' + Configuration.TABLE_EXTENSION

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
	query_table = Query.gaia_cone(ra, dec, Configuration.QUERY_RADIUS, tbl_query_gaia_path)

	# create an extracted source catalog
	source_table = Photometry.extract_sources(stack_data, stack_header, tbl_source_path)

	# match the catalogs
	match_table = Photometry.match_catalogs(source_table, query_table, tbl_match_path)

	# perform photometry on the stack
	master_table = Photometry.frame_aperture_photometry(date, field, stack_data, stack_header, match_table, tbl_master_path, output_name=out_name)

	# perform timeseries on clean frames
	test_ra = 148.2712725
	test_dec = -6.5137330
	timeseries_table = Photometry.timeseries(field, date, test_ra, test_dec)

	# difference frames
	Photometry.difference_frames(field, date)

# create global



fn = time.time()
dt = np.around(fn - st, decimals=2)
print('\nMain pipeline finished in', dt, 's')