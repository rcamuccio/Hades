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

import logging
logger = logging.getLogger('astroquery')
logger.setLevel(logging.INFO)

os.system('clear')
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

# grab the dates
field = Configuration.FIELD
ra = float(Survey.get_field(field)[0])
dec = float(Survey.get_field(field)[1])
date_list = Survey.get_field_information(field)
num_dates = len(date_list)

# run the main pipeline
for dte in range(num_dates):
	date = date_list[dte]
	out_name = 'FIELD_' + str(field) + '_' + str(date)

	# define the output directories
	date_dark_directory = dark_directory + date + '/FIELD_' + field + '/'
	date_flat_directory = flat_directory + date + '/FIELD_' + field + '/'
	date_clean_frame_directory = clean_frame_directory + date + '/FIELD_' + field + '/'
	bkg_directory = date_clean_frame_directory + 'bkg/'
	hst_directory = date_clean_frame_directory + 'hst/'
	img_directory = date_clean_frame_directory + 'img/'
	pho_directory = date_clean_frame_directory + 'pho/'
	plt_directory = date_clean_frame_directory + 'plt/'
	res_directory = date_clean_frame_directory + 'res/'
	tbl_directory = date_clean_frame_directory + 'tbl/'

	# create the output directories
	directory_list = [date_dark_directory, date_flat_directory, date_clean_frame_directory, bkg_directory, hst_directory, img_directory, pho_directory, plt_directory, res_directory, tbl_directory]
	for dr in directory_list:
		if not os.path.exists(dr):
			os.makedirs(dr)

	# define the output tables
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

	# select stars for photometry
	filtered_table = Photometry.select_stars(master_table)

	# perform photometry on the clean frames
	Photometry.frame_timeseries(date, field, match_table)

	# perform timeseries on clean frames
	#timeseries_table = Photometry.timeseries(field, date, Configuration.SOURCE_RA, Configuration.SOURCE_DEC)

	# difference frames
	#Photometry.difference_frames(field, date)

# create global


fn = time.time()
dt = np.around(fn - st, decimals=2)
print('\n[\033[1m' + 'ᾍδης ἀπέρχεται' + '\033[0m] - Main pipeline finished (' + str(dt) + 's)\n')