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
input_dark_directory = input_data_directory + 'dark/'
input_flat_directory = input_data_directory + 'flat/'
raw_frame_directory = input_data_directory + 'raw/'

# configure output directories
output_data_directory = Configuration.OUTPUT_DATA_DIRECTORY
output_dark_directory = output_data_directory + 'dark/'
output_flat_directory = output_data_directory + 'flat/'
clean_frame_directory = output_data_directory + 'clean/'
output_data_directory_list = [output_dark_directory, output_flat_directory, clean_frame_directory]

# grab the dates
field = Configuration.FIELD
ra = float(Survey.get_field(field)[0])
dec = float(Survey.get_field(field)[1])
date_list = Survey.get_field_information(field)
n_dates = len(date_list)

# run the main pipeline
for dte in range(n_dates):
	date = date_list[dte]
	date_extension = date + '/FIELD_' + field + '/'

	# define the output directories
	date_dark_directory = output_dark_directory + date_extension
	date_flat_directory = output_flat_directory + date_extension
	date_clean_frame_directory = clean_frame_directory + date_extension
	bkg_directory = date_clean_frame_directory + 'background/'
	#gal_directory = date_clean_frame_directory + 'gal/'
	flx_directory = date_clean_frame_directory + 'flux/'
	hst_directory = date_clean_frame_directory + 'histogram/'
	img_directory = date_clean_frame_directory + 'image/'
	#pho_directory = date_clean_frame_directory + 'pho/'
	pht_directory = date_clean_frame_directory + 'photometry/'
	#res_directory = date_clean_frame_directory + 'res/'
	tbl_directory = date_clean_frame_directory + 'table/'
	var_directory = date_clean_frame_directory + 'variable/'

	# create the output directories
	directory_list = [date_dark_directory, date_flat_directory, date_clean_frame_directory, bkg_directory, flx_directory, hst_directory, img_directory, pht_directory, tbl_directory, var_directory]
	for dr in directory_list:
		if not os.path.exists(dr):
			os.makedirs(dr)

	#
	# PREPROCESSING
	#

	# create the calibration frames
	Photometry.make_dark(date, 'flat')
	Photometry.make_dark(date, 'light')
	Photometry.make_flat(date)

	# reduce the light frames
	frame_table = Photometry.clean_raw_frames(date, field)

	# create a master stack frame
	stack_name = 'stack_FIELD_' + str(field) + '_' + str(date)
	stack_data, stack_header = Photometry.make_stack(date, field, frame_table)

	# extract sources
	tbl_source_path = tbl_directory + 'source' + Configuration.TABLE_EXTENSION
	source_table = Photometry.extract_sources(stack_data, stack_header, tbl_source_path)
	
	#
	# GAIA PHOTOMETRY
	# 

	# perform a source query
	tbl_query_gaia_path = tbl_directory + 'query_gaia' + Configuration.TABLE_EXTENSION
	query_gaia_table = Query.gaia_cone(ra, dec, Configuration.QUERY_RADIUS, tbl_query_gaia_path)

	# match the catalogs
	tbl_match_gaia_path = tbl_directory + 'match_gaia' + Configuration.TABLE_EXTENSION
	match_gaia_table = Photometry.match_gaia_catalog(source_table, query_gaia_table, tbl_match_gaia_path, 'gaia_dr3')

	# perform photometry on the stack
	tbl_master_path = tbl_directory + 'master_gaia' + Configuration.TABLE_EXTENSION
	master_table = Photometry.single_frame_gaia_aperture_photometry(date, field, stack_data, stack_header, match_gaia_table, tbl_master_path, stack_name)

	# perform photometry on the series
	Photometry.timeseries_gaia_aperture_photometry(date, field, match_gaia_table, 'gaia_dr3')

	#
	# VARIABLE PHOTOMETRY
	#

	# perform a source query
	tbl_query_aavso_path = tbl_directory + 'query_aavso' + Configuration.TABLE_EXTENSION
	query_aavso_table = Query.aavso_vsx_cone(ra, dec, Configuration.QUERY_RADIUS, tbl_query_aavso_path)

	# match the catalogs
	tbl_match_aavso_path = tbl_directory + 'match_aavso' + Configuration.TABLE_EXTENSION
	match_aavso_table = Photometry.match_aavso_catalog(source_table, query_aavso_table, tbl_match_aavso_path)

	# perform photometry on the series
	Photometry.timeseries_aavso_aperture_photometry(date, field, match_aavso_table, 'aavso_vsx')

	# 
	# GALAXY PHOTOMETRY
	# 

	# perform a galaxy query
	#tbl_query_glade_path = tbl_directory + 'query_glade' + Configuration.TABLE_EXTENSION
	#query_glade_table = Query.glade_plus_cone(ra, dec, Configuration.QUERY_RADIUS, tbl_query_glade_path)

	# select stars for photometry
	#filtered_table = Photometry.select_stars(master_table)

	# perform timeseries on clean frames
	#timeseries_table = Photometry.timeseries(field, date, Configuration.SOURCE_RA, Configuration.SOURCE_DEC)

	# difference frames
	#Photometry.difference_frames(field, date)

# create global

fn = time.time()
dt = np.around(fn - st, decimals=2)
print('\n[\033[1m' + 'ᾍδης ἀπέρχεται' + '\033[0m] - Main pipeline finished (' + str(dt) + 's)\n')