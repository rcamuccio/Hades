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
	tbl_query_aavso_path = tbl_directory + 'query_aavso' + Configuration.TABLE_EXTENSION
	tbl_query_gaia_path = tbl_directory + 'query_gaia' + Configuration.TABLE_EXTENSION
	tbl_query_glade_path = tbl_directory + 'query_glade' + Configuration.TABLE_EXTENSION

	tbl_source_path = tbl_directory + 'source' + Configuration.TABLE_EXTENSION
	tbl_match_path = tbl_directory + 'match' + Configuration.TABLE_EXTENSION
	tbl_master_path = tbl_directory + 'master' + Configuration.TABLE_EXTENSION

	# create the master dark frames
	Photometry.make_dark(date, 'flat')
	Photometry.make_dark(date, 'light')

	# create the flatfield frame
	Photometry.make_flat(date)

	# reduce the light frames
	frame_table = Photometry.clean_raw_frames(date, field)

	# create the master stack frame
	stack_name = 'stack_FIELD_' + str(field) + '_' + str(date)
	stack_data, stack_header = Photometry.make_stack(date, field, frame_table)

	# create a queried source catalog
	query_table = Query.gaia_cone(ra, dec, Configuration.QUERY_RADIUS, tbl_query_gaia_path)

	# create an extracted source catalog
	source_table = Photometry.extract_sources(stack_data, stack_header, tbl_source_path)

	# match the catalogs
	match_table = Photometry.match_catalogs(source_table, query_table, tbl_match_path)

	# perform photometry on the stack
	master_table = Photometry.frame_aperture_photometry(date, field, stack_data, stack_header, match_table, tbl_master_path, output_name=stack_name)

	# select stars for photometry
	#filtered_table = Photometry.select_stars(master_table)

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