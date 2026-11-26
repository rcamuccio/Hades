import glob
import numpy as np
import os
import time
from astropy.utils.exceptions import AstropyUserWarning
from astropy.io import ascii, fits
from astropy.table import Table
from astropy.time import Time

from conf import Configuration
from lib.boethos import Boethos
from lib.metretes import Metretes
from lib.skopeus import Skopeus

from astropy import log
log.setLevel('WARNING')
import logging
logger = logging.getLogger('astroquery')
logger.setLevel(logging.INFO)
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.simplefilter('ignore', AstropyUserWarning)

# write the preamble
os.system('clear')
script_preamble = Boethos.get_script_preamble('preprocessing_pipeline')
print(script_preamble)
line_preamble = Boethos.get_line_preamble('script', 'preprocessing_pipeline')

# start the clock
st = time.time()

# set the output directories
output_directory = Configuration.OUTPUT_DIRECTORY
output_dark_directory = output_directory + 'dark/'
output_flat_directory = output_directory + 'flat/'
output_clean_directory = output_directory + 'clean/'

# get the field information
field = Configuration.FIELD
date_list = Skopeus.get_field_dates(field)
dates = len(date_list)
field_ra = float(Skopeus.get_field_position(field)[0])
field_dec = float(Skopeus.get_field_position(field)[1])

# run the main pipeline
for dte in range(dates):

	# get the date (e.g. 'yyyy-mm-dd')
	date = date_list[dte]

	# set the date extension (e.g. 'yyyy-mm-dd/FIELD_xx.xxx/')
	date_extension = date + '/FIELD_' + field + '/'

	# set the directories
	full_dark_directory = output_dark_directory + date
	full_flat_directory = output_flat_directory + date
	full_clean_directory = output_clean_directory + date_extension

	background_directory = full_clean_directory + 'bkg/'
	flux_directory = full_clean_directory + 'flx/'
	histogram_directory = full_clean_directory + 'hst/'
	image_directory = full_clean_directory + 'img/'
	photometry_directory = full_clean_directory + 'pht/'
	table_directory = full_clean_directory + 'tbl/'
	variable_directory = full_clean_directory + 'var/'

	# create the output directories
	directory_list = [full_dark_directory, full_flat_directory, full_clean_directory, background_directory, flux_directory, histogram_directory, image_directory, photometry_directory, table_directory, variable_directory]

	for dir in directory_list:
		if not os.path.exists(dir):
			os.makedirs(dir)

	# create the master dark frame (flat exposure)
	Metretes.make_dark(date, 'flat')

	# create the master dark frame (light exposure)
	Metretes.make_dark(date, 'light')

	# create the flatfield frame
	Metretes.make_flat(date)

	# reduce the light frames
	
	raw_frame_directory = Configuration.INPUT_DIRECTORY + 'raw/' + date + '/FIELD_' + field + '/'

	os.chdir(raw_frame_directory)
	raw_frame_list = sorted(glob.glob('*' + Configuration.DATA_EXTENSION))
	
	raw_frame_number = len(raw_frame_list)
	frame_idx = 1

	# set the calibration frame paths
	dark_frame_path = full_dark_directory + '/mdl_' + date + Configuration.DATA_EXTENSION
	flat_frame_path = full_flat_directory + '/ff_' + date + Configuration.DATA_EXTENSION

	# generate a frame table
	frame_table_path = full_clean_directory + 'stat_FIELD_' + field + '_' + date + Configuration.TABLE_EXTENSION
	tbl_col = ('frm', 'dte', 'jd', 'tex', 'bmd', 'bsd', 'rmn', 'rmd', 'rsd')
	tbl_dt = (str, str, float, float, float, float, float, float, float)
	frame_table = Table(names=tbl_col, dtype=tbl_dt)

	# set the reference frame for alignment
	reference_frame = 'cln_' + raw_frame_list[0]
	reference_frame_path = full_clean_directory + reference_frame

	# clean each frame
	for frm in range(raw_frame_number):

		os.chdir(raw_frame_directory)
		
		# set the raw frame name and path
		raw_frame_file = raw_frame_list[frm]
		raw_frame_name = raw_frame_file.split('.')[0] + '.' + raw_frame_file.split('.')[1]

		# set the clean frame name and path
		clean_frame = 'cln_' + raw_frame_file
		clean_frame_name = 'cln_' + raw_frame_name
		clean_frame_path = full_clean_directory + clean_frame

		# set the diagram paths
		bkg_image_path = full_clean_directory + 'bkg/bkg_' + raw_frame_name + Configuration.IMAGE_EXTENSION
		res_image_path = full_clean_directory + 'bkg/res_' + raw_frame_name + Configuration.IMAGE_EXTENSION

		# check if the clean frame exists
		if not os.path.exists(clean_frame_path):

			# open the raw frame
			raw_frame = fits.open(raw_frame_file)
			raw_frame_data = raw_frame[0].data
			raw_frame_header = raw_frame[0].header

			# get the exposure time
			exposure_time = Skopeus.get_exposure_time('light')

			# get the observation time
			obs_time = Time(raw_frame_header['DATE'], format='isot', scale='utc')
			jd = obs_time.jd

			# remove the overscan region
			clean_frame_data, clean_frame_header = Metretes.remove_overscan(raw_frame_data, raw_frame_header)

			# subtract the master dark
			clean_frame_data, clean_frame_header = Metretes.subtract_dark(clean_frame_data, clean_frame_header, dark_frame_path)

			# divide the flatfield
			clean_frame_data, clean_frame_header = Metretes.divide_flat(clean_frame_data, clean_frame_header, flat_frame_path)

			# subtract the background
			clean_frame_data, clean_frame_header, bkg2d_ini_md, bkg2d_ini_sd, bkgsc_ini_mn, bkgsc_ini_md, bkgsc_ini_sd, bkg2d_res_md, bkg2d_res_sd, bkgsc_res_mn, bkgsc_res_md, bkgsc_res_sd = Metretes.subtract_background(clean_frame_data, clean_frame_header, bkg_image_path, res_image_path)

			# solve the field
			fits.writeto(clean_frame_path, clean_frame_data, clean_frame_header)
			Metretes.solve_field(clean_frame, full_clean_directory, field, verbose=Configuration.VERBOSE)

			# add a row to the frame table
			tbl_row = (raw_frame_file, date, jd, exposure_time, bkg2d_ini_md, bkg2d_ini_sd, bkgsc_res_mn, bkgsc_res_md, bkgsc_res_sd)
			frame_table.add_row(tbl_row)

		frame_idx += 1

	# read or save the frame table
	if os.path.exists(frame_table_path):
		frame_table = ascii.read(frame_table_path, format=Configuration.TABLE_FORMAT, delimiter='|')
	else:
		ascii.write(frame_table, frame_table_path, format=Configuration.TABLE_FORMAT)

	# create a master stack frame
	stack_data, stack_header = Metretes.make_stack(date, field, frame_table)

# write the postamble
script_postamble = Boethos.get_script_postamble('preprocessing_pipeline')
print(script_postamble)
fn = time.time()
dt = np.around(fn - st, decimals=2)
print('[' + str(dt) + ' s]')