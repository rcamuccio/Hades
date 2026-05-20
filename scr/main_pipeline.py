from config import Configuration
from lib.photometry import Photometry
from lib.query import Query
import glob
import numpy as np
import os
import time
import warnings
warnings.simplefilter('ignore')

os.system('clear')
print('[\033[1m' + 'Πλούτων' + '\033[0m]')
print('\nRunning main pipeline\n')
st = time.time()

# configure input directories
input_dir = '/media/epimetheus/692e5e1e-6b16-4928-a25f-46dbabe207e6/toros_data/'
dark_dir = input_dir + 'dark/'
flat_dir = input_dir + 'flat/'
raw_dir = input_dir + 'raw/'

# configure output directories
output_dir = '/media/epimetheus/ExtremeSSD/'
md_dir = output_dir + 'dark/'
mf_dir = output_dir + 'flat/'
cf_dir = output_dir + 'clean/'
output_dir_list = [md_dir, mf_dir, cf_dir]
for dr in output_dir_list:
	if not os.path.exists(dr):
		os.mkdir(dr)

# grab the dates
field = Configuration.FIELD
date_list = Configuration.DATES
num_dates = len(date_list)

# run the main pipeline
for dte in range(num_dates):
	date = date_list[dte]
	output_dir = Configuration.OUTPUT_DATA_DIRECTORY + 'clean/' + date + '/' + field + '/'
	qry_tbl_path = output_dir + 'tbl_qry' + Configuration.TABLE_EXTENSION
	src_tbl_path = output_dir + 'tbl_src' + Configuration.TABLE_EXTENSION
	msr_tbl_path = output_dir + 'tbl_msr' + Configuration.TABLE_EXTENSION

	# make a master (flat) dark
	Photometry.make_dark(date, 'flat')

	# make a master (light) dark
	Photometry.make_dark(date, 'light')

	# make a normalized flatfield
	Photometry.make_flat(date)

	# reduce the frames
	frm_tbl = Photometry.clean_raw_frames(date, field)

	# make a stack
	stk_data, stk_header = Photometry.make_stack(date, field, frm_tbl)

	# create a queried source catalog
	qry_tbl = Query.gaia_cone(stk_header['CRVAL1'], stk_header['CRVAL2'], Configuration.QUERY_RADIUS, qry_tbl_path)

	# create an extracted source catalog
	src_tbl = Photometry.extract_sources(stk_data, stk_header, src_tbl_path)

	# match the catalogs
	msr_tbl = Photometry.match_catalogs(src_tbl, qry_tbl, msr_tbl_path)

fn = time.time()
dt = np.around(fn - st, decimals=2)
print('\nMain pipeline finished in', dt, 's')