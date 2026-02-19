from config import Configuration
from astropy.coordinates import AltAz, SkyCoord
from scipy.optimize import curve_fit
import logging
import numpy as np
import os
import pandas as pd
import sys
logging.getLogger('gcn').setLevel(logging.WARNING)
logging.getLogger('healpy').setLevel(logging.WARNING)
logging.getLogger('requests').setLevel(logging.WARNING)
logging.getLogger('urllib3').setLevel(logging.WARNING)

class Utils:

	@staticmethod
	def angular_distance(ra1, de1, ra2, de2):
		'''This function determines the angular distance between two sky coordinates.

		:parameter ra1 - The right ascension of point 1
		:parameter de1 - The declination of point 1
		:parameter ra2 - The right ascension of point 2
		:parameter de2 - The declination of point 2

		:return ang_dist_deg - The angular distance between the two points [deg]
		'''

		# convert to radians
		ra1_rad = np.deg2rad(ra1)
		de1_rad = np.deg2rad(de1)
		ra2_rad = np.deg2rad(ra2)
		de2_rad = np.deg2rad(de2)

		# the angular distance in radians
		ang_dist_rad = np.arccos((np.sin(de1_rad) * np.sin(de2_rad)) + (np.cos(de1_rad) * np.cos(de2_rad) * np.cos(ra1_rad - ra2_rad)))

		# the angular distance in degrees
		ang_dist_deg = np.rad2deg(ang_dist_rad)

		return ang_dist_deg

	@staticmethod
	def calc_airmass(h, method='ky1994'):

		# convert altitude to zenith distance
		z = 90 - h

		# Pickering (2002)
		if method == 'p2002':
			x = 1 / np.sin(h + (244/(47*h**1.1)))

		# Kasten and Young (1994)
		else:
			x = 1 / np.cos(np.deg2rad(z)) * 0.50572*((6.07995 + 90 - z)**-1.6364)

		return x

	@staticmethod
	def calc_altitude(location, time, ra, de):

		aa = AltAz(location=location, obstime=time)
		coord = SkyCoord(str(ra), str(de), unit='deg')
		coord_transform = coord.transform_to(aa)

		altitude = coord_transform.alt.degree

		return altitude

	@staticmethod
	def config_data_dir(create=True):

		dirs = {}
		data_directory = Configuration.DATA_DIRECTORY

		raw_directory = data_directory + 'raw/'
		clean_directory = data_directory + 'clean/'
		master_main_directory = data_directory + 'master/'
		master_directory = master_main_directory + Configuration.FIELD + '/'
		master_tmp_directory = master_directory + 'tmp_master/'
		centroid_directory = master_directory + 'centroids/'
		calibration_directory = data_directory + 'calibration/'
		bias_directory = calibration_directory + 'tmp_bias/'
		flat_directory = calibration_directory + 'tmp_flat/'
		dark_directory = calibration_directory + 'tmp_dark/'
		lightcurve_directory = data_directory + 'lc/'
		lightcurve_field_directory = lightcurve_directory + Configuration.FIELD + '/'
		differenced_directory = data_directory + 'diff/'
		flux_directory = data_directory + 'flux/'
		review_directory = data_directory + 'review/'

		dirs['raw'] = raw_directory
		dirs['clean'] = clean_directory
		dirs['master_main'] = master_main_directory
		dirs['master'] = master_directory
		dirs['master_tmp'] = master_tmp_directory
		dirs['centroid'] = centroid_directory
		dirs['calibration'] = calibration_directory
		dirs['bias'] = bias_directory
		dirs['flat'] = flat_directory
		dirs['dark'] = dark_directory
		dirs['lightcurve'] = lightcurve_directory
		dirs['lightcurve_field'] = lightcurve_field_directory
		dirs['differenced'] = differenced_directory
		dirs['flux'] = flux_directory
		dirs['review'] = review_directory

		dir_list = [raw_directory, clean_directory, master_main_directory, master_directory, master_tmp_directory, centroid_directory, calibration_directory, bias_directory, flat_directory, dark_directory, lightcurve_directory, lightcurve_field_directory, differenced_directory, flux_directory, review_directory]

		if create:
			Utils.create_directories(dir_list)

		return dirs

	@staticmethod
	def config_field(field_id):
		'''This function returns the right ascension and declination of the given field ID.

		:parameter field_id [str] - The hexadecimal code for the field

		:return ra [float] - The right ascension of the field [deg]
		:return de [float] - The declination of the field [deg]
		'''

		output_dirs = Utils.config_output_dir(create=False)

		# read the field list
		field_path = output_dirs['analysis'] + 'toros_fields.dat'
		field_list = pd.read_csv(field_path, header=0, sep=' ')

		# find the field
		for i in range(len(field_list.field_id)):
			if field_id == field_list.field_id[i]:
				indx = i

		# return the field position
		ra = float(field_list.iloc[indx].ra)
		de = float(field_list.iloc[indx].dec)

		return ra, de

	@staticmethod
	def config_observatory(name):
		'''This function returns a dictionary of parameters related to a given observatory.

		:parameter name [str] - The name of the observatory

		:return params [dict] - The parameters of the observatory
		'''

		params = {}

		if name == 'toros':
			axs_x = 10560
			axs_y = 10560
			axs_x_rw = 12000
			axs_y_rw = 10560
			bp = [147, 141, 147, 147]
			cwl = [473.5, 638.5, 775.5, 922.5]
			dec_limit = 26.66
			elevation = 2420
			exp_light = 300.
			exp_dark = 300.
			exp_flat = 5.
			gain = 0.380
			latitude = -31.8023
			longitude = -69.3265
			mirror_d = 0.610
			num_exp = 1
			num_pix = 10560
			overhead = 30.
			ovs_x = 180
			ovs_y = 20
			peakmax = 45000.
			pixel_scale = 0.4959
			qe_atm = [0.8, 0.9, 0.9, 0.9]
			qe_ccd = 0.85
			qe_fil = 0.9
			qe_pri = 0.96
			qe_sec = 0.96
			read_noise = 5.0
			read_time = 90.
			seeing = 0.93
			sky = [22.1, 21.1, 20.1, 18.7]
			tap_x = 1320
			tap_y = 5280
			utc = -3
			vignetting = 0.756

		exp_day = exp_light / 60. / 60. / 24.
		exp_tot = (exp_light + read_time) * num_exp + overhead
		fov = (pixel_scale * num_pix) / 3600.
		mirror_r = mirror_d / 2
		search_dist = fov
		throughput = qe_ccd * qe_fil * qe_pri * qe_sec * vignetting

		params['axs_x'] = axs_x
		params['axs_y'] = axs_y
		params['axs_x_rw'] = axs_x_rw
		params['axs_y_rw'] = axs_y_rw
		params['bp'] = bp
		params['cwl'] = cwl
		params['dec_limit'] = dec_limit
		params['elevation'] = elevation
		params['exp_day'] = exp_day
		params['exp_light'] = exp_light
		params['exp_dark'] = exp_dark
		params['exp_flat'] = exp_flat
		params['exp_tot'] = exp_tot
		params['fov'] = fov
		params['gain'] = gain
		params['latitude'] = latitude
		params['longitude'] = longitude
		params['mirror_d'] = mirror_d
		params['mirror_r'] = mirror_r
		params['num_exp'] = num_exp
		params['num_pix'] = num_pix
		params['overhead'] = overhead
		params['ovs_x'] = ovs_x
		params['ovs_y'] = ovs_y
		params['peakmax'] = peakmax
		params['pixel_scale'] = pixel_scale
		params['qe_atm'] = qe_atm
		params['qe_ccd'] = qe_ccd
		params['qe_fil'] = qe_fil
		params['qe_pri'] = qe_pri
		params['qe_sec'] = qe_sec
		params['read_noise'] = read_noise
		params['read_time'] = read_time
		params['search_dist'] = search_dist
		params['seeing'] = seeing
		params['sky'] = sky
		params['tap_x'] = tap_x
		params['tap_y'] = tap_y
		params['throughput'] = throughput
		params['utc'] = utc
		params['vignetting'] = vignetting

		return params

	@staticmethod
	def config_output_dir(create=True):

		dirs = {}
		output_directory = Configuration.OUTPUT_DIRECTORY

		ale_dir = output_directory + 'alerts/'
		ana_dir = output_directory + 'analysis/'
		cat_dir = output_directory + 'catalogs/'
		dif_dir = output_directory + 'difference/'
		log_dir = output_directory + 'logs/'
		sou_dir = output_directory + 'sources/'

		dirs['alerts'] = ale_dir
		dirs['analysis'] = ana_dir
		dirs['catalogs'] = cat_dir
		dirs['difference'] = dif_dir
		dirs['logs'] = log_dir
		dirs['sources'] = sou_dir

		dir_list = [ale_dir, ana_dir, cat_dir, dif_dir, log_dir, sou_dir]

		if create:
			Utils.create_directories(dir_list)

		return dirs

	@staticmethod
	def create_directories(path_list):
		'''This function checks for each directory in the list, and creates it, if it doesn't already exist.

		:parameter path_list - The list of directories to create

		:return - Nothing is returned, but directories are created, if necessary
		'''

		for path in path_list:
			if not os.path.exists(path):
				os.mkdir(path)
				Utils.log(path + ' created.', 'info')
			else:
				Utils.log(path + ' already exists. Skipping.', 'info')

	@staticmethod
	def get_all_files_per_field(path, field, step, file_ext):
		'''This function will return all files for a given field without their extension.

		:parameter path [str] - The path of the file
		:parameter field [str] - The field
		:parameter step [str] - The step of the process
		:parameter file_ext [str] - The type of file

		:return files_no_ext [list] - A list of files without their extension
		:return uni_dte_dir [list] - A list of unique dates
		'''

		# get the files in the path with the given file extension
		files_no_ext = []
		dte_dir = []
		no_dte_list = []
		
		sub_dir = os.listdir(path)

		for f in sub_dir:
			if step == 'raw':
				z = path + f + '/' + field + '/'
			else:
				z = path + f + '/' + field + '/'
			try:
				for x in os.listdir(z):
					if ((x.split('_')[0] + '_' + x.split('_')[1] == field) | (x.split('_')[0] == field)):
						if x.endswith(file_ext):
							files_no_ext.append(z + x)
			except:
				no_dte_list.append(f)
		no_dte_list = sorted(no_dte_list)
		
		# get the unique dates with the field
		for file in files_no_ext:
			dte_dir.append(file.split('/')[-3])
		uni_dte_dir = np.unique(dte_dir).tolist()

		'''
		Utils.log(str(field) + ' not observed on:', 'info')
		for dte in no_dte_list:
			Utils.log('  ' + str(dte), 'info')

		Utils.log(str(field) + ' observed on:', 'info')
		for dte in uni_dte_dir:
			Utils.log('  ' + str(dte), 'info')
		'''

		return files_no_ext, uni_dte_dir

	@staticmethod
	def get_file_list(path, file_ext):
		'''This function will return the files in a given directory without their extension.

		:parameter path [str] - The path of the file
		:parameter file_ext [str] - The file type to make a list (e.g. *.fits)

		:return file_list [list] - A sorted list of files without extensions
		'''

		# get the files in the path with the given file extension
		file_list = [f for f in os.listdir(path) if f.endswith(file_ext)]

		# sort based on the number of the image
		file_list.sort(key=len)

		return file_list

	@staticmethod
	def log(statement, level):
		'''This function logs all activity from the program to both screen and file.

		:parameter statement [str] - The statement to log
		:parameter level [str] - The type of statement

		:return - Nothing is returned, but the log is updated and printed to the screen
		'''

		output_dirs = Utils.config_output_dir(create=False)

		# create the log
		log_path = output_dirs['logs'] + 'toros.log'
		if not os.path.isfile(log_path):
			with open(log_path, 'a') as f:
				f.write('|--- This is the beginning of the log. ---|\n')
		else:
			pass

		# create the logger
		logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s', filename=log_path, filemode='a')
		logger = logging.getLogger()

		if not getattr(logger, 'handler_set', None):
			logger.setLevel(logging.DEBUG)

			# create console handler and set level to 'debug'
			ch = logging.StreamHandler()
			ch.setLevel(logging.DEBUG)

			# create the formatter
			formatter = logging.Formatter('%(asctime)s - %(levelname)s: %(message)s')

			# add formatter to console handler
			ch.setFormatter(formatter)

			# add console handler to logger
			logger.addHandler(ch)

			# set the handler
			logger.handler_set = True

		# define log statement levels
		if level == 'info':
			logger.info(statement)

		if level == 'debug':
			logger.debug(statement)

		if level == 'warning':
			logger.warning(statement)

		if level == 'error':
			logger.error(statement)

		if level == 'critical':
			logger.critical(statement)

	@staticmethod
	def sigma_clip(xlist, ylist, sigma):

		xarr = np.asarray(xlist)
		yarr = np.asarray(ylist)

		xm = np.mean(xarr)
		ym = np.mean(yarr)

		xs = np.std(xarr)
		ys = np.std(yarr)

		new_xlist = []
		new_ylist = []

		for i in range(len(xlist)):
			if ylist[i] < (ym - (sigma * ys)):
				pass
			else:
				new_xlist.append(xlist[i])
				new_ylist.append(ylist[i])

		return new_xlist, new_ylist

	@staticmethod
	def unweighted_fit(xlist, ylist):

		xarr = np.asarray(xlist)
		yarr = np.asarray(ylist)

		def f(x, m, b):
			y = (m*x) + b
			return y

		popt, pcov = curve_fit(f, xarr, yarr)
		yfit = f(xarr, *popt)

		m = popt[0]
		dm = np.sqrt(pcov[0][0])

		b = popt[1]
		db = np.sqrt(pcov[1][1])

		return yfit, m, dm, b, db

	@staticmethod
	def write_txt(path, typ, line):
		'''This function will either write or append a text file, primarily used for results or logging.

		:parameter path: The location for the file
		:parameter typ: The write type either 'w', or 'a'
		:parameter line: The line you want to write to the file

		:return: Nothing is returned, but the file is written or appended.
		'''

		# open file
		file = open(path, typ)

		# write line to file
		file.write(line)

		# close file
		file.close()