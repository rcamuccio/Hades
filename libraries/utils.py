from config import Configuration

import glob
import logging
import numpy as np
import os
import shutil

class Utils:

	@staticmethod
	def config_observatory(name):
		''' This function returns a dictionary of parameters related to a given observatory.

		:parameter name [string] - The name of the observatory

		:return params [dictionary] - The parameters of the observatory

			atmosphere_qe
			bandpass 				[nm]
			ccd_qe
			central_wavelength 		[nm]
			collecting_area			[m^2]
			declination_limit		[deg]
			elevation				[m]
			etendue					[m^2 deg^2]
			field_separation		[deg]
			field_size				[deg]
			filter_qe
			focal_length			[m]
			fov						[m^2]
			gain					[e/ADU]
			latitude				[deg]
			longitude				[deg]
			pixel_number
			pixel_scale				[arcsec/px]
			primary_area			[m^2]
			primary_diameter		[m]
			primary_qe
			readout_noise			[e]
			secondary_area			[m^2]
			secondary_diameter		[m]
			secondary_qe
			seeing					[arcsec]
			sky						[mag/arcsec^2]
			telescope
			total_throughput
			vignetting

		'''

		params = {}

		if name == 'CTMO':
			atmosphere_qe = [0.8, 0.8, 0.9, 0.9, 0.9]					# better estimate
			bandpass = [65, 149, 133, 149, 280]
			ccd_qe = 0.9
			central_wavelength = [352.5, 475.5, 628.5, 769.5, 960.0]
			declination_limit = -34.00
			elevation = 11.5
			filter_qe = 0.9
			focal_length = 2.939
			gain = 1.385
			latitude = 25.995789
			longitude = -97.568956
			pixel_number = 4096
			pixel_scale = 0.6305
			primary_diameter = 0.610
			primary_qe = 0.96
			readout_noise = 15.8
			secondary_diameter = 0.190
			secondary_qe = 0.96
			seeing = 4.0												# better estimate
			sky = [19.81, 19.81, 19.81, 19.81, 19.81]					# better estimate
			telescope = 'reflector'

		elif name == 'Macon':
			atmosphere_qe = [0.8, 0.9, 0.9, 0.9]
			bandpass = [147, 141, 147, 147]
			ccd_qe = 0.9
			central_wavelength = [473.5, 638.5, 775.5, 922.5]
			declination_limit = 33.84
			elevation = 4650
			filter_qe = 0.9
			focal_length = 3.974
			gain = 2.18
			latitude = -24.62055556
			longitude = -67.32833333
			pixel_number = 10560
			pixel_scale = 0.468
			primary_diameter = 0.610
			primary_qe = 0.96
			readout_noise = 5.0
			secondary_diameter = 0.280
			secondary_qe = 0.96
			seeing = 0.93
			sky = [22.1, 21.1, 20.1, 18.7]
			telescope = 'reflector'

		elif name == 'OAFA':
			atmosphere_qe = [0.9]
			bandpass = [50]
			ccd_qe = 0.9
			central_wavelength = [530.0]
			declination_limit = 26.66
			elevation = 2420
			filter_qe = 1.0
			focal_length = 3.7
			gain = 2.18
			latitude = -31.8023
			longitude = -69.3265
			pixel_number = 10560
			pixel_scale = 0.4959
			primary_diameter = 0.508
			primary_qe = 0.92
			readout_noise = 5.0
			secondary_diameter = 0.0
			secondary_qe = 1.0
			seeing = 2.5												# better estimate
			sky = [21.96]
			telescope = 'refractor'										# better estimate

		field_size = pixel_scale * pixel_number / 3600

		fov = field_size ** 2

		primary_area = np.pi * (primary_diameter / 2) ** 2
		secondary_area = np.pi * (secondary_diameter / 2) ** 2

		if telescope == 'reflector':
			collecting_area = primary_area - secondary_area
		else:
			collecting_area = primary_area

		etendue = collecting_area * fov

		vignetting_circle = np.pi * (field_size/2) ** 2

		if name == 'Macon':
			vignetting = 0.756
			field_separation = 1.19
		else:
			vignetting = vignetting_circle / fov
			field_separation = vignetting * np.sqrt((field_size/2)**2 + field_size**2)

		total_throughput = ccd_qe * filter_qe * primary_qe * secondary_qe * vignetting

		params['atmosphere_qe'] = atmosphere_qe
		params['bandpass'] = bandpass
		params['ccd_qe'] = ccd_qe
		params['central_wavelength'] = central_wavelength
		params['collecting_area'] = collecting_area
		params['declination_limit'] = declination_limit
		params['elevation'] = elevation
		params['etendue'] = etendue
		params['field_separation'] = field_separation
		params['field_size'] = field_size
		params['filter_qe'] = filter_qe
		params['focal_length'] = focal_length
		params['fov'] = fov
		params['gain'] = gain
		params['latitude'] = latitude
		params['longitude'] = longitude
		params['pixel_number'] = pixel_number
		params['pixel_scale'] = pixel_scale
		params['primary_area'] = primary_area
		params['primary_diameter'] = primary_diameter
		params['primary_qe'] = primary_qe
		params['readout_noise'] = readout_noise
		params['secondary_area'] = secondary_area
		params['secondary_diameter'] = secondary_diameter
		params['secondary_qe'] = secondary_qe
		params['seeing'] = seeing
		params['sky'] = sky
		params['total_throughput'] = total_throughput
		params['vignetting'] = vignetting
		params['vignetting_circle'] = vignetting_circle

		return params

	@staticmethod
	def create_directories(main_dir):

		raw_dir = os.path.join(main_dir, 'raw')
		if not os.path.exists(raw_dir):
			os.mkdir(raw_dir)

		cal_dir = os.path.join(main_dir, 'cal')
		if not os.path.exists(cal_dir):
			os.mkdir(cal_dir)

		wcs_dir = os.path.join(main_dir, 'wcs')
		if not os.path.exists(wcs_dir):
			os.mkdir(wcs_dir)

		align_dir = os.path.join(main_dir, 'align')
		if not os.path.exists(align_dir):
			os.mkdir(align_dir)

		os.chdir(main_dir)

		for item in glob.glob('*.fit'):

			old_item_path = os.path.join(main_dir, item)
			new_item_path = os.path.join(raw_dir, item)

			shutil.move(old_item_path, new_item_path)

		return raw_dir, cal_dir, wcs_dir, align_dir

	@staticmethod
	def log(statement, level):

		log_path = Configuration.MAIN_DIR + 'logs/main.log'

		# create the logger
		logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s', filename=log_path, filemode='a')
		logger = logging.getLogger()

		if not getattr(logger, 'handler_set', None):
			logger.setLevel(logging.DEBUG)

			# create console handler and set level to debug
			ch = logging.StreamHandler()
			ch.setLevel(logging.DEBUG)

			# create the formatter
			formatter = logging.Formatter('%(asctime)s - %(levelname)s: %(message)s')

			# add formatter to ch
			ch.setFormatter(formatter)

			# add ch to logger
			logger.addHandler(ch)

			# 'set' Handler
			logger.handler_set = True

		if level == 'info':
			logger.info(statement)

		if level == 'debug':
			logger.debug(statement)

		if level == 'warning':
			logger.warning(statement)

		if level == 'error':
			logger.error(statement)

		if level == 'critical':
			logger.error(statement)

	@staticmethod
	def setup_hades(main_dir):

		# alerts
		alerts_dir = os.path.join(main_dir, 'alerts/')
		if not os.path.exists(alerts_dir):
			os.mkdir(alerts_dir)

		# alerts > gw
		alerts_gw_dir = os.path.join(alerts_dir, 'gw/')
		if not os.path.exists(alerts_gw_dir):
			os.mkdir(alerts_gw_dir)

		# alerts > gw > json
		alerts_gw_json_dir = os.path.join(alerts_gw_dir, 'json/')
		if not os.path.exists(alerts_gw_json_dir):
			os.mkdir(alerts_gw_json_dir)

		# alerts > gw > moc
		alerts_gw_moc_dir = os.path.join(alerts_gw_dir, 'moc/')
		if not os.path.exists(alerts_gw_moc_dir):
			os.mkdir(alerts_gw_moc_dir)

		# alerts > gw > skymap
		alerts_gw_skymap_dir = os.path.join(alerts_gw_dir, 'skymap/')
		if not os.path.exists(alerts_gw_skymap_dir):
			os.mkdir(alerts_gw_skymap_dir)

		# analysis
		analysis_dir = os.path.join(main_dir, 'analysis/')
		if not os.path.exists(analysis_dir):
			os.mkdir(analysis_dir)

		# analysis > survey
		analysis_survey_dir = os.path.join(analysis_dir, 'survey/')
		if not os.path.exists(analysis_survey_dir):
			os.mkdir(analysis_survey_dir)

		# analysis > gw
		analysis_gw_dir = os.path.join(analysis_dir, 'gw/')
		if not os.path.exists(analysis_gw_dir):
			os.mkdir(analysis_gw_dir)

		# analysis > gw > fields
		analysis_gw_fields_dir = os.path.join(analysis_gw_dir, 'fields/')
		if not os.path.exists(analysis_gw_fields_dir):
			os.mkdir(analysis_gw_fields_dir)

		# analysis > gw > galaxies
		analysis_gw_galaxies_dir = os.path.join(analysis_gw_dir, 'galaxies/')
		if not os.path.exists(analysis_gw_galaxies_dir):
			os.mkdir(analysis_gw_galaxies_dir)

		# logs
		log_dir = os.path.join(main_dir, 'logs/')
		if not os.path.exists(log_dir):
			os.mkdir(log_dir)