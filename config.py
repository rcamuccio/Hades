class Configuration:

	# configuration
	MACHINE = 'epimetheus' # epimetheus, qu
	DATE = ''

	# directories
	HADES_DIRECTORY = '/home/' + MACHINE + '/Documents/hades/'
	TARTARUS_DIRECTORY = '/home/' + MACHINE + '/Documents/tartarus/'
	PERSEPHONE_DIRECTORY = '/home/' + MACHINE + '/Documents/persephone/'

	ALERTS_DIRECTORY = TARTARUS_DIRECTORY + 'alerts/'
	ANALYSIS_DIRECTORY = TARTARUS_DIRECTORY + 'analysis/'
	CATALOGS_DIRECTORY = TARTARUS_DIRECTORY + 'catalogs/'
	LOGS_DIRECTORY = TARTARUS_DIRECTORY + 'logs/'

	ALERTS_GW_DIRECTORY = ALERTS_DIRECTORY + 'gw/'
	ALERTS_GW_JSON_DIRECTORY = ALERTS_GW_DIRECTORY + 'json/'
	ALERTS_GW_MOC_DIRECTORY = ALERTS_GW_DIRECTORY + 'moc/'
	ALERTS_GW_PARAM_DIRECTORY = ALERTS_GW_DIRECTORY + 'param/'
	ALERTS_GW_SKYMAP_DIRECTORY = ALERTS_GW_DIRECTORY + 'skymap/'

	ANALYSIS_SURVEY_DIRECTORY = ANALYSIS_DIRECTORY + 'survey/'
	ANALYSIS_GW_DIRECTORY = ANALYSIS_DIRECTORY + 'gw/'
	ANALYSIS_GW_FIELDS_DIRECTORY = ANALYSIS_GW_DIRECTORY + 'fields/'
	ANALYSIS_GW_GALAXIES_DIRECTORY = ANALYSIS_GW_DIRECTORY + 'galaxies/'

	TARTARUS_DIRECTORIES = [TARTARUS_DIRECTORY, ALERTS_DIRECTORY, ANALYSIS_DIRECTORY, 
							CATALOGS_DIRECTORY, LOGS_DIRECTORY, ALERTS_GW_DIRECTORY, 
							ALERTS_GW_JSON_DIRECTORY, ALERTS_GW_MOC_DIRECTORY, ALERTS_GW_PARAM_DIRECTORY, 
							ALERTS_GW_SKYMAP_DIRECTORY, ANALYSIS_SURVEY_DIRECTORY, ANALYSIS_GW_DIRECTORY, 
							ANALYSIS_GW_FIELDS_DIRECTORY, ANALYSIS_GW_GALAXIES_DIRECTORY]

	DARKS_DIRECTORY = PERSEPHONE_DIRECTORY + 'darks/' + DATE + '/'
	BIAS_DIRECTORY = PERSEPHONE_DIRECTORY + 'bias/' + DATE + '/'
	FLATS_DIRECTORY = PERSEPHONE_DIRECTORY + 'flats/' + DATE + '/'
	RAW_DIRECTORY = PERSEPHONE_DIRECTORY + 'raw/'
	CLEAN_DIRECTORY = PERSEPHONE_DIRECTORY + 'clean/'
	MASTER_DIRECTORY = PERSEPHONE_DIRECTORY + 'master/'
	CALIBRATION_DIRECTORY = PERSEPHONE_DIRECTORY + 'calibration/'
	CENTROID_DIRECTORY = MASTER_DIRECTORY + 'centroids/'
	LIGHTCURVE_DIRECTORY = PERSEPHONE_DIRECTORY + 'lc/'
	DIFFERENCED_DIRECTORY = PERSEPHONE_DIRECTORY + 'diff/'
	CLEAN_DATE_DIRECTORY = CLEAN_DIRECTORY + DATE + '/'

	PERSEPHONE_DIRECTORIES = [PERSEPHONE_DIRECTORY, DARKS_DIRECTORY, BIAS_DIRECTORY, FLATS_DIRECTORY,
								RAW_DIRECTORY, CLEAN_DIRECTORY, MASTER_DIRECTORY, CALIBRATION_DIRECTORY,
								CENTROID_DIRECTORY, LIGHTCURVE_DIRECTORY, DIFFERENCED_DIRECTORY,
								CLEAN_DATE_DIRECTORY]


	# survey
	CATALOG_TARGET = 'glade+' # glade24, glade+
	EXPOSURE_TIME = 300.
	FIELD_PRIORITY = ['commissioning', 'ligo', 'science', 'main_survey']
	GALACTIC_PLANE = 15.
	MOON_DISTANCE = 60.
	OBSERVATORY = 'OAFA' # CTMO, Macon, OAFA
	OVERHEAD = 30.
	READ_TIME = 100.

	# photometry
	AIRMASS_METHOD = 'ky1998'
	BKG_METHOD = 'flat'
	BOX_SIZE = (50, 50)
	CATALOG = 'gaia-cone'
	CATALOG_F1 = 'r'
	CATALOG_F2 = 'i'
	COMBINE_METHOD = 'median'
	DILATE_SIZE = 25
	DTYPE = 'float32'
	FIELD_RA = 326.1413
	FIELD_DEC = 38.5948
	FILTER_SIZE = (3, 3)
	MEM_LIMIT = 32e9
	NPIXELS = 3
	OBJECT = 'grb240615a'
	PHOT_F1 = 'r'
	PHOT_F2 = 'i'
	RAD_AN_IN = 17
	RAD_AN_OUT = 20
	RAD_AP = 14
	RAD_SOLVE = 1
	RAD_QUERY = 1
	SIGMA_BKG = 3.0
	SIGMA_SRC = 5.0
	WORKING_DIR = '/home/epimetheus/Downloads/2024-06-16/grb240615a/align/'

	# alerts
	CLIENT_ID = 'j1hrnlcur5mc2mkch0a037s98'
	CLIENT_SECRET = '1ojbdiie5bl6ccsf369b74hk0v0lcuemaeq32r4qa9dor9up4cp0'
	EMAIL = 'toros.alerts@yahoo.com'
	MAILING_LIST = ['rcamuccio@gmail.com', 'moises.castillo01@utrgv.edu', 'moemyself3@gmail.com']
	PAS = 'ncnrlqwdthofhoch'
	PORT = 465
	SMTP = 'smtp.mail.yahoo.com'

	AVAILABLE_TOPICS = ['igwn.gwalert',										# #
						'gcn.notices.icecube.lvk_nu_track_search',			# #
						'gcn.notices.swift.bat.guano',						# 
						'gcn.classic.text.FERMI_GBM_ALERT',					# #
						'gcn.classic.text.FERMI_GBM_FIN_POS',				# #
						'gcn.classic.text.FERMI_GBM_FLT_POS',				# #
						'gcn.classic.text.FERMI_GBM_GND_POS',				# #
						'gcn.classic.text.FERMI_GBM_POS_TEST',				# #
						'gcn.classic.text.FERMI_GBM_SUBTHRESH',				# #
						'gcn.classic.text.FERMI_LAT_MONITOR',				# 
						'gcn.classic.text.FERMI_LAT_OFFLINE',
						'gcn.classic.text.FERMI_LAT_POS_TEST',				# 
						'gcn.classic.text.FERMI_POINTDIR',					# #
						'gcn.classic.text.ICECUBE_ASTROTRACK_BRONZE',
						'gcn.classic.text.ICECUBE_ASTROTRACK_GOLD',
						'gcn.classic.text.ICECUBE_CASCADE',
						'gcn.classic.text.SWIFT_ACTUAL_POINTDIR',			# #
						'gcn.classic.text.SWIFT_BAT_GRB_LC',
						'gcn.classic.text.SWIFT_BAT_GRB_POS_ACK',
						'gcn.classic.text.SWIFT_BAT_GRB_POS_TEST',			# #
						'gcn.classic.text.SWIFT_BAT_QL_POS',
						'gcn.classic.text.SWIFT_BAT_SCALEDMAP',
						'gcn.classic.text.SWIFT_BAT_TRANS',
						'gcn.classic.text.SWIFT_FOM_OBS',
						'gcn.classic.text.SWIFT_POINTDIR',					# #
						'gcn.classic.text.SWIFT_SC_SLEW',
						'gcn.classic.text.SWIFT_TOO_FOM',					# 
						'gcn.classic.text.SWIFT_TOO_SC_SLEW',				# 
						'gcn.classic.text.SWIFT_UVOT_DBURST',
						'gcn.classic.text.SWIFT_UVOT_DBURST_PROC',
						'gcn.classic.text.SWIFT_UVOT_EMERGENCY',
						'gcn.classic.text.SWIFT_UVOT_FCHART',
						'gcn.classic.text.SWIFT_UVOT_FCHART_PROC',
						'gcn.classic.text.SWIFT_UVOT_POS',
						'gcn.classic.text.SWIFT_UVOT_POS_NACK',
						'gcn.classic.text.SWIFT_XRT_CENTROID',
						'gcn.classic.text.SWIFT_XRT_IMAGE',
						'gcn.classic.text.SWIFT_XRT_IMAGE_PROC',
						'gcn.classic.text.SWIFT_XRT_LC',
						'gcn.classic.text.SWIFT_XRT_POSITION',
						'gcn.classic.text.SWIFT_XRT_SPECTRUM',
						'gcn.classic.text.SWIFT_XRT_SPECTRUM_PROC',
						'gcn.classic.text.SWIFT_XRT_SPER',
						'gcn.classic.text.SWIFT_XRT_SPER_PROC',
						'gcn.classic.text.SWIFT_XRT_THRESHPIX',
						'gcn.classic.text.SWIFT_XRT_THRESHPIX_PROC']

	UNAVAILABLE_TOPICS = ['gcn.classic.text.FERMI_GBM_LC',
						'gcn.classic.text.FERMI_GBM_TRANS',
						'gcn.classic.text.FERMI_LAT_GND',
						'gcn.classic.text.FERMI_LAT_POS_DIAG',
						'gcn.classic.text.FERMI_LAT_POS_INI',
						'gcn.classic.text.FERMI_LAT_POS_UPD',
						'gcn.classic.text.FERMI_LAT_TRANS',	
						'gcn.classic.text.FERMI_SC_SLEW',
						'gcn.classic.text.SWIFT_BAT_ALARM_LONG',
						'gcn.classic.text.SWIFT_BAT_ALARM_SHORT',
						'gcn.classic.text.SWIFT_BAT_GRB_ALERT',
						'gcn.classic.text.SWIFT_BAT_GRB_LC_PROC',
						'gcn.classic.text.SWIFT_BAT_GRB_POS_NACK',
						'gcn.classic.text.SWIFT_BAT_KNOWN_SRC',
						'gcn.classic.text.SWIFT_BAT_MONITOR',
						'gcn.classic.text.SWIFT_BAT_SLEW_POS',
						'gcn.classic.text.SWIFT_BAT_SUB_THRESHOLD',
						'gcn.classic.text.SWIFT_BAT_SUBSUB',
						'gcn.classic.text.SWIFT_FOM_PPT_ARG_ERR',
						'gcn.classic.text.SWIFT_FOM_SAFE_POINT',
						'gcn.classic.text.SWIFT_FOM_SLEW_ABORT',
						'gcn.classic.text.SWIFT_XRT_EMERGENCY']
