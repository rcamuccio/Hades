class Configuration:

	WORKING_DIR = '/home/epimetheus/Downloads/2024-06-16/grb240615a/align/'

	MACHINE = 'qu' # epimetheus, qu
	MAIN_DIR = '/home/' + MACHINE + '/Documents/hades-main/'
	CATALOG_DIR = '/home/' + MACHINE + '/Documents/catalogs/'

	OBJECT = 'grb240615a'
	FIELD_RA = 326.1413
	FIELD_DEC = 38.5948
	
	OBSERVATORY = 'OAFA' # CTMO, Macon, OAFA

	EXPOSURE_TIME = 300
	FIELD_PRIORITY = ['commissioning', 'ligo', 'science', 'main_survey']
	GALACTIC_PLANE = 15
	MOON_DISTANCE = 60
	OVERHEAD = 30.
	READ_TIME = 100.

	AIRMASS_METHOD = 'ky1998'
	BKG_METHOD = 'flat'
	BOX_SIZE = (50, 50)
	CATALOG = 'gaia-cone'
	CATALOG_TARGET = 'glade24' # glade24, glade+
	CATALOG_F1 = 'r'
	CATALOG_F2 = 'i'
	COMBINE_METHOD = 'median'
	DILATE_SIZE = 25
	DTYPE = 'float32'
	FILTER_SIZE = (3, 3)
	MEM_LIMIT = 32e9
	NPIXELS = 3
	PHOT_F1 = 'r'
	PHOT_F2 = 'i'
	RAD_AN_IN = 17
	RAD_AN_OUT = 20
	RAD_AP = 14
	RAD_SOLVE = 1
	RAD_QUERY = 1
	SIGMA_BKG = 3.0
	SIGMA_SRC = 5.0

	CLIENT_ID = 'client_id'
	CLIENT_SECRET = 'client_secret'

	AVAILABLE_TOPICS = ['igwn.gwalert',
						'gcn.notices.icecube.lvk_nu_track_search',
						'gcn.notices.swift.bat.guano',
						'gcn.classic.text.FERMI_GBM_ALERT',
						'gcn.classic.text.FERMI_GBM_FIN_POS',
						'gcn.classic.text.FERMI_GBM_FLT_POS',
						'gcn.classic.text.FERMI_GBM_GND_POS',
						'gcn.classic.text.FERMI_GBM_POS_TEST',
						'gcn.classic.text.FERMI_GBM_SUBTHRESH',
						'gcn.classic.text.FERMI_LAT_MONITOR',
						'gcn.classic.text.FERMI_LAT_OFFLINE',
						'gcn.classic.text.FERMI_LAT_POS_TEST',
						'gcn.classic.text.FERMI_POINTDIR',
						'gcn.classic.text.ICECUBE_ASTROTRACK_BRONZE',
						'gcn.classic.text.ICECUBE_ASTROTRACK_GOLD',
						'gcn.classic.text.ICECUBE_CASCADE',
						'gcn.classic.text.SWIFT_ACTUAL_POINTDIR',
						'gcn.classic.text.SWIFT_BAT_GRB_LC',
						'gcn.classic.text.SWIFT_BAT_GRB_POS_ACK',
						'gcn.classic.text.SWIFT_BAT_GRB_POS_TEST',
						'gcn.classic.text.SWIFT_BAT_QL_POS',
						'gcn.classic.text.SWIFT_BAT_SCALEDMAP',
						'gcn.classic.text.SWIFT_BAT_TRANS',
						'gcn.classic.text.SWIFT_FOM_OBS',
						'gcn.classic.text.SWIFT_POINTDIR',
						'gcn.classic.text.SWIFT_SC_SLEW',
						'gcn.classic.text.SWIFT_TOO_FOM',
						'gcn.classic.text.SWIFT_TOO_SC_SLEW',
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