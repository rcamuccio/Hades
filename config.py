class Configuration:

	# hardware
	MACHINE = 'epimetheus'
	FILE_EXTENSION = '.fits'
	IMAGE_EXTENSION = '.png'
	TABLE_EXTENSION = '.dat'
	TABLE_FORMAT = 'ascii.fixed_width'
	OBSERVATORY = 'toros'
	VERBOSE = False

	# products
	FIELD = 'FIELD_42.074'
	#DATES = ['2024-11-10', '2024-11-11', '2024-11-12', '2024-11-13', '2024-11-14', '2024-11-15', '2024-11-16', '2024-11-17', '2024-11-18', '2024-11-19', '2024-11-20']
	DATES = ['2024-11-10']

	# directories
	MAIN_DIRECTORY = '/home/epimetheus/Documents/'
	CODE_DIRECTORY = MAIN_DIRECTORY + 'hades/'
	OUTPUT_DIRECTORY = MAIN_DIRECTORY + 'output/'
	DATA_DIRECTORY = '/media/epimetheus/Harbor/'

	INPUT_DATA_DIRECTORY = '/media/epimetheus/692e5e1e-6b16-4928-a25f-46dbabe207e6/toros_data/'
	OUTPUT_DATA_DIRECTORY = '/media/epimetheus/ExtremeSSD/'

	# photometry
	AIRMASS_METHOD = 'ky1994'
	ANNULI_INNER = 18
	ANNULI_OUTER = 20
	APER_SIZE = 16
	BKG_METHOD = 'flat'
	BOX_SIZE = (60, 60)
	COMBINE_METHOD = 'median'
	DATA_TYPE = 'float32'
	FILTER_SIZE = (3, 3)
	FOOTPRINT_RADIUS = 10
	MEMORY_LIMIT = 8e9
	NPIXELS = 3
	SIG_BKG = 2.5
	SIG_SRC = 5.0
	SOLVE_RADIUS = 1.
	UNIT = 'adu'

	# differencing
	AXS_LIMIT = 100 # number of pixels close to edge of frame
	BRIGHT_STARS = 20000 # top stars to search for in kernel stars
	FWHM = 15. # FWHM of the image
	KERNEL_LIMIT = 0.5 # maximum allowable offset in zeropoint in magnitudes
	KRNL = 2 # kernel size = 2 * KRNL + 1
	NRSTARS = 500 # number of stars used to solve for kernel
	ORDR = 1 # order of kernel; 0 = stationary, 1/2 = spatially varying
	PIX = 220 # sky subtraction
	RMS_LOW_LIMIT = 0.005 # lower limit on precision to use for kernel stars
	RMS_UPP_LIMIT = 0.02 # upper limit on precision to use for kernel stars
	STMP = 10 # stamp size = 2 * STMP + 1
	THRESHOLD = 5. # threshold for a source above the background

	# querying
	QUERY_RADIUS = 1.
	QUERY_SOURCE = 'gaia-cone'
	ROW_LIMIT = -1

	# plotting
	CMAP = 'gray'
	DPI = 300
	FIGURE_SIZE = (10, 10)
	FONT_NAME = 'Monospace'
	FONT_SIZE = 12
	HISTOGRAM_BINS = 'fd'
	HISTOGRAM_LIMIT = 5
	HISTOGRAM_SCALE = 'log'
	HISTOGRAM_TYPE = 'step'
	INTERVAL = 'zscale'
	SAVE_FIGURE = True

	# listener
	CATALOG = 'glade+'
	CLIENT_ID = '30mai3vn918g150vbdm63crtut'
	CLIENT_SECRET = '16895pg64pqb556jmsqdmj1nsjvvm37risqvmfn1329tcorpb6ht'
	EMAIL = 'toros.alerts@yahoo.com'
	FIELD_GENERATION = 'N'
	FIELD_PRIORITY = ['commissioning', 'ligo', 'science', 'main_survey']
	GALACTIC_PLANE = 15.
	GALAXY_LIST = 50
	LISTEN_NED_WAIT = 1
	MAILING_LIST = ['rcamuccio@gmail.com']
	MOON_DISTANCE = 60.
	PAS = 'ncnrlqwdthofhoch'
	PORT = 465
	SMTP = 'smtp.mail.yahoo.com'

	AVAILABLE_TOPICS = ['gcn.circulars', 'gcn.heartbeat', 'gcn.notices.chime.frb', 'gcn.notices.dsa110.frb', 'gcn.notices.einstein_probe.wxt.alert', 'gcn.notices.icecube.gold_bronze_track_alerts', 'gcn.notices.icecube.lvk_nu_track_search', 'gcn.notices.superk.sn_alert', 'gcn.notices.swift.bat.guano', 'igwn.gwalert']
	
	# ARCHIVE

	#KNOWN_VARIABLES = 'Y'
	#BAD_DATES = '2024-10-03'
	#NSKY_STARS = 20000

	# pipeline toggles
	#CLEAN_SKIP = 'N'
	#WRITE_SKY = 'N'
	#MASTER_SKIP = 'N'
	#DIFFERENCE_SKIP = 'N'
	#PHOTOMETRY_SKIP = 'N'
	#LIGHTCURVE_SKIP = 'N'
	#SYSTEMATICS_SKIP = 'N'

	# reduction toggles
	#SUBTRACT_BIAS = 'Y'
	#SUBTRACT_DARK = 'Y'
	#DIVIDE_FLAT = 'Y'
	#CLIP_IMAGE = 'Y'
	#SUBTRACT_SKY = 'Y'
	#PLATE_SOLVE = 'Y'
	#WRITE_SKY = 'Y'
	#SOURCE_EXTRACTOR_PATH = '/usr/bin/source-extractor'