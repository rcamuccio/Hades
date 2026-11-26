class Configuration:

	# hardware
	MACHINE = 'epimetheus'
	DATA_EXTENSION = '.fits'
	IMAGE_EXTENSION = '.png'
	TABLE_EXTENSION = '.dat'
	TABLE_FORMAT = 'fixed_width'
	OBSERVATORY = 'toros'
	VERBOSE = False

	# products
	FIELD = '36.007'
	SOURCE_RA = 10.1958012
	SOURCE_DEC = -21.9258309

	# paths
	HADES_DIRECTORY = '/home/epimetheus/Documents/hades/'
	PLOUTON_DIRECTORY = '/home/epimetheus/Documents/plouton/'
	INPUT_DIRECTORY = '/media/epimetheus/692e5e1e-6b16-4928-a25f-46dbabe207e6/toros_data/'
	OUTPUT_DIRECTORY = '/media/epimetheus/ExtremeSSD/'

	# photometry
	AIRMASS_METHOD = 'ky1989'
	ANNULUS_INNER = 22
	ANNULUS_OUTER = 25
	APERTURE_SIZE = 20
	BACKGROUND_METHOD = '2d'
	BOX_SIZE = (60, 60)
	COMBINE_METHOD = 'median'
	DATA_TYPE = 'float32'
	FILTER_SIZE = (3, 3)
	FLAT_MAXIMUM = 40000.
	FLAT_MINIMUM = 3000.
	FOOTPRINT_RADIUS = 10
	KEEP_ALIGN_FRAMES = False
	MATCH_TOLERANCE = 5
	MEMORY_LIMIT = 8e9
	MINIMUM_PIXELS = 3
	SIGMA_BACKGROUND = 2.5
	SIGMA_SOURCE = 5.
	SOLVE_RADIUS = 1.
	UNIT = 'adu'

	# differencing
	AXIS_LIMIT = 100
	BRIGHT_STARS = 20000
	FWHM = 15.
	KERNEL_LIMIT = 0.5
	KERNEL_ORDER = 1
	KERNEL_STARS = 500
	KRNL = 2
	PIX = 220
	RMS_LOWER_LIMIT = 0.005
	RMS_UPPER_LIMIT = 0.02
	STMP = 10
	THRESHOLD = 5.

	# querying
	QUERY_RADIUS = 1.
	QUERY_SOURCE = 'gaia'
	ROW_LIMIT = -1

	# plotting
	COLORMAP = 'gray'
	DPI = 250
	FIGURE_SIZE = (10, 10)
	FONT_NAME = 'Monospace'
	FONT_SIZE = 12
	HISTOGRAM_BINS = 'fd'
	HISTOGRAM_LIMIT = 5
	HISTOGRAM_SCALE = 'log'
	HISTOGRAM_TYPE = 'step'
	INTERVAL = 'zscale'
	PLOT_GALACTIC_CENTER = True
	PLOT_GALACTIC_CONTOUR = True
	PLOT_GALACTIC_PLANE = True
	PROJECTION = 'hammer'
	SAVE_FIGURE = True
	SKYMAP_MODE = 'comm_fields' # comm_fields, survey_fields

	# listener
	CLIENT_ID = '7vjhnkra2bveb93u0d2dqdmtks'
	CLIENT_SECRET = '1v840t3ljqgob81m8kk0903edss4ul4huf5h3tq931f15lsjhoab'
	EMAIL_ID = 'toros.alerts@yahoo.com'
	EMAIL_PASSWORD = 'ncnrlqwdthofhoch'
	EMAIL_PORT = 465
	EMAIL_SMTP = 'smtp.mail.yahoo.com'
	FIELD_PRIORITY = ['commissioning', 'ligo', 'science', 'main_survey']
	GALACTIC_PLANE = 15.
	GALAXY_CATALOG = 'glade+'
	GALAXY_LIST = 50
	GENERATE_FIELDS = 'N'
	LISTEN_NED_WAIT = 1
	MAILING_LIST = ['2674218771@tmomail.net', 'richardcamuccio@aol.com']
	MOON_DISTANCE = 60.
	SEND_ALERT = True