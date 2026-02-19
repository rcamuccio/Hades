class Configuration:

	# computer
	MACHINE = 'epimetheus'
	FILE_EXTENSION = '.fits'
	OBSERVATORY = 'toros'
	FIELD = 'FIELD_0b.001'
	KNOWN_VARIABLES = 'Y'

	# directories
	MAIN_DIRECTORY = '/home/epimetheus/Documents/'
	CODE_DIRECTORY = MAIN_DIRECTORY + 'hades/'
	OUTPUT_DIRECTORY = MAIN_DIRECTORY + 'output/'
	DATA_DIRECTORY = '/media/epimetheus/Harbor/'

	# pipeline toggles
	CLEAN_SKIP = 'N'
	WRITE_SKY = 'N'
	MASTER_SKIP = 'N'
	DIFFERENCE_SKIP = 'N'
	PHOTOMETRY_SKIP = 'N'
	LIGHTCURVE_SKIP = 'N'

	# reduction toggles
	SUBTRACT_BIAS = 'Y'
	SUBTRACT_DARK = 'Y'
	DIVIDE_FLAT = 'Y'
	CLIP_IMAGE = 'Y'
	SUBTRACT_SKY = 'Y'
	PLATE_SOLVE = 'Y'
	SOURCE_EXTRACTOR_PATH = '/usr/bin/source-extractor'

	# photometry
	AIRMASS_METHOD = 'ky1994'
	ANNULI_INNER = 18
	ANNULI_OUTER = 20
	APER_SIZE = 16 # circular aperture
	AXS_LIMIT = 100 # number of pixels close to edge of frame
	BRIGHT_STARS = 20000 # top stars to search for in kernel stars
	FWHM = 15. # fwhm of the image
	KERNEL_LIMIT = 0.5 # maximum allowable offset in zeropoint in magnitudes
	KRNL = 2 # kernel size = 2 * KRNL + 1
	NRSTARS = 500 # number of stars used to solve for kernel
	ORDR = 1 # order of kernel; 0 = stationary, 1/2 = spatially varying
	PIX = 220 # sky subtraction
	RMS_LOW_LIMIT = 0.005 # lower limit on precision to use for kernel stars
	RMS_UPP_LIMIT = 0.02 # upper limit on precision to use for kernel stars
	STMP = 10 # stamp size = 2 * STMP + 1
	THRESHOLD = 5. # threshold for a source above the background

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
	QUERY = 'gaia-cone'
	SMTP = 'smtp.mail.yahoo.com'

	AVAILABLE_TOPICS = ['gcn.circulars', 'gcn.heartbeat', 'gcn.notices.icecube.lvk_nu_track_search', 'igwn.gwalert', 'gcn.notices.swift.bat.guano', 'gcn.notices.einstein_probe.wxt.alert', 'gcn.notices.superk.sn_alert']