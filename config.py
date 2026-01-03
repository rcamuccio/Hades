from libraries.utils import Utils

class Configuration:

    # computer
    MACHINE = 'epimetheus'
    RAW_FILE_EXTENSION = '.fits'
    FILE_EXTENSION = '.fits'

    # field
    FIELD = 'FIELD_30.000'
    RA = 0.723
    DEC = -29.286

    # transient
    TRANSIENT_LC = 'N'
    TRANSIENT_NAME = 'AT2024xsq'
    TRANSIENT_RA = 10.1958012
    TRANSIENT_DEC = -21.9258309

    # reduction toggles
    CLEAN_SKIP = 'N'
    WRITE_SKY = 'N'
    MASTER_SKIP = 'N'
    DIFFERENCE_SKIP = 'N'
    PHOTOMETRY_SKIP = 'N'
    LIGHTCURVE_SKIP = 'N'

    # pipeline toggles
    SUBTRACT_BIAS = 'Y'
    SUBTRACT_DARK = 'Y'
    DIVIDE_FLAT = 'Y'
    CLIP_IMAGE = 'Y'
    SUBTRACT_SKY = 'Y'
    PLATE_SOLVE = 'Y'
    SOURCE_EXTRACTOR_PATH = '/usr/bin/source-extractor'

    # telescope
    PIXEL_SIZE = 0.4959 # arcsec/px
    NUM_PIXELS = 10560 # px/side
    TOROS_DEC_LIMIT = 26.66 # deg
    FOV = (PIXEL_SIZE * NUM_PIXELS) / 3600. # deg
    SEARCH_DIST = FOV # deg
    EXP_TIME = 300. # s
    DARK_EXP = 300. # s
    FLAT_EXP = 5. # s
    GAIN = 0.380  # e-/ADU
    PEAKMAX = 45000  # peakmax used in DAOStarFinder (use None if unsure)

    # image
    AXS_X_RW = 12000
    AXS_Y_RW = 10600
    AXS_X = 10560
    AXS_Y = 10560
    AXS = 10560
    OVS_X = 180
    OVS_Y = 20
    TAP_X = 1320
    TAP_Y = 5280

    # differencing
    KRNL = 2 # kernel size = 2 * KNRL + 1
    STMP = 10 # stamp size = 2 * STMP + 1
    ORDR = 1 # order of kernel; 0 = stationary, 1/2 = spatially varying
    NRSTARS = 500 #number of stars used to solve for kernel
    BRIGHT_STARS = 20000 # top stars to search for in kernel stars
    KERNEL_LIMIT = 0.5 # maximum allowable offset in zeropoint in magnitudes
    AXS_LIMIT = 100 # number of pixels close to edge of frame
    RMS_LOW_LIMIT = 0.005 # lower limit on precision to use for kernel stars
    RMS_UP_LIMIT = 0.02  # upper limit on precision to use for kernel stars

    # sky subtraction
    PIX = 220

    # photometry
    APER_SIZE = 16 # circular aperture
    ANNULI_INNER = APER_SIZE + 2
    ANNULI_OUTER = APER_SIZE + 4
    FWHM = 15.  # fwhm of the image
    STAR_LIST_MAX = 30000
    THRESHOLD = 5. # the threshold for a source above the background

    # main directory
    MAIN_DIRECTORY = '/home/epimetheus/Documents/'

    # codebase directory
    CODE_DIRECTORY = MAIN_DIRECTORY + 'hades/'

    # output directories
    OUTPUT_DIRECTORY = MAIN_DIRECTORY + 'output/'
    ALERTS_DIRECTORY = OUTPUT_DIRECTORY + 'alerts/'
    ANALYSIS_DIRECTORY = OUTPUT_DIRECTORY + 'analysis/'
    LOG_DIRECTORY = OUTPUT_DIRECTORY + 'logs/'
    CODE_DIFFERENCE_DIRECTORY = OUTPUT_DIRECTORY + 'difference/'

    SECRET_DIRECTORY = '/home/epimetheus/Documents/'

    # data directories
    DATA_DIRECTORY = '/media/epimetheus/Harbor/hades-data/'
    RAW_DIRECTORY = DATA_DIRECTORY + 'raw/'
    CLEAN_DIRECTORY = DATA_DIRECTORY + 'clean/'
    MASTER_MAIN_DIRECTORY = DATA_DIRECTORY + 'master/'
    MASTER_DIRECTORY = MASTER_MAIN_DIRECTORY + FIELD + '/'
    MASTER_TMP_DIRECTORY = MASTER_DIRECTORY + 'tmp_master/'
    CENTROID_DIRECTORY = MASTER_DIRECTORY + 'centroids/'
    CALIBRATION_DIRECTORY = DATA_DIRECTORY + 'calibration/'
    BIAS_DIRECTORY = CALIBRATION_DIRECTORY + 'tmp_bias/'
    FLAT_DIRECTORY = CALIBRATION_DIRECTORY + 'tmp_flat/'
    DARK_DIRECTORY = CALIBRATION_DIRECTORY + 'tmp_dark/'
    LIGHTCURVE_DIRECTORY = DATA_DIRECTORY + 'lc/'
    LIGHTCURVE_FIELD_DIRECTORY = LIGHTCURVE_DIRECTORY + FIELD + "/"
    DIFFERENCED_DIRECTORY = DATA_DIRECTORY + 'diff/'
    FLUX_DIRECTORY = DATA_DIRECTORY + 'flux/'
    REVIEW_DIRECTORY = DATA_DIRECTORY + 'review/'

    # directory list
    DIRECTORIES = [ANALYSIS_DIRECTORY, DATA_DIRECTORY, LOG_DIRECTORY, CALIBRATION_DIRECTORY, FLUX_DIRECTORY, CLEAN_DIRECTORY, MASTER_MAIN_DIRECTORY, MASTER_DIRECTORY, MASTER_TMP_DIRECTORY, LIGHTCURVE_DIRECTORY, CENTROID_DIRECTORY, RAW_DIRECTORY, BIAS_DIRECTORY, DARK_DIRECTORY, FLAT_DIRECTORY, DIFFERENCED_DIRECTORY,CODE_DIFFERENCE_DIRECTORY, LIGHTCURVE_FIELD_DIRECTORY, REVIEW_DIRECTORY]

    # listener
    CLIENT_ID = Utils.get_credentials('id')
    CLIENT_SECRET = Utils.get_credentials('secret')
    EMAIL = 'toros.alerts@yahoo.com'
    MAILING_LIST = ['rcamuccio@gmail.com']
    PAS = 'ncnrlqwdthofhoch'
    PORT = 465
    SMTP = 'smtp.mail.yahoo.com'

    AVAILABLE_TOPICS = ['gcn.circulars', 'gcn.heartbeat', 'gcn.notices.icecube.lvk_nu_track_search', 'igwn.gwalert', 'gcn.notices.swift.bat.guano', 'gcn.notices.einstein_probe.wxt.alert', 'gcn.notices.superk.sn_alert']

    # BROKER CONFIGURATION SPECIFICS
    LISTEN_NED_WAIT = 1

    # observing conditions
    SEEING = 0.93  # assumes 2 pix FWHM

    # sky brightness at TOLAR in SDSS griz
    SKY = [22.1, 21.1, 20.1, 18.7]

    # SDSS griz bandpass values in nm (width of the filter)
    BP = [147, 141, 147, 147]

    # SDSS griz central wavelength
    CWL = [473.5, 638.5, 775.5, 922.5]

    # telescope information
    TOROS_MIRROR_D = 0.610  # m
    TOROS_MIRROR_R = TOROS_MIRROR_D / 2  # nm

    # CCD information
    READOUT_NOISE = 5.0  # electrons / pixel
    CCD_QE = 0.85
    FILTER_QE = 0.9
    TELESCOPE_SEC_QE = 0.96
    TELECSCOPE_PRI_QE = 0.96
    VIGNETTING = 0.756
    ATMOSPHERE_QE = [0.8, 0.9, 0.9, 0.9]

    # total throughput needs to be multiplied by atmospheric quantum efficiency
    TOTAL_THROUGHPUT = CCD_QE * FILTER_QE * TELESCOPE_SEC_QE * TELECSCOPE_PRI_QE * VIGNETTING

    # force toros field generation?
    FIELD_GENERATION = 'N'

    # The Felix Aguilar Observatory is more Southern than Tolar
    TOROS_LONGITUDE = -69.3265  # -67.32833333
    TOROS_LATITUDE = -31.8023  # -24.62055556
    TOROS_ELEVATION = 2420
    UTC = -3
    MOON_DISTANCE = 60

    EXPOSURE_TIME = 300
    EXPOSURE_TIME_DAY = EXPOSURE_TIME / 60. / 60. / 24.
    NUM_EXPOSURES = 1
    READ_TIME = 90
    OVERHEAD = 30
    TOTAL_EXPOSURE = (EXPOSURE_TIME + READ_TIME) * NUM_EXPOSURES + OVERHEAD