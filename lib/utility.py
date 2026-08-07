from config import Configuration
import logging
import matplotlib
import os
import warnings

logging.getLogger('ccdproc').setLevel(logging.WARNING)
logging.getLogger('gcn').setLevel(logging.WARNING)
logging.getLogger('healpy').setLevel(logging.WARNING)
logging.getLogger('requests').setLevel(logging.WARNING)
logging.getLogger('urllib3').setLevel(logging.WARNING)
matplotlib.set_loglevel(level='warning')
pil_logger = logging.getLogger('PIL')
pil_logger.setLevel(logging.INFO)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=Warning)

class Utility:

	@staticmethod
	def create_directories(path_list):
		'''This function checks for each directory in the list, and creates it, if it doesn't already exist.

		:parameter path_list - The list of directories to create

		:return - Nothing is returned, but directories are created, if necessary
		'''

		for path in path_list:
			if not os.path.exists(path):
				os.mkdir(path)
				Utility.log(path + ' created.', 'info')
			else:
				Utility.log(path + ' already exists. Skipping.', 'info')

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
		#for file in files_no_ext:
		#	dte_dir.append(file.split('/')[-3])
		#uni_dte_dir = np.unique(dte_dir).tolist()

		# get the unique dates with images
		for file in files_no_ext:
			dte_idx = file.split('/').index(field)
			dte_dir.append(file.split('/')[dte_idx - 1])
		uni_dte_dir = np.unique(dte_dir).tolist()

		file_srt = np.argsort(files_no_ext)

		'''
		Utility.log(str(field) + ' not observed on:', 'info')
		for dte in no_dte_list:
			Utility.log('  ' + str(dte), 'info')

		Utility.log(str(field) + ' observed on:', 'info')
		for dte in uni_dte_dir:
			Utility.log('  ' + str(dte), 'info')
		'''

		return np.array(files_no_ext)[file_srt], uni_dte_dir

	@staticmethod
	def get_directories(location, create=True):

		dirs = {}

		if location == 'data':
			directory = Configuration.DATA_DIRECTORY

			raw_directory = directory + 'raw/'
			clean_directory = directory + 'clean/'
			master_main_directory = directory + 'master/'
			master_directory = master_main_directory + Configuration.FIELD + '/'
			master_tmp_directory = master_directory + 'tmp_master/'
			centroid_directory = master_directory + 'centroids/'
			calibration_directory = directory + 'calibration/'
			bias_directory = calibration_directory + 'tmp_bias/'
			flat_directory = calibration_directory + 'tmp_flat/'
			dark_directory = calibration_directory + 'tmp_dark/'
			lightcurve_directory = directory + 'lc/'
			lightcurve_field_directory = lightcurve_directory + Configuration.FIELD + '/'
			lightcurve_field_raw_directory = lightcurve_directory + Configuration.FIELD + '/raw/'
			lightcurve_field_detrend_directory = lightcurve_directory + Configuration.FIELD + '/detrend/'
			differenced_directory = directory + 'diff/'
			flux_directory = directory + 'flux/'
			review_directory = directory + 'review/'

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
			dirs['lightcurve_field_raw'] = lightcurve_field_raw_directory
			dirs['lightcurve_field_detrend'] = lightcurve_field_detrend_directory
			dirs['differenced'] = differenced_directory
			dirs['flux'] = flux_directory
			dirs['review'] = review_directory

			dir_list = [raw_directory, clean_directory, master_main_directory, master_directory, master_tmp_directory, centroid_directory, calibration_directory, bias_directory, flat_directory, dark_directory, lightcurve_directory, lightcurve_field_directory, lightcurve_field_raw_directory, lightcurve_field_detrend_directory, differenced_directory, flux_directory, review_directory]

		elif location == 'output':
			directory = Configuration.OUTPUT_DIRECTORY

			dir_alerts = directory + 'alerts/'
			dir_analysis = directory + 'analysis/'
			dir_catalogs = directory + 'catalogs/'
			dir_difference = directory + 'difference/'
			dir_logs = directory + 'logs/'
			dir_sources = directory + 'sources/'

			dirs['alerts'] = dir_alerts
			dirs['analysis'] = dir_analysis
			dirs['catalogs'] = dir_catalogs
			dirs['difference'] = dir_difference
			dirs['logs'] = dir_logs
			dirs['sources'] = dir_sources

			dir_list = [dir_alerts, dir_analysis, dir_catalogs, dir_difference, dir_logs, dir_sources]

		if create:
			Utility.create_directories(dir_list)

		return dirs

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

		dirs = Utility.get_directories(location='output', create=False)

		# create the log
		log_path = dirs['logs'] + 'toros.log'
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