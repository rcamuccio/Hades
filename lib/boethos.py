'''
Βοηθός - The Helper
'''

import os

from conf import Configuration

class Boethos:

	@staticmethod
	def create_directories(path_list):
		'''This function checks for each directory in the list, and creates it if it does not exist.

		:p path_list [list][str] - the list of directories to create

		:r - nothing
		'''

		preamble = Boethos.get_line_preamble('library', 'boethos', 'create_directories')

		for path in path_list:
			if not os.path.exists(path):
				os.makedirs(path)
				print(preamble, path, 'created')
			else:
				print(preamble, path, 'already exists')

	@staticmethod
	def get_greek(string):

		if string == 'angelos':
			greek = 'Ἄγγελος'

		elif string == 'boethos':
			greek = 'Βοηθός'

		elif string == 'grapheus':
			greek = 'Γραφεύς'

		elif string == 'logistes':
			greek = 'Λογιστής'

		elif string == 'metretes':
			greek = 'Μετρητής'

		elif string == 'skopeus':
			greek = 'Σκοπεύς'

		elif string == 'zetetes':
			greek = 'Ζητητής'

		return greek

	@staticmethod
	def get_line_preamble(type, module, function=None):

		if type == 'library':
			greek = Boethos.get_greek(module)
			preamble = '\033[1m[' + greek + '|' + function + ']\033[0m :'

		elif type == 'script':
			preamble = '\033[1m[' + module + ']\033[0m :'

		return preamble

	@staticmethod
	def get_script_postamble(script):

		if script == 'gcn_listener':
			postamble = '\n\033[1m' + 'ᾍδης ἀπέρχεται' + '\033[0m | Script terminated [gcn_listener]'

		elif script == 'preprocessing_pipeline':
			postamble = '\n\033[1m' + 'ᾍδης ἀπέρχεται' + '\033[0m | Script terminated [preprocessing_pipeline]'

		elif script == 'skymap_plotter':
			postamble = '\n\033[1m' + 'ᾍδης ἀπέρχεται' + '\033[0m | Script terminated [skymap_plotter]'

		return postamble

	@staticmethod
	def get_script_preamble(script):

		if script == 'gcn_listener':
			preamble = '\033[1m' + 'ᾍδης ζῇ' + '\033[0m | Running script [gcn_listener]\n'

		elif script == 'preprocessing_pipeline':
			preamble = '\033[1m' + 'ᾍδης ζῇ' + '\033[0m | Running script [preprocessing_pipeline]\n'

		elif script == 'skymap_plotter':
			preamble = '\033[1m' + 'ᾍδης ζῇ' + '\033[0m | Running script [skymap_plotter]\n'

		return preamble

	@staticmethod
	def setup_plouton():

		preamble = Boethos.get_line_preamble('library', 'boethos', 'setup_plouton')
		print(preamble, 'Setting up Πλούτων')

		directories = {}
		plouton_directory = Configuration.PLOUTON_DIRECTORY

		alert_directory = plouton_directory + 'alert/'
		camera_directory = plouton_directory + 'camera/'
		catalog_directory = plouton_directory + 'catalog/'
		field_directory = plouton_directory + 'field/'
		skymap_directory = plouton_directory + 'skymap/'

		directories['alert'] = alert_directory
		directories['camera'] = camera_directory
		directories['catalog'] = catalog_directory
		directories['field'] = field_directory
		directories['skymap'] = skymap_directory

		directory_list = [alert_directory, camera_directory, catalog_directory, field_directory, skymap_directory]

		Boethos.create_directories(directory_list)