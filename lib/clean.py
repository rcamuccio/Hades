from config import Configuration
from lib.preprocessing import Preproc
from lib.utilities import Utils
from astropy.io import fits
import numpy as np
import os
import sys
import time

class Clean:

	@staticmethod
	def clean_images():
		'''This is the main script to clean multiple images.

		:return - Nothing is returned
		'''

		data_dirs = Utils.config_data_dir(create=False)

		# get file list for all dates field was observed
		Utils.log('Getting file list.', 'info')
		files, dates = Utils.get_all_files_per_field(data_dirs['raw'], Configuration.FIELD, 'raw', Configuration.FILE_EXTENSION)

		# make output directories for data files
		output_dirs = []
		for dte in dates:
			output_dirs.append(data_dirs['clean'] + dte)
			output_dirs.append(data_dirs['clean'] + dte + '/' + Configuration.FIELD)
			output_dirs.append(data_dirs['review'] + dte)
			output_dirs.append(data_dirs['review'] + dte + '/' + Configuration.FIELD)
		Utils.create_directories(output_dirs)

		# break if there are no files
		if len(files) == 0:
			Utils.log('No ' + Configuration.FILE_EXTENSION + ' files found for ' + Configuration.FIELD + '! Breaking.', 'debug')
			sys.exit()

		# begin file cleaning
		Utils.log('Starting to clean ' + str(len(files)) + ' images.', 'info')

		for idx, file in enumerate(files):
			# make a new name for the file based on which actions are taken
			file_name = Preproc.mk_nme(file, 'N', Configuration.SUBTRACT_BIAS, Configuration.SUBTRACT_DARK, Configuration.DIVIDE_FLAT, Configuration.CLIP_IMAGE, Configuration.SUBTRACT_SKY, Configuration.PLATE_SOLVE)

			raw_file_name = file.split('/')[-1]
			mod_file_name = file_name.split('/')[-1]

			# check if file exists
			if os.path.isfile(file_name):
				Utils.log('Image ' + mod_file_name + ' already exists. Skipping.', 'info')
			# clean the image
			else:
				clean_img, header = Clean.clean_img(file, file_name)

			remain = len(files) - idx - 1
			if remain == 1:
				Utils.log(str(remain) + ' image remains to be cleaned.', 'info')
			else:
				Utils.log(str(remain) + ' images remain to be cleaned.', 'info')

	@staticmethod
	def clean_img(file, file_name=None):
		'''This function is the primary script to clean an image.

		:parameter file - The file path of the unprocessed image
		:parameter file_name - The file name of the unprocessed image

		:return img - 
		:return header - 
		'''

		raw_file_nme = file.split('/')[-1]
		Utils.log('Cleaning ' + raw_file_nme + '.', 'info')

		# load image
		img, header = fits.getdata(file, header=True)

		# remove bias and dark
		if (Configuration.SUBTRACT_BIAS == 'Y') & (Configuration.SUBTRACT_DARK == 'Y'):
			bias, dark = Preproc.mk_combined_bias_and_dark(image_overwrite='N')
			img, header = Preproc.subtract_scaled_bias_dark(img, header)
		else:
			Utils.log('Skipping bias and dark subtraction.', 'info')

		# flat divide
		if Configuration.DIVIDE_FLAT == 'Y':
			flat = Preproc.mk_flat()
			img, header = Preproc.flat_divide(img, header)
		else:
			Utils.log('Skipping image flattening.', 'info')

		# clip image
		if Configuration.CLIP_IMAGE == 'Y':
			img, header = Preproc.clip_image(img, header)
		else:
			Utils.log('Skipping overscan removal.', 'info')

		# sky subtract
		if Configuration.SUBTRACT_SKY == 'Y':
			img, header = Preproc.sky_subtract(img, header, Configuration.WRITE_SKY)
		else:
			Utils.log('Skipping sky subtraction.', 'info')

		# plate solve
		if Configuration.PLATE_SOLVE == 'Y':
			img, header = Preproc.correct_header(img, header, file_name)

		Utils.log('Cleaning finished on ' + raw_file_nme + '.', 'info')

		return img, header