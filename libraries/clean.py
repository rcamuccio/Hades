from config import Configuration
from libraries.preprocessing import Preprocessing
from libraries.utils import Utils
import numpy as np
import os
import time
from astropy.io import fits

class Clean:

	@staticmethod
	def clean_images(clip_image='N', subtract_bias='N', subtract_dark='N', divide_flat='N', subtract_sky='N', plate_solve='N'):
		''' This is the main function script to clean multiple images; alternatively clean_img can be used to clean a single image.

		:parameter clip_image - Y/N if you want to clip any excess from the images (default = N)
		:parameter subtract_bias - Y/N if you want to subtract the bias from the images (default = N)
		:parameter subtract_dark - Y/N if you want to subtract a dark frame from the science images (default = N)
		:parameter divide_flat - Y/N if you want to flatten the images (default = N)
		:parameter subtract_sky - Y/N if you want to subtract the sky from the images (default = N)
		:parameter plate_solve - Y/N if you want to plate solve the image (default = N)

		:return - Nothing is returned, but the images from in_path are cleaned and deposited in out_path
		'''

		st = time.time()

		Utils.log('Getting file list...', 'info')
		files, date_dirs = Utils.get_all_files_per_field(Configuration.RAW_DIRECTORY, Configuration.FIELD, Configuration.RAW_FILE_EXTENSION)

		output_dirs = []

		for dte in date_dirs:
			output_dirs.append(Configuration.DATA_DIRECTORY + 'clean/' + dte)
			output_dirs.append(Configuration.DATA_DIRECTORY + 'clean/' + dte + '/' + Configuration.FIELD)
			output_dirs.append(Configuration.DATA_DIRECTORY + 'diff/' + dte)
			output_dirs.append(Configuration.DATA_DIRECTORY + 'diff/' + dte + '/' + Configuration.FIELD)

		Utils.create_directories(output_dirs)

		if len(files) == 0:
			Utils.log('No .fits files found for ' + Configuration.FIELD + '!' + ' Breaking...', 'debug')
			return

		Utils.log('Starting to clean ' + str(len(files)) + ' images.', 'info')
		for idx, file in enumerate(files):

			file_name = Preprocessing.mk_nme(file, 'N', clip_image, subtract_sky, subtract_bias, divide_flat, subtract_dark, plate_solve)

			if os.path.isfile(file_name) == 1:
				Utils.log('Image ' + file_name + ' already exists. Skipping for now...', 'info')

			if os.path.isfile(file_name) == 0:
				clean_img, header, bd_flag = Clean.clean_img(file, clip_image, subtract_bias, subtract_dark, divide_flat, subtract_sky, plate_solve)

				if bd_flag == 0:
					fits.writeto(file_name, clean_img, header, overwrite=True)
					Utils.log('Cleaned image written as ' + file_name + '.', 'info')
				else:
					Utils.log(file_name + ' is a bad image. Not written.', 'info')

			Utils.log(str(len(files) - idx - 1) + ' images remain to be cleaned.', 'info')

		fn = time.time()
		Utils.log('Image cleaning complete in ' + str(np.around((fn - st), decimals=2)) + 's.', 'info')

	@staticmethod
	def clean_img(file, clip_image='Y', subtract_bias='N', subtract_dark='N', divide_flat='N', subtract_sky='N', plate_solve='N'):
		''' This function is the primary script to clean the image. Various other functions found in this class can be found in the various libraries imported.

		:parameter file - The file name of the image you would like to clean
		:parameter clip_image - Y/N if you want to clip the image (default = Y)
		:parameter subtract_sky - Y/N if you want to subtract the sky from the image (default = N)
		:parameter subtract_bias - Y/N if you want to remove a bias frame (default = N)
		:parameter divide_flat - Y/N if you want to flatten the image (default = N)
		:parameter subtract_dark - Y/N if you want to subtract the dark frame (default = N)
		:parameter plate_solve - Y/N if you want to plate solve the image (default = N)
		'''

		Utils.log('Now cleaning ' + file + '.', 'info')

		img, header = fits.getdata(file, header=True)

		if (subtract_bias == 'Y') & (subtract_dark == 'Y'):
			st = time.time()
			bias, dark = Preprocessing.mk_combined_bias_and_dark(image_overwrite='N')
			img, header = Preprocessing.subtract_scaled_bias_dark(img, header)
			fn = time.time()
			Utils.log('Image bias and dark corrected in ' + str(np.around((fn - st), decimals=2)) + 's.', 'info')
		else:
			Utils.log('Skipping bias and darks subtraction...', 'info')

		if divide_flat == 'Y':
			st = time.time()
			flat = Preprocessing.mk_flat(5, 300)
			img, header = Preprocessing.divide_flat(img, header)
			fn = time.time()
			Utils.log('Image flattened in ' + str(np.around((fn - st), decimals=2)) + 's.', 'info')
		else:
			Utils.log('Skipping image flattening...', 'info')

		if clip_image == 'Y':
			st = time.time()
			img, header = Preprocessing.clip_image(img, header)
			fn = time.time()
			Utils.log('Image overscan removed in ' + str(np.around(fn - st, decimals=2)) + 's.', 'info')
		else:
			Utils.log('Skipping overscan removal...', 'info')

		if subtract_sky == 'Y':
			st = time.time()
			Utils.log('A background box of ' + str(Configuration.PIX) + ' x ' + str(Configuration.PIX) + ' will be used for background subtraction.', 'info')
			img, header = Preprocessing.subtract_sky(img, header, Configuration.WRITE_SKY)
			fn = time.time()
			Utils.log('Sky subtracted in ' + str(np.around((fn - st), decimals=2)) + 's.', 'info')
		else:
			Utils.log('Skipping sky subtraction...', 'info')

		if plate_solve == 'Y':
			st = time.time()
			Utils.log('Now plate solving and correcting the header...', 'info')
			img, header = Preprocessing.correct_header(img, header)
			fn = time.time()
			Utils.log('The image has been plate solved in ' + str(np.around((fn - st), decimals=2)) + 's.', 'info')

		Utils.log('Cleaning finished.', 'info')

		bd_flag = 0

		return img, header, bd_flag