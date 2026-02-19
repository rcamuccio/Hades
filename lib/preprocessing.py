from config import Configuration
from lib.utilities import Utils
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.wcs.utils import proj_plane_pixel_scales
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder
import astropy.wcs as pywcs
import numpy as np
import os
import pandas as pd
import scipy.ndimage
import shutil
import subprocess

class Preproc:

	@staticmethod
	def mk_mask(img_data):

		mask = np.zeros(img)

		return mask, boxes

	@staticmethod
	def align_img(image, header1, header2, preserve_bad_pixels=False):
		'''This function is based on the FITS_tools utility hcongrid and interpolates an image from one FITS header to another.

		:parameter image - The image to transform
		:parameter header1 - The header of the image
		:parameter header2 - The header to transform to
		:parameter preserve_bad_pixels - Try to set NAN pixels to NAN in the zoomed image (otherwise bad pixels are set to zero)

		:return new_image - The image which has been transformed to the reference header
		'''

		# convert the headers to WCS objects
		wcs1 = pywcs.WCS(header1)
		wcs1.naxis1 = header1['NAXIS1']
		wcs1.naxis2 = header1['NAXIS2']
		wcs2 = pywcs.WCS(header2)
		wcs2.naxis1 = header2['NAXIS1']
		wcs2.naxis2 = header2['NAXIS2']

		# get the shape
		outshape = [wcs2.naxis2, wcs2.naxis1]
		yy2, xx2 = np.indices(outshape)

		# get the world coordinates of the output image
		lon2, lat2 = wcs2.wcs_pix2world(xx2, yy2, 0)
		xx1, yy1 = wcs1.wcs_world2pix(lon2, lat2, 0)

		# make a grid for the image to transform to
		grid1 = np.array([yy1.reshape(outshape), xx1.reshape(outshape)])

		# identify bad pixels
		bad_pixels = np.isnan(image) + np.isinf(image)
		image[bad_pixels] = 0

		# make the new image
		new_image = scipy.ndimage.map_coordinates(image, grid1)

		# replace bad pixels
		if preserve_bad_pixels:
			newbad = scipy.ndimage.map_coordinates(bad_pixels, grid1, order=0, mode='constant', cval=np.nan)
			new_image[newbad] = np.nan

		return new_image

	@staticmethod
	def clip_image(img, header):
		'''This function will clip the image and remove the overscan region. This function is written for TOROS specifically, and will need to be updated for any given CCD.

		:parameter img - The image to clip
		:parameter header - The header of the image

		:return image_clip - The clipped image
		:return header - The updated header
		'''

		observatory = Utils.config_observatory(Configuration.OBSERVATORY)

		# make the clipped image
		image_clip = np.zeros((observatory['axs_y'], observatory['axs_x']))

		# size of the overscan
		ovs_x_size = observatory['ovs_x']
		ovs_y_size = observatory['ovs_y']

		# size of the taps
		chip_x_size = observatory['tap_x']
		chip_y_size = observatory['tap_y']

		# size of the full chip
		full_chip_x = chip_x_size + ovs_x_size
		full_chip_y = chip_y_size + ovs_y_size

		# move through x and y
		idx = 0
		for x in range(0, observatory['axs_x_rw'], full_chip_x):
			idy = 0
			for y in range(0, observatory['axs_y_rw'], full_chip_y):
				# put the clipped image into the holder image
				if y == 0:
					image_clip[idy:idy + chip_y_size, idx:idx + chip_x_size] = img[y:y + chip_y_size, x:x + chip_x_size]
				else:
					image_clip[idy:idy + chip_y_size, idx:idx + chip_x_size] = img[y + ovs_y_size:y + ovs_y_size + chip_y_size, x:x + chip_x_size]

				# increase the size of the yclip
				idy += chip_y_size

			# increase the size of the xclip
			idx += chip_x_size

		# update the header
		header['OVERSCAN'] = 'removed'
		header['X_CLIP'] = ovs_x_size
		header['Y_CLIP'] = ovs_y_size

		return image_clip, header

	@staticmethod
	def correct_header(img, header, clean_file_name):
		'''This function will plate solve the image, add the time stamp, and exposure time to the header if need be.

		:parameter img - This is the image you want to plate solve
		:parameter header - The header of the image you want to correct
		:parameter clean_file_name - 

		:return img - 
		:return header - 
		'''

		output_dirs = Utils.config_output_dir(create=False)

		working_directory_file = output_dirs['analysis'] + os.path.basename(clean_file_name) + '.temp'
		fits.writeto(working_directory_file, img, header, overwrite=True)

		output_file = os.path.splitext(working_directory_file)[0]
		output_dir = os.path.dirname(working_directory_file)

		astrometry_command = 'solve-field --use-source-extractor --source-extractor-path ' + Configuration.SOURCE_EXTRACTOR_PATH + ' --scale-units arcsecperpix --scale-low 0.48 --scale-high 0.5 --no-plots --temp-axy --index-xyls none --match none --rdls none --solved none --corr none --dir ' + output_dir + ' --new-fits ' + output_file + ' ' + working_directory_file

		subprocess.run(astrometry_command, shell=True)

		# clean up astrometry output, remove .temp and .wcs files
		astrometry_extensions = ['.temp', '.wcs']

		# check that astrometry generated plate solved file (.fits)
		if os.path.exists(output_file):
			shutil.move(output_file, clean_file_name)
			Utils.log(f'Cleaned image written as ' + clean_file_name + '.', 'info')
			# if file exists then remove extra files
			for extension in astrometry_extensions:
				file_to_remove = output_file + extension
				if os.path.exists(file_to_remove):
					os.remove(file_to_remove)
					Utils.log(f'File {file_to_remove} deleted.', 'info')
		else:
			# if file does not exist then move temp file from clean for review
			temp_file_name = clean_file_name.replace('/clean/', '/review/')
			shutil.move(working_directory_file, temp_file_name)
			Utils.log(f'Plate solve failed. Moving file to review directory for human inspection.', 'info') 

		return img, header

	@staticmethod
	def flat_divide(img, header):
		'''This function will divide an image with a flatfield.

		:parameter img - The image to flatten
		:parameter header - The image header

		:return flat_div, header - The flattened image
		:return header - The updated header
		'''

		data_dirs = Utils.config_data_dir(create=False)

		# read in the flat frame
		flat = fits.getdata(data_dirs['calibration'] + 'flat' + Configuration.FILE_EXTENSION)

		# subtract the bias from the image
		flat_div = img / flat

		# update the header
		header['FLATTEN'] = 'Y'

		return flat_div, header

	@staticmethod
	def mk_combined_bias_and_dark(image_overwrite):
		'''This function will make a master bias frame (sans a mean pixel value) and a master dark frame. The master dark frame uses the bias frames from the same night to generate. The master bias frame has its median pixel value removed so it can be scaled to the value of the overscan region later. This is because the TOROS bias level seems to be changing with time.

		:parameter image_overwrite - Y/N if you want to force an overwrite for each master frame

		:return bias - The master bias frame
		:return dark - The master dark frame
		'''

		observatory = Utils.config_observatory(Configuration.OBSERVATORY)
		data_dirs = Utils.config_data_dir(create=False)

		# set the names of the master frames
		bias_name = 'bias' + Configuration.FILE_EXTENSION
		dark_name = 'dark' + Configuration.FILE_EXTENSION

		# check to see if the master bias frame exists
		if ((os.path.isfile(data_dirs['calibration'] + dark_name) == 0) or (os.path.isfile(data_dirs['calibration'] + bias_name) == 0) or (image_overwrite == 'Y')):

			# get the image lists
			biases = pd.read_csv(data_dirs['calibration'] + 'bias_list.csv', sep=',')
			darks = pd.read_csv(data_dirs['calibration'] + 'dark_list.csv', sep=',')
			total_bias = 0
			total_dark = 0

			# group based on the date to get the unique dates
			bias_dates = biases.groupby('Date').agg({'count'}).reset_index()['Date'].to_list()
			ndates = len(bias_dates)

			# make holders for the total number of nights
			bias_bulk = np.ndarray(shape=(ndates, observatory['axs_y_rw'], observatory['axs_x_rw']))
			dark_bulk = np.ndarray(shape=(ndates, observatory['axs_y_rw'], observatory['axs_x_rw']))

			bias_bulk_filepath = []
			dark_bulk_filepath = []
			zdx = 0
			for dte in bias_dates:
				# determine how many bias and dark frames exist on this date
				bias_frames = biases[biases.Date == dte]
				bias_list = bias_frames.apply(lambda x: Configuration.DATA_DIRECTORY + 'bias/'  + x.Date + x.File, axis=1).to_list()
				nbias = len(bias_frames)

				dark_frames = darks[darks.Date == dte]
				ndarks = len(dark_frames)

				# if dark frames exist, then move forward
				if ndarks > 0:
					Utils.log('Making mini bias and dark files for ' + dte + '.', 'info')
					dark_list = dark_frames.apply(lambda x: Configuration.DATA_DIRECTORY + 'darks/' + x.Date + x.Files, axis=1).to_list()

					# generate the frame holder
					bias_hld = np.ndarray(shape=(nbias, observatory['axs_y_rw'], observatory['axs_x_rw']))
					dark_hld = np.ndarray(shape=(ndarks, observatory['axs_y_rw'], observatory['axs_x_rw']))

					# read in the bias file
					idx = 0
					for ii in range(0, nbias):
						bias_tmp, bias_head = fits.getdata(bias_list[ii], header=True)
						bias_hld[idx] = bias_tmp
						idx += 1
						total_bias += 1
						bias_tmp = None
					bias_tmp_mdn = np.median(bias_hld, axis=0)
					bias_hld = None
					Utils.log('Bias done.', 'info')

					# read in the dark file
					jdx = 0
					for jj in range(0, ndarks):
						dark_tmp, dark_head = fits.getdata(dark_list[jj], header=True)
						dark_hld[jdx] = dark_tmp - bias_tmp_mdn
						jdx += 1
						total_dark += 1
						dark_tmp = None
					dark_tmp_mdn = np.median(dark_hld, axis=0)
					dark_hld = None
					dark_tmp_filename = data_dirs['dark'] + str(zdx).zfill(2) + '_scale_tmp_dark.fits'
					fits.writeto(dark_tmp_filename, dark_tmp_mdn, overwrite=True)
					Utils.log('Dark done.', 'info')

					# size of the overscan
					ovs_x_size = observatory['ovs_x']
					ovs_y_size = observatory['ovs_y']

					# size of the taps
					chip_x_size = observatory['tap_x']
					chip_y_size = observatory['tap_y']

					# size of the full chip
					full_chip_x = chip_x_size + ovs_x_size
					full_chip_y = chip_y_size + ovs_y_size

					# move through x, y to mask the 'image' parts of the image
					for x in range(0, observatory['axs_x_rw'], full_chip_x):
						for y in range(0, observatory['axs_y_rw'], full_chip_y):
							# put the clipped image into the holder image
							bias_tmp_mdn[y:y + full_chip_y, x:x + full_chip_x] = (
									bias_tmp_mdn[y:y + full_chip_y, x:x + full_chip_x] -
									np.median(bias_tmp_mdn[y:y + full_chip_y, x:x + full_chip_x]))
					bias_tmp_filename = data_dirs['bias'] + str(zdx).zfill(2) + '_scale_tmp_bias.fits'
					fits.writeto(bias_tmp_filename, bias_tmp_mdn, overwrite=True)

					# update the bulk holder
					bias_bulk_filepath.append(bias_tmp_filename)
					dark_bulk_filepath.append(dark_tmp_filename)
					zdx +=1 
					bias_tmp_mdn = None
					dark_tmp_mdn = None

			# load bias bulk
			for index, tmp_bias_file in enumerate(bias_bulk_filepath):
				bias_bulk[index] = fits.getdata(tmp_bias_file)

			# update the header with relevant information
			bias_hdu = fits.PrimaryHDU()
			bias_header = bias_hdu.header
			bias_header['BIAS_COMB'] = 'median'
			bias_header['NUM_BIAS'] = total_bias
			bias = np.median(bias_bulk[0:zdx], axis=0)
			bias_bulk = None
			fits.writeto(data_dirs['calibration'] + bias_name, bias, bias_header, overwrite=True)

			# load dark bulk
			for index, tmp_dark_file in enumerate(dark_bulk_filepath):
				dark_bulk[index] = fits.getdata(tmp_dark_file)

			# update the header with relevant information
			dark_hdu = fits.PrimaryHDU()
			dark_header = dark_hdu.header
			dark_header['DARK_COMB'] = 'median'
			dark_header['NUM_DARK'] = total_dark
			dark_header['EXPTIME'] = observatory['exp_dark']
			dark_header['BIAS_SUB'] = 'Y'
			dark = np.median(dark_bulk[0:zdx], axis=0)
			dark_bulk = None
			# write the image out to the master directory
			fits.writeto(data_dirs['calibration'] + dark_name, dark, dark_header, overwrite=True)

			Utils.log('Scalable bias and dark frames generated.', 'info')

		else:
			bias = fits.getdata(data_dirs['calibration'] + 'bias' + Configuration.FILE_EXTENSION)
			dark = fits.getdata(data_dirs['calibration'] + 'dark' + Configuration.FILE_EXTENSION)

		return bias, dark

	@staticmethod
	def mk_flat():
		'''This function will make the master flat frame using the provided image list.

		:return nflat_image - The flat field for the given date is returned and written to the calibration directory
		'''

		observatory = Utils.config_observatory(Configuration.OBSERVATORY)
		data_dirs = Utils.config_data_dir(create=False)

		if not os.path.isfile(data_dirs['calibration'] + 'flat' + Configuration.FILE_EXTENSION):

			# read in the bias and dark frames
			bias = fits.getdata(data_dirs['calibration'] + 'bias' + Configuration.FILE_EXTENSION)
			dark = fits.getdata(data_dirs['calibration'] + 'dark' + Configuration.FILE_EXTENSION)

			# size of the overscan
			ovs_x_size = observatory['ovs_x']
			ovs_y_size = observatory['ovs_y']

			# size of the taps
			chip_x_size = observatory['tap_x']
			chip_y_size = observatory['tap_y']

			# size of the full chip
			full_chip_x = chip_x_size + ovs_x_size
			full_chip_y = chip_y_size + ovs_y_size

			# get the image list
			images = pd.read_csv(data_dirs['calibration'] + 'flat_list.csv', sep=',')
			image_list = images.apply(lambda x: Configuration.DATA_DIRECTORY + 'flats/' + x.Date + x.Files, axis=1).to_list()

			# determine number of loops needed to move through for each image
			nfiles = len(image_list)
			nbulk = 20

			# get the integer and remainder for the combination
			full_bulk = nfiles // nbulk
			part_bulk = nfiles % nbulk

			if part_bulk > 0:
				hold_bulk = full_bulk + 1
			else:
				hold_bulk = full_bulk

			# here is the 'holder'
			hold_data = np.ndarray(shape=(hold_bulk, observatory['axs_y_rw'], observatory['axs_x_rw']))
			hold_data_filepath = []

			# update the log
			Utils.log('Generating a master flat field from multiple files in bulks of ' + str(nbulk) + ' images. There are ' + str(nfiles) + ' images to combine, which means there should be ' + str(hold_bulk) + ' mini-files to combine.', 'info')

			tmp_num = 0
			for kk in range(0, hold_bulk):
				# loop through the images in sets of nbulk
				if kk < full_bulk:
					# generate the image holder
					block_hold = np.ndarray(shape=(nbulk, observatory['axs_y_rw'], observatory['axs_x_rw']))
					# generate the max index
					mx_index = nbulk
				else:
					# generate the image holder
					block_hold = np.ndarray(shape=(part_bulk, observatory['axs_y_rw'], observatory['axs_x_rw']))
					# generate the max index
					mx_index = part_bulk

				# make the starting index
				loop_start = kk * nbulk
				idx_cnt = 0

				Utils.log('Making mini-flat frame number ' + str(kk) + '.', 'info')

				# loop through the images
				for jj in range(loop_start, mx_index + loop_start):
					# make a copy of the bias frame
					bias_scl = bias.copy()

					# read in the flat file
					flat_tmp, flat_head = fits.getdata(image_list[jj], header=True)

					# get the scale factor for the dark frame
					dark_scale = observatory['exp_dark'] / observatory['exp_flat']

					# move through x, y to mask the 'image' parts of the image
					for x in range(0, observatory['axs_x_rw'], full_chip_x):
						for y in range(0, observatory['axs_y_rw'], full_chip_y):
							# pull out overscan from raw image and set the 'image' to 0
							if y == 0:
								img_slice = flat_tmp[y:y + full_chip_y, x:x + full_chip_x].copy()
								img_slice[0:chip_y_size, 0:chip_x_size] = 0
							else:
								img_slice = flat_tmp[y:y + full_chip_y, x:x + full_chip_x].copy()
								img_slice[ovs_y_size:ovs_y_size + chip_y_size, 0:chip_x_size] = 0
							# put the clipped image into the holder image
							bias_scl[y:y + full_chip_y, x:x + full_chip_x] = (bias_scl[y:y + full_chip_y, x:x + full_chip_x] + np.median(img_slice[img_slice > 0]))
							img_slice = None

					# remove the bias from the image
					flat_bias = flat_tmp - bias_scl

					# remove the dark current from the image
					block_hold[idx_cnt] = flat_bias - (dark / dark_scale)

					if (jj % 20) == 0:
						Utils.log(str(jj + 1) + ' files read in. ' + str(nfiles - jj - 1) + ' files remain.', 'info')

					# increase the iteration
					idx_cnt += 1
					flat_tmp = None
					flat_bias = None
					bias_scl = None

				# median combine the data into a single file
				block_hold_median = np.median(block_hold, axis=0)
				block_hold = None

				# write out the temporary file
				flat_tmp_filename = data_dirs['flat'] + str(tmp_num).zfill(2) + '_tmp_flat' + Configuration.FILE_EXTENSION
				fits.writeto(flat_tmp_filename, block_hold_median, overwrite=True)
				hold_data_filepath.append(flat_tmp_filename)
				block_hold_median = None
				tmp_num += 1

			# load temporary files into hold_data
			for index, tmp_file in enumerate(hold_data_filepath):
				hold_data[index] = fits.getdata(tmp_file)

			# median combine the mini-images into one large image
			flat_image = np.median(hold_data, axis=0)
			nflat_image = flat_image / np.median(flat_image[5950:5950 + 3900, 3200:3200 + 900])

			# pull the header information from the first file of the set
			flat_header = fits.getheader(image_list[0])
			flat_header['comb_typ'] = 'median'
			flat_header['median_val'] = np.median(flat_image[5950:5950 + 3900, 3200:3200 + 900])
			flat_header['norm_pix'] = 'median'
			flat_header['num_comb'] = len(image_list)
			flat_header['mean_pix'] = np.mean(nflat_image)
			flat_header['std_pix'] = np.std(nflat_image)
			flat_header['max_pix'] = np.max(nflat_image)
			flat_header['min_pix'] = np.min(nflat_image)
			flat_header['mean'] = np.mean(flat_image)
			flat_header['std'] = np.std(flat_image)
			flat_header['max'] = np.max(flat_image)
			flat_header['min'] = np.min(flat_image)

			# write the image out to the master directory
			fits.writeto(data_dirs['calibration'] + 'flat' + Configuration.FILE_EXTENSION, nflat_image, flat_header, overwrite=True)

		else:
			nflat_image = fits.getdata(data_dirs['calibration'] + 'flat' + Configuration.FILE_EXTENSION)

		return nflat_image

	@staticmethod
	def mk_nme(file, difference_image='N', bias_subtract='N', dark_subtract='N', flat_divide='N', image_clip='N', sky_subtract='N',  plate_solve='N'):
		'''This function will create the appropriate name for the file based on while steps are taken.

		:parameter file - The string with the file name
		:parameter difference_image - Y/N for image subtraction
		:parameter bias_subtract - Y/N for bias subtraction
		:parameter dark_subtract - Y/N for dark subtraction
		:parameter flat_divide - Y/N for flat division
		:parameter image_clip - Y/N for image clipping
		:parameter sky_subtract - Y/N for sky subtraction
		:parameter plate_solve - Y/N for plate solution

		:return file_name - A string with the new file name
		'''

		data_dirs = Utils.config_data_dir(create=False)

		# if everything is N then the file name is the original filename
		file_name = file

		# update the file name with a 'd' if at the differencing step
		if difference_image == 'Y':
			file = file_name.replace('/clean/', '/diff/')
			nme_hld = file.split(Configuration.FILE_EXTENSION)
			file_name = nme_hld[0] + 'ad' + Configuration.FILE_EXTENSION

		# otherwise...
		if difference_image == 'N':
			# first replace the 'raw' directory with the 'clean' directory
			file_hld = file_name.split('/')
			file = data_dirs['clean'] + file_hld[-3] + '/' + Configuration.FIELD + '/' + file_hld[-1]
			nme_hld = file.split('.fits')

			# update the name to be appropriate for what was done to the file
			# nothing occurs
			if (bias_subtract == 'N') & (flat_divide == 'N') & (sky_subtract == 'N') & (dark_subtract == 'N') & (image_clip == 'N'):
				file_name =  nme_hld[0] + Configuration.FILE_EXTENSION
			# bias
			if (bias_subtract == 'Y') & (flat_divide == 'N') & (sky_subtract == 'N') & (dark_subtract == 'N') & (image_clip =='N'):
				file_name =  nme_hld[0] + '_b' + Configuration.FILE_EXTENSION
			# flat
			if (bias_subtract == 'N') & (flat_divide == 'Y') & (sky_subtract == 'N') & (dark_subtract == 'N') & (image_clip =='N'):
				file_name =  nme_hld[0] + '_f' + Configuration.FILE_EXTENSION
			# sky
			if (bias_subtract == 'N') & (flat_divide == 'N') & (sky_subtract == 'Y') & (dark_subtract == 'N') & (image_clip =='N'):
				file_name =  nme_hld[0] + '_s' + Configuration.FILE_EXTENSION
			# dark
			if (bias_subtract == 'N') & (flat_divide == 'N') & (sky_subtract 	== 'N') & (dark_subtract == 'Y') & (image_clip =='N'):
				file_name =  nme_hld[0] + '_k' + Configuration.FILE_EXTENSION
			# clip
			if (bias_subtract == 'N') & (flat_divide == 'N') & (sky_subtract == 'N') & (dark_subtract == 'N') & (image_clip =='Y'):
				file_name =  nme_hld[0] + '_c' + Configuration.FILE_EXTENSION
			# bias/clip
			if (bias_subtract == 'Y') & (flat_divide == 'N') & (sky_subtract == 'N') & (dark_subtract == 'N') & (image_clip =='Y'):
				file_name =  nme_hld[0] + '_bc' + Configuration.FILE_EXTENSION
			# bias/flat
			if (bias_subtract == 'Y') & (flat_divide == 'Y') & (sky_subtract == 'N') & (dark_subtract == 'N') & (image_clip =='N'):
				file_name =  nme_hld[0] + '_bf' + Configuration.FILE_EXTENSION
			# bias/sky
			if (bias_subtract == 'Y') & (flat_divide == 'N') & (sky_subtract == 'Y') & (dark_subtract == 'N') & (image_clip =='N'):
				file_name =  nme_hld[0] + '_bs' + Configuration.FILE_EXTENSION
			# bias/dark
			if (bias_subtract == 'Y') & (flat_divide == 'N') & (sky_subtract == 'Y') & (dark_subtract == 'N') & (image_clip =='N'):
				file_name =  nme_hld[0] + '_bk' + Configuration.FILE_EXTENSION
			# bias/flat/sky
			if (bias_subtract == 'Y') & (flat_divide == 'Y') & (sky_subtract == 'Y') & (dark_subtract == 'N') & (image_clip =='N'):
				file_name =  nme_hld[0] + '_bfs' + Configuration.FILE_EXTENSION
			# bias/flat/dark
			if (bias_subtract == 'Y') & (flat_divide == 'Y') & (sky_subtract == 'N') & (dark_subtract == 'Y') & (image_clip =='N'):
				file_name =  nme_hld[0] + '_bkf' + Configuration.FILE_EXTENSION
			# bias/flat/clip
			if (bias_subtract == 'Y') & (flat_divide == 'Y') & (sky_subtract == 'N') & (dark_subtract == 'N') & (image_clip =='Y'):
				file_name =  nme_hld[0] + '_bfc' + Configuration.FILE_EXTENSION
			# bias/flat/sky/dark
			if (bias_subtract == 'Y') & (flat_divide == 'Y') & (sky_subtract == 'Y') & (dark_subtract == 'Y') & (image_clip =='N'):
				file_name =  nme_hld[0] + '_bkfs' + Configuration.FILE_EXTENSION
			# bias/flat/sky/clip
			if (bias_subtract == 'Y') & (flat_divide == 'Y') & (sky_subtract == 'Y') & (dark_subtract == 'N') & (image_clip =='Y'):
				file_name =  nme_hld[0] + '_bfcs' + Configuration.FILE_EXTENSION
			# bias/flat/dark/sky/clip
			if (bias_subtract == 'Y') & (flat_divide == 'Y') & (sky_subtract == 'Y') & (dark_subtract == 'Y') & (image_clip =='Y'):
				file_name =  nme_hld[0] + '_bkfcs' + Configuration.FILE_EXTENSION
			# bias/flat/dark/sky/clip/plate
			if (bias_subtract == 'Y') & (flat_divide == 'Y') & (sky_subtract == 'Y') & (dark_subtract == 'Y') & (image_clip =='Y') & (plate_solve == 'Y'):
				file_name =  nme_hld[0] + '_bkfcsp' + Configuration.FILE_EXTENSION

		return file_name

	@staticmethod
	def sky_subtract(img, header, sky_write='N'):
		'''The function has been updated to include the photutils background subtraction routine.

		:parameter img - The image to be cleaned
		:parameter header - The header object to be updated
		:parameter sky_write - Y/N if you want to write the residual background for de-bugging

		:return fin_img - The cleaned image
		:return header - The updated header file
		'''

		observatory = Utils.config_observatory(Configuration.OBSERVATORY)
		output_dirs = Utils.config_output_dir(create=False)

		# set up the background clipping functions
		sigma_clip = SigmaClip(sigma=3)
		bkg_estimator = MedianBackground()

		mask_img = None
		# do the 2D background estimation, if there is no mask, then remove mask_img
		bkg = Background2D(img, (Configuration.PIX, Configuration.PIX), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator, mask=mask_img)
		sky = bkg.background

		# subtract the sky gradient and add back the median background
		img_sub = img - sky
		fin_img = img_sub + bkg.background_median

		# now correct the pixel line near the overscan
		x_skip = observatory['tap_x']
		y_skip = 2 * observatory['ovs_y'] + observatory['ovs_x'] # 220
		y_full = observatory['tap_y']
		for x in range(0, observatory['axs_x'], x_skip):
			for y in range(0, observatory['axs_y'], y_skip):

				# get the sky from the slice and full y column
				mdn_slc = np.median(fin_img[y:y + y_skip, x])
				mdn = np.median(fin_img[y:y + y_full, x])

				# update the column with the new background, unless something chonky is thowing off the calculation
				if mdn_slc < (mdn + 0.1 * mdn):
					fin_img[y:y + y_skip, x] = (fin_img[y:y + y_skip, x] - np.median(fin_img[y:y + y_skip, x]) + bkg.background_median)
				else:
					fin_img[y:y + y_skip, x] = (fin_img[y:y + y_skip, x] - np.median(fin_img[y:y + y_full, x]) + bkg.background_median)

		# update the header
		header['sky'] = bkg.background_median
		header['sky_sig'] = bkg.background_rms_median
		header['sky_sub'] = 'yes'

		# if desired, write out the sky background to the working directory
		if sky_write == 'Y':
			fits.writeto(output_dirs['analysis'] + 'sky_background' + Configuration.FILE_EXTENSION, sky, overwrite=True)
			fits.writeto(output_dirs['analysis'] + 'img' + Configuration.FILE_EXTENSION, img, overwrite=True)
			fits.writeto(output_dirs['analysis'] + 'img_sub' + Configuration.FILE_EXTENSION, fin_img, header=header, overwrite=True)

		return fin_img, header

	@staticmethod
	def subtract_scaled_bias_dark(img, header):
		'''This function will scale the bias frame to match the overscan and then remove the dark level.

		:parameter img - The image to bias and dark correct
		:parameter header - The header of the image

		:return img_bias_dark - The corrected image
		:return header - The updated header
		'''

		observatory = Utils.config_observatory(Configuration.OBSERVATORY)
		data_dirs = Utils.config_data_dir(create=False)

		# read in the bias frame and dark frame
		bias = fits.getdata(data_dirs['calibration'] + 'bias' + Configuration.FILE_EXTENSION)
		dark = fits.getdata(data_dirs['calibration'] + 'dark' + Configuration.FILE_EXTENSION)

		# size of the overscan
		ovs_x_size = observatory['ovs_x']
		ovs_y_size = observatory['ovs_y']

		# size of the taps
		chip_x_size = observatory['tap_x']
		chip_y_size = observatory['tap_y']

		# size of the full chip
		full_chip_x = chip_x_size + ovs_x_size
		full_chip_y = chip_y_size + ovs_y_size

		# make a copy of the bias frame
		bias_scl = bias.copy()
		del bias

		# move through x and y to mask the 'image' parts of the image
		for x in range(0, observatory['axs_x_rw'], full_chip_x):
			for y in range(0, observatory['axs_y_rw'], full_chip_y):
				# pull out the overscan from the raw image and set the 'image' to 0
				if y == 0:
					img_slice = img[y:y + full_chip_y, x:x + full_chip_x].copy()
					img_slice[0:chip_y_size, 0:chip_x_size] = 0
				else:
					img_slice = img[y:y + full_chip_y, x:x + full_chip_x].copy()
					img_slice[ovs_y_size:ovs_y_size + chip_y_size, 0:chip_x_size] = 0

				# put the clipped image into the holder image
				bias_scl[y:y + full_chip_y, x:x + full_chip_x] = (bias_scl[y:y + full_chip_y, x:x + full_chip_x] + np.median(img_slice[img_slice > 0]))

		# remove the bias from the image
		img_bias = img - bias_scl

		# remove the dark current from the image
		img_bias_dark = img_bias - dark

		# update the header
		header['BIAS_SUBT'] = 'Y'
		header['BIAS_TYPE'] = 'SCALE'
		header['DARK_SUBT'] = 'Y'

		return img_bias_dark, header