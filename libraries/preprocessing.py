from config import Configuration
from libraries.utils import Utils
import astroalign as aa
import numpy as np
import os
import pandas as pd
import twirl # deprecate
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.time import Time
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder

class Preprocessing:

	@staticmethod
	def align_img(img, header, ref_path):
		''' This function will align an image to a reference image.

		:parameter img - The image to align
		:parameter header - The header of the image
		:parameter ref_path - The path to the reference image

		:return align_img, header - The updated image and header
		'''

		ref, ref_header = fits.getdata(ref_path, header=True)

		ref_files = Utils.get_file_list(Configuration.CLEAN_DIRECTORY + Configuration.REF_DATE + '/', Configuration.FILE_EXTENSION)
		ref1, ref1_header = fits.getdata(Configuration.CLEAN_DIRECTORY + Configuration.REF_DATE + '/' + ref_files[0], header=True)

		transf, (source_list, target_list) = aa.find_transform(ref, ref1)
		align_img = aa.apply_transform(transf, img, ref1)

		return align_img[0], header

	@staticmethod
	def clip_image(img, header):
		''' This function will clip the image and remove any overscan regions.

		:parameter img - The image to clip
		:parameter header - The header of the image

		:return image_clip, header - The clipped image and the new header
		'''

		image_clip = np.zeros((Configuration.AXS_Y, Configuration.AXS_X))

		ovs_x_size = 180
		ovs_y_size = 40

		chip_x_size = 1320
		chip_y_size = 5280

		full_chip_x = chip_x_size + ovs_x_size
		full_chip_y = chip_y_size + ovs_y_size

		idx = 0
		for x in range(0, Configuration.AXS_X_RW, full_chip_x):

			idy = 0
			for y in range(0, Configuration.AXS_Y_RW, full_chip_y):

				image_clip[idy:idy + chip_y_size, idx:idx + chip_x_size] = img[y:y + chip_y_size, x:x + chip_x_size]

				idy += chip_y_size

			idx += chip_x_size

		header['OVERSCAN'] = 'removed'
		header['X_CLIP'] = ovs_x_size
		header['Y_CLIP'] = ovs_y_size

		return image_clip, header

	@staticmethod
	def correct_header(img, header):
		''' This function will plate solve the image, add the time stamp, and add the exposure time to the header

		:parameter img - This is the image you want to plate solve
		:parameter header - The header of the image

		:return img, header - The corrected image and header
		'''

		center = SkyCoord(Configuration.RA, Configuration.DEC, unit=['deg', 'deg'])
		pixel = Configuration.PIXEL_SIZE * u.arcsec
		fov = np.max(np.shape(img)) * pixel.to(u.deg)

		all_stars = twirl.gaia_radecs(center, 1.25*fov)

		all_stars = tiwlr.geometry.sparsify(all_stars, 0.01)[0:30]

		xy = twirl.find_peaks(img[500:,500:])[0:30] + 500

		if len(xy) > 30:

			wcs = twirl.compute_wcs(xy, all_stars)

			h = wcs.to_header()

			for idx, v in enumerate(h):
				header[v] = (h[idx], h.comments[idx])

			try:
				header['DATE']
			except:
				header['EXPTIME'] = Configuration.EXP_TIME
				header['DATE'] = Time.now().iso

		else:
			Utils.log('Bad image!', 'info')
			header['BADIMAGE'] = 'Y'

		return img, header

	@staticmethod
	def divide_flat(img, header):
		''' This function will divide a flat field.

		:parameter img - The image to flatten
		:parameter header - The header of the image

		:return flat_div, header - The updated image and header
		'''

		flat = fits.getdata(Configuration.CALIBRATION_DIRECTORY + 'flat.fits')

		flat_div = img / flat

		header['FLATTEN'] = 'Y'

		return flat_div, header

	@staticmethod
	def mk_bias(overwrite_bias, combine_type):
		''' This function will make the master bias frame using the provided image list and desired method.

		:parameter overwrite_bias - Y/N if you want to force an overwrite for the bias frame
		:parameter combine_type - Right now the combination type is mean, but it can be updated for whatever method you desire

		:return - The bias frame is returned and written to the calibration directory
		'''

		file_name = 'bias' + Configuration.FILE_EXTENSION

		if (os.path.isfile(Configuration.CALIBRATION_DIRECTORY + file_name) == 0) | (overwrite_bias == 'Y'):
			
			if combine_type == 'mean':
				
				chk_bias_list = Utils.get_file_list(Configuration.BIAS_DIRECTORY, Configuration.FILE_EXTENSION)

				images = pd.read_csv(Configuration.CALIBRATION_DIRECTORY + 'bias_list.csv', sep=',')
				image_list = images.apply(lambda x: Configuration.RAW_DIRECTORY + x.Date + x.File, axis=1).to_list()

				nfiles = len(image_list)

				if len(chk_bias_list) == 0:
					
					Utils.log('Generating a master bias frame from multiple files using a mean combination. There are ' + str(nfiles) + ' images to combine.', 'info')

					tmp_num = 0

					bias_header = fits.getheader(image_list[0])

					for kk in range(0, nfiles):

						bias_tmp = fits.getdata(image_list[kk]).astype('float')

						try:
							bias_image += bias_tmp
						except:
							bias_image = bias_tmp

						if (kk % 20) == 0:
							Utils.log(str(kk + 1) + ' files read in.' + str(nfiles - kk - 1) + ' files remain.', 'info')

						if ((kk % 100 == 0) & (kk > 0)) | (kk == nfiles - 1):
							Utils.log('Writing temporary file ' + str(tmp_num) + '. ' + str(nfiles - kk - 1) + ' files remain.', 'info')

							fits.writeto(Configuration.BIAS_DIRECTORY + str(tmp_num) + '_tmp_bias.fits', bias_image, bias_header, overwrite=True)

							tmp_num += 1

							del bias_image

						del bias_tmp

				else:
					Utils.log('Temporary bias files found. Skipping generating files. Delete them if you need to!', 'info')

				tmp_bias_list = Utils.get_file_list(Configuration.BIAS_DIRECTORY, Configuration.FILE_EXTENSION)

				for kk in range(0, len(tmp_bias_list)):

					bias_tmp = fits.getdata(Configuration.BIAS_DIRECTORY + tmp_bias_list[kk]).astype('float')

					try:
						bias_image_fin += bias_tmp
					except:
						bias_image_fin = bias_tmp

				del bias_tmp

				bias_image_mean = bias_image_fin / nfiles
				del bias_image_fin

				bias_header['BIAS_COMB'] = 'mean'
				bias_header['NUM_BIAS'] = nfiles

				fits.writeto(Configuration.CALIBRATION_DIRECTORY + file_name, bias_image_mean, bias_header, overwrite=True)

			else:
				Utils.log('Specific bias-combination method is unavailable, try again.', 'info')

		else:
			bias_image_mean = fits.getdata(Configuration.CALIBRATION_DIRECTORY + file_name, 0)

		return bias_image_mean

	@staticmethod
	def mk_dark(time_scale, overwrite_dark, combine_type):
		''' This function will make the master dark frame using the provided image list.

		:parameter time_scale - The length of the exposure for the dark frames
		:parameter overwrite_dark - Y/N if you want to force the current file to be overwritten
		:parameter combine_type - Either median or mean depending on how you want to combine the files

		:return - The dark frame is returned and written to the calibration directory
		'''

		file_name = str(time_scale) + 's_dark_' + Configuration.FILE_EXTENSION

		if (os.path.isfile(Configuration.CALIBRATION_DIRECTORY + file_name) == 0) | (overwrite_dark == 'Y'):

			if combine_type == 'mean':

				chk_dark_list = Utils.get_file_list(Configuration.DARK_DIRECTORY, str(time_scale) + Configuration.FILE_EXTENSION)

				images = pd.read_csv(Configuration.CALIBRATION_DIRECTORY + 'dark_list.csv', sep=',')
				image_list = images.apply(lambda x: Configuration.RAW_DIRECTORY + x.Date + x.Files, axis=1).to_list()

				nfiles = len(image_list)
				if len(chk_dark_list) == 0:

					bias = Preprocessing.mk_bias(Configuration.BIAS_DIRECTORY, combine_type='mean')

					Utils.log('Generating a master dark frame for time scale ' + str(time_scale) + 's from multiple files using a mean combination. There are ' + str(nfiles) + ' images to combine.', 'info')

					tmp_num = 0

					dark_header = fits.getheader(image_list[0])

					for kk in range(0, nfiles):

						dark_tmp = fits.getdata(image_list[kk]).astype('float')

						try:
							dark_image = (dark_image - bias) + dark_tmp
						except:
							dark_image = dark_tmp - bias

						if (kk % 20) == 0:
							Utils.loget(str(kk + 1) + ' files read in. ' + str(nfiles - kk - 1) + ' files remain.', 'info')

						if ((kk % 100 == 0) & (kk > 0)) | (kk == nfiles - 1):
							Utils.log('Writing temporary file ' + str(tmp_num) + '. ' + str(nfiles - kk - 1) + ' files remain.', 'info')

							fits.writeto(Configuration.DARK_DIRECTORY + str(tmp_num) + '_' + str(time_scale) + 's' + '_tmp_dark.fits', dark_image, dark_header, overwrite=True)
							
							del dark_image
						
						del dark_tmp

				else:
					Utils.log('Temporary dark files found. Skipping generating files.', 'info')

				tmp_dark_list = Utils.get_file_list(Configuration.DARK_DIRECTORY, Configuration.FILE_EXTENSION)

				for kk in range(0, len(tmp_dark_list)):

					dark_tmp = fits.getdata(Configuration.BIAS_DIRECTORY + tmp_dark_list[kk]).astype('float')

					try:
						dark_image_fin += dark_tmp
					except:
						dark_image_fin = dark_tmp

				del dark_tmp

				dark_image_mean = dark_image_fin / nfiles

				del dark_image_fin

				dark_header['DARK_COMB'] = 'mean'
				dark_header['NUM_DARK'] = nfiles
				dark_header['EXPTIME'] = time_scale

				fits.writeto(Configuration.CALIBRATION_DIRECTORY + file_name, dark_image_mean, dark_header, overwrite=True)

			else:
				Utils.log('Specific dark combination method is unavailable, try again.', 'info')

		else:
			dark_image_mean = fits.getdata(Configuration.CALIBRATION_DIRECTORY + file_name, 0)

		return dark_image_mean

	@staticmethod
	def mk_flat(flat_exp, dark_exp):

		return nflat_image

	@staticmethod
	def mk_nme(file, difference_image='N', clip_image='N', subtract_sky='N', subtract_bias='N', divide_flat='N', subtract_dark='N', plate_solve='N'):

		return file_name

	@staticmethod
	def subtract_bias(img, header):
		''' This function will subtract a bias frame.

		:parameter img - The image to de-bias
		:parameter header - The header of the image

		:return bias_sub, header - The updated image and header
		'''

		bias = fits.getdata(Configuration.CALIBRATION_DIRECTORY + 'bias.fits')

		bias_sub = img - bias

		header['BIAS_SUB'] = 'Y'

		return bias_sub, header

	@staticmethod
	def subtract_dark(img, header):
		''' This function will subtract a dark frame.

		:parameter img - The image to dark subtract
		:parameter header - The header of the image

		:return dark_sub, header - The updated image and header
		'''

		dark = fits.getdata(Configuration.CALIBRATION_DIRECTORY + 'dark.fits')

		dark_sub = img - dark

		header['DARK_SUB'] = 'Y'

		return dark_sub, header

	@staticmethod
	def subtract_sky(img, header, write_sky):
		''' This function will subtract the background sky from the image.

		:parameter img - The image to be cleaned
		:parameter header - The header of the image
		:parameter write_sky - Y/N if you want to write the residual background for debugging

		:return img_sub, header - The updated image and header
		'''

		sigma_clip = SigmaClip(sigma=3)
		bkg_estimator = MedianBackground()

		daofind = DAOStarFinder(fwhm=3.0, threshold=50)
		sources = daofind(img[3000:, 3000:])

		x_cen = (np.sum(sources[sources['flux'] > 0]['xcentroid'] * sources[sources['flux'] > 0]['flux']) / np.sum(sources[sources['flux'] > 0]['flux'])) + 3000
		y_cen = (np.sum(sources[sources['flux'] > 0]['ycentroid'] * sources[sources['flux'] > 0]['flux']) / np.sum(sources[sources['flux'] > 0]['flux'])) + 3000

		mask_img = np.zeros((Configuration.AXS_Y, Configuration.AXS_X))

		mask_img[int(y_cen - 1000):int(y_cen + 1000), int(x_cen - 1000):int(x_cen + 1000)] = 1
		if ((x_cen - 5200) - 200 > 0) & ((x_cen - 5200) > 0):
			mask_img[int((y_cen - 500) - 200):int((y_cen - 500) + 200),
			int((x_cen - 5200) - 200):int((x_cen - 5200) + 200)] = 1

		bkg = Background2D(img, (Configuration.PIX, Configuration.PIX), filter_size=(3,3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator, mask=mask_img)
		sky = bkg.background

		img_sub = img - sky

		header['sky_medn'] = bkg.background_median
		header['sky_sig'] = bkg.background_rms_median
		header['sky'] = np.quantile(sky, 0.25)
		header['sky_sub'] = 'yes'

		if write_sky == 'Y':
			fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'sky_background' + Configuration.FILE_EXTENSION, sky, overwrite=True)
			fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'img' + Configuration.FILE_EXTENSION, img, overwrite=True)
			fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'img_sub' + Configuration.FILE_EXTENSION, img_sub, header=header, overwrite=True)

		return img_sub, header