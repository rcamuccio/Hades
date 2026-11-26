'''
Μετρητής - The Measurer
'''

import astropy.wcs as pywcs
import ccdproc
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy.ndimage
import shutil
import subprocess
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import ascii, fits
from astropy.nddata.utils import Cutout2D
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.table import hstack, Table
from astropy.time import Time
from astropy.visualization import ImageNormalize, LinearStretch, MinMaxInterval, SqrtStretch, ZScaleInterval
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord, proj_plane_pixel_scales, skycoord_to_pixel
from photutils.aperture import aperture_photometry, ApertureStats, CircularAperture, CircularAnnulus, RectangularAperture
from photutils.background import Background2D, MedianBackground
from photutils.centroids import centroid_sources
from photutils.detection import DAOStarFinder
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint
from reproject import reproject_interp

from conf import Configuration
from lib.grapheus import Grapheus
from lib.logistes import Logistes
from lib.skopeus import Skopeus
from lib.zetetes import Zetetes

class Metretes:

	@staticmethod
	def align_frame(target_frame_path, reference_frame_path, aligned_frame_path):
		'''This function aligns a target image with a reference image.

		:p target_frame_path [str] - the target frame path
		:p reference_frame_path [str] - the reference frame path
		:p aligned_frame_path [str] - the aligned frame path

		:r - nothing
		'''

		# open the reference frame
		reference_frame = fits.open(reference_frame_path)
		reference_frame_header = reference_frame[0].header

		# open the target frame
		target_frame = fits.open(target_frame_path)
		target_frame_data = target_frame[0].data
		target_frame_header = target_frame[0].header
		target_frame_hdu = fits.PrimaryHDU(target_frame_data, header=target_frame_header)

		# align the target frame with the reference frame
		aligned_frame_data, aligned_frame_footprint = reproject_interp(target_frame_hdu, reference_frame_header)

		# correct the bad pixels
		aligned_frame_data, aligned_frame_header = Metretes.correct_pixels(aligned_frame_data, target_frame_header, report_pixels=False)

		# set the aligned frame data type
		aligned_frame_data = aligned_frame_data.astype(Configuration.DATA_TYPE)

		# modify the aligned frame header
		aligned_frame_header['ALNIM'] = ('yes', 'Status: image aligned')

		# convert the aligned frame data to an HDU object
		aligned_frame_hdu = fits.PrimaryHDU(aligned_frame_data, header=aligned_frame_header)

		# save the aligned frame HDU object
		if not os.path.exists(aligned_frame_path):
			aligned_frame_hdu.writeto(aligned_frame_path)

	@staticmethod
	def correct_pixels(frame_data, frame_header, report_pixels=True):
		'''This function corrects the bad pixels of a data frame, and reports the individual bad pixels, if desired.

		:p frame_data [numpy.ndarray] - the image data array
		:p frame_header [numpy.ndarray] - the image header 
		:p report_pixels [bool] - a toggle for reporting bad pixels

		:r correct_data [numpy.ndarray] - the corrected image data array
		:r correct_header [astropy.io.fits.header] - the corrected image header
		'''

		# set the array dimensions
		try:
			overscan_present = frame_header['OVSRM']
			x_axis = Skopeus.get_axis_x_raw()
			y_axis = Skopeus.get_axis_y_raw()
		except KeyError:
			x_axis = Skopeus.get_axis_x()
			y_axis = Skopeus.get_axis_y()

		# report the bad pixels
		if report_pixels:
			nan_idx = 0
			inf_idx = 0
			for y in range(0, y_axis):
				for x in range(0, x_axis):
					if np.isnan(frame_data[y][x]):
						nan_idx += 1
					elif np.isinf(frame_data[y][x]):
						inf_idx += 1

		# calculate the background statistics
		background_mn, background_md, background_sd = sigma_clipped_stats(frame_data, sigma=Configuration.SIGMA_BACKGROUND)

		# create a holder array for the bad pixels
		bad_pixels = np.isnan(frame_data) + np.isinf(frame_data)

		# replace the bad pixels with the background mean
		frame_data[bad_pixels] = background_mn
		correct_data = frame_data

		# update the header
		correct_header = frame_header
		correct_header['BPCOR'] = ('yes', 'Status: pixels corrected')

		return correct_data, correct_header

	@staticmethod
	def divide_flat(raw_frame_data, raw_frame_header, flat_frame_path):
		'''This function corrects an image data array with a normalized flatfield data array, saves the file, and returns both the corrected data array and header.

		:p raw_frame_data [numpy.ndarray] - the raw image data array
		:p raw_frame_header [astropy.io.fits.header] - the raw image header
		:p flat_frame_path [str] - the flatfield image path

		:r clean_frame_data [numpy.ndarray] - the cleaned image data array
		:r clean_frame_header [astropy.io.fits.header] - the updated image header
		'''

		# open the flatfield frame
		flat_frame = fits.open(flat_frame_path)
		flat_frame_data = flat_frame[0].data
		flat_frame_header = flat_frame[0].header

		# divide the flat
		clean_frame_data = raw_frame_data / flat_frame_data

		# update the header
		clean_frame_header = raw_frame_header
		clean_frame_header['FLTDV'] = ('yes', 'Status: flatfield corrected')

		return clean_frame_data, clean_frame_header

	@staticmethod
	def extract_sources(frame_data, frame_header, table_path):
		'''This function extracts sources from a FITS image data array and returns a source table.

		:p frame_data [numpy.ndarray] - the image data array
		:p frame_header [astropy.io.fits.header] - the image header
		:p source_table_path [str] - the source table path

		:r source_table [astropy.io] - the source table
		'''

		line_preamble = '\t\033[1m[Μετρητής | extract_sources]\033[0m :'

		# extract the sources
		if not os.path.exists(source_table_path):
			print(line_preamble, 'Extracting sources')

			# read the WCS information
			wcs = WCS(frame_header)

			# create a mask
			mask, boxes = Metretes.make_mask(frame_data)

			# calculate the background statistics
			frame_mn, frame_md, frame_sd = sigma_clipped_stats(frame_data, mask=mask, sigma=Configuration.SIGMA_BACKGROUND)

			# extract the sources into a table
			daofind = DAOStarFinder(fwhm=Configuration.FWHM, threshold=Configuration.SIGMA_SOURCE)
			source_table = daofind(frame_data, mask=mask)

			ra_list = []
			dec_list = []

			for ln in table:
				# get the pixel coordinate
				x = ln['xcentroid']
				y = ln['ycentroid']

				# transform the pixel coordinate to sky coordinate
				sky_pos = pixel_to_skycoord(x, y, wcs=wcs)
				ra = sky_pos.ra.deg
				dec = sky_pos.dec.deg

				# add the sky coordinate to list
				ra_list.append(ra)
				dec_list.append(dec)

			# add the sky coordinates to the table
			source_table['ra'] = ra_list
			source_table['dec'] = dec_list

			# save the source table
			source_table.write(source_table_path, format=Configuration.TABLE_FORMAT)

		else:
			print(line_preamble, 'Reading existing source table')
			source_table = Table.read(source_table_path, format=Configuration.TABLE_FORMAT)

		return source_table

	@staticmethod
	def make_dark(date, dark_type):
		'''This function creates a combined master dark frame from a series of individual dark frames taken on a given night and with a given exposure time.

		:p date [str] - the date string [e.g. yyyy-mm-dd]
		:p dark_type [str] - the dark frame type based on exposure time [flat, light]

		:r master_data [numpy.ndarray] - the master dark frame
		:r master_header [astropy.io.fits.header] - the master dark header
		'''

		line_preamble = '\t\033[1m[Μετρητής | make_dark]\033[0m :'

		# set the master directory
		master_directory = Configuration.OUTPUT_DIRECTORY + 'dark/' + date + '/'

		# create the master dark frame directory
		if not os.path.exists(master_directory):
			os.mkdir(master_directory)

		# set the master dark frame and table paths
		if dark_type == 'flat':
			master_path = master_directory + 'mdf_' + date + Configuration.DATA_EXTENSION
			dark_table_path = master_directory + 'stat_mdf_' + date + Configuration.TABLE_EXTENSION
		elif dark_type == 'light':
			master_path = master_directory + 'mdl_' + date + Configuration.DATA_EXTENSION
			dark_table_path = master_directory + 'stat_mdl_' + date + Configuration.TABLE_EXTENSION

		# create the master dark frame
		if not os.path.exists(master_path):
			print(line_preamble, 'Generating a master dark (' + dark_type + ') for', date)

			dark_ccddata_list = []
			keep_idx = 0

			# load the dark frames to combine
			flat_directory = Configuration.INPUT_DIRECTORY + 'dark/' + date + '/'
			os.chdir(flat_directory)

			# set the search string
			if dark_type == 'flat':
				search_str = 'dark_' + str(int(Skopeus.get_exposure_time('dark-f'))) + 's*' + Configuration.DATA_EXTENSION
			elif dark_type == 'light':
				search_str = 'dark_' + str(int(Skopeus.get_exposure_time('dark-l'))) + 's*' + Configuration.DATA_EXTENSION

			# get the number of dark frames
			dark_frame_list = sorted(glob.glob(search_str))
			dark_number = len(dark_frame_list)

			print(line_preamble, 'Found', str(dark_number), 'frames')

			# create the dark frame table
			tbl_col = ('frm', 'dte', 'jd', 'tex', 'mn', 'md', 'sd', 'kp')
			tbl_dt = (str, str, float, float, float, float, float, str)
			dark_table = Table(names=tbl_col, dtype=tbl_dt)

			# process each dark frame
			for frm in range(dark_number):
				# open the dark frame
				dark_frame = fits.open(dark_frame_list[frm])
				dark_frame_data = dark_frame[0].data
				dark_frame_header = dark_frame[0].header

				# get the time information
				time = Time(dark_frame_header['DATE'], format='isot', scale='utc')
				jd = time.jd

				# calculate the background statistics
				dark_mn, dark_md, dark_sd = sigma_clipped_stats(dark_frame_data, sigma=Configuration.SIGMA_BACKGROUND)

				# remove the overscan
				dark_frame_data, dark_frame_header = Metretes.remove_overscan(dark_frame_data, dark_frame_header)

				# add the dark frame object to list
				dark_ccddata = ccdproc.CCDData(dark_frame_data, unit=Configuration.UNIT)
				dark_ccddata_list.append(dark_ccddata)

				# count the dark
				keep_idx += 1
				keep_flag = 'Y'
				print('\t\t+ Frame:', dark_frame_list[frm])

				# add a row to the dark table
				tbl_row = (dark_frame_list[frm], date, jd, 300., dark_mn, dark_md, dark_sd, keep_flag)
				dark_table.add_row(tbl_row)

			# combine the dark frames
			if dark_number > 0:
				# save the dark table
				if os.path.exists(dark_table_path):
					dark_table = ascii.read(dark_table_path, format=Configuration.TABLE_FORMAT)
				else:
					ascii.write(dark_table, dark_table_path, format=Configuration.TABLE_FORMAT)

				print(line_preamble, 'Combining', keep_idx, 'frames into a master (' + Configuration.COMBINE_METHOD + ') dark')

				# combine the dark frames
				master_ccddata = ccdproc.combine(dark_ccddata_list, method=Configuration.COMBINE_METHOD, unit=Configuration.UNIT, mem_limit=Configuration.MEMORY_LIMIT, dtype=Configuration.DATA_TYPE)

				# convert the ccddata object to an array
				master_data = np.asarray(master_ccddata)

				# calculate the background statistics
				master_mn, master_md, master_sd = sigma_clipped_stats(master_data, sigma=Configuration.SIGMA_BACKGROUND)

				# modify the master header
				master_hdu = fits.PrimaryHDU(master_data)
				master_header = master_hdu.header
				master_header['SIMPLE'] = ('T', 'Conform to FITS standard')
				master_header['BITPIX'] = (32, 'Number of bits per data pixel')
				master_header['NAXIS'] = (2, 'Number of axes')
				master_header['NAXIS1'] = (10560, 'Image width')
				master_header['NAXIS2'] = (10560, 'Image height')
				master_header['BZERO'] = (32768, 'Offset for unsigned short')
				master_header['BSCALE'] = (1, 'Default scaling factor')
				master_header['DATE'] = (date, 'Date of combined darks')
				master_header['EXPTIME'] = (300, 'Exposure time')
				master_header['IMAGETYP'] = ('DARK', 'Type of image')
				master_header['DRKCB'] = ('yes', 'Status: darks combined')
				master_header['DRKNM'] = (keep_idx, 'Number of combined dark frames')
				master_header['MDMD'] = (master_md, 'Master dark median [ADU]')
				master_header['MDMN'] = (master_mn, 'Master dark mean [ADU]')
				master_header['MDSD'] = (master_sd, 'Master dark standard deviation [ADU]')
				master_header['OVSRM'] = ('yes', 'Status: overscan removed')

				# save the master frame
				fits.writeto(master_path, master_data, master_header)

				print(line_preamble, 'Master dark generated for', date)

			else:
				print(line_preamble, 'No dark frames available')

		else:
			print(line_preamble, 'Master dark frame already exists for', date)

	@staticmethod
	def make_flat(date):
		'''This function creates a normalized flatfield frame from a series of individual flat frames taken on a given night and with a given exposure time. The function uses a master dark frame which has the same exposure time as the flat frames.

		:p date [str] - the date string (e.g. yyyy-mm-dd)

		:r - nothing
		'''

		line_preamble = '\t\033[1m[Μετρητής | make_flat]\033[0m :'

		# set the flat path
		flatfield_directory = Configuration.OUTPUT_DIRECTORY + 'flat/' + date + '/'
		flatfield_path = flatfield_directory + 'ff_' + date + Configuration.DATA_EXTENSION

		# create the output directory
		if not os.path.exists(flatfield_directory):
			os.mkdir(flatfield_directory)

		# set the dark path
		dark_directory = Configuration.OUTPUT_DIRECTORY + 'dark/' + date + '/'
		dark_path = dark_directory + 'mdf_' + date + Configuration.DATA_EXTENSION

		# set the table path
		flat_table_path = flatfield_directory + 'stat_ff_' + date + Configuration.TABLE_EXTENSION

		# create the normalized flatfield
		if not os.path.exists(flatfield_path):
			print(line_preamble, 'Generating a normalized flatfield for', date)

			flat_ccddata_list = []
			keep_idx = 0

			# load the master dark
			dark_ccddata = ccdproc.fits_ccddata_reader(dark_path, unit=Configuration.UNIT)
			dark_exptime = Skopeus.get_exposure_time('dark-f') * u.second
			flat_exptime = Skopeus.get_exposure_time('flat') * u.second

			# load the flat frames to combine
			flat_frame_directory = Configuration.INPUT_DIRECTORY + 'flat/' + date + '/'
			os.chdir(flat_frame_directory)
			search_str = 'flat_' + str(int(Skopeus.get_exposure_time('flat'))) + 's*' + Configuration.DATA_EXTENSION
			flat_frame_list = sorted(glob.glob(search_str))

			# get the number of flat frames
			flat_number = len(flat_frame_list)

			print(line_preamble, 'Found', str(flat_number), 'frames')

			# create the flat frame table
			tbl_col = ('frm', 'dte', 'jd', 'tex', 'mn', 'md', 'sd', 'kp')
			tbl_dt = (str, str, float, float, float, float, float, str)
			flat_table = Table(names=tbl_col, dtype=tbl_dt)

			# process each flat frame
			for frm in range(flat_number):
				# open the flat frame
				flat_frame = fits.open(flat_frame_list[frm])
				flat_frame_data = flat_frame[0].data
				flat_frame_header = flat_frame[0].header

				# get the time information
				time = Time(flat_frame_header['DATE'], format='isot', scale='utc')
				jd = time.jd

				# calculate the background statistics
				flat_mn, flat_md, flat_sd = sigma_clipped_stats(flat_frame_data, sigma=Configuration.SIGMA_BACKGROUND)

				# filter out flats with backgrounds outside the limits
				if (flat_mn > Configuration.FLAT_MINIMUM) and (flat_mn < Configuration.FLAT_MAXIMUM):
					# remove the overscan
					flat_frame_data, flat_frame_header = Metretes.remove_overscan(flat_frame_data, flat_frame_header)

					# subtract the master dark
					flat_ccddata = ccdproc.CCDData(flat_frame_data, unit=Configuration.UNIT)
					flat_ccddata = ccdproc.subtract_dark(flat_ccddata, dark_ccddata, dark_exposure=dark_exptime, data_exposure=flat_exptime)

					# add the flat frame to list
					flat_ccddata_list.append(flat_ccddata)

					# count the flat
					keep_idx += 1
					keep_flag = 'Y'
					print('\t\t+ Frame:', flat_frame_list[frm])

				else:
					keep_flag = 'N'
					print('\t\t- Frame:', flat_frame_list[frm])

				# add the row to the flat table
				tbl_row = (flat_frame_list[frm], date, jd, 5., flat_mn, flat_md, flat_sd, keep_flag)
				flat_table.add_row(tbl_row)

			# combine the flat frames
			if flat_number > 0:
				# save the flat table
				if os.path.exists(flat_table_path):
					flat_table = ascii.read(flat_table_path, format=Configuration.TABLE_FORMAT, delimiter='|')
				else:
					ascii.write(flat_table, flat_table_path, format=Configuration.TABLE_FORMAT)

				print(line_preamble, 'Combining', keep_idx, 'frames into a normalized flatfield')

				# combine the flat frames
				combined_ccddata = ccdproc.combine(flat_ccddata_list, method=Configuration.COMBINE_METHOD, unit=Configuration.UNIT, mem_limit=Configuration.MEMORY_LIMIT, dtype=Configuration.DATA_TYPE)

				# convert the ccddata object to an array
				combined_data = np.asarray(combined_ccddata)

				# calculate the combined flat background statistics
				combined_mn, combined_md, combined_sd = sigma_clipped_stats(combined_data, sigma=Configuration.SIGMA_BACKGROUND)

				# normalize the combined flat
				print(line_preamble, 'Normalizing the combined flat')
				flatfield_data = combined_data / combined_mn
				flatfield_data = flatfield_data.astype(Configuration.DATA_TYPE)

				# calculate the normalized flatfield background statistics
				flatfield_mn, flatfield_md, flatfield_sd = sigma_clipped_stats(flatfield_data, sigma=Configuration.SIGMA_BACKGROUND)

				# modify the flatfield header
				flatfield_hdu = fits.PrimaryHDU(flatfield_data)
				flatfield_header = flatfield_hdu.header
				flatfield_header['SIMPLE'] = ('T', 'Conform to FITS standard')
				flatfield_header['BITPIX'] = (32, 'Number of bits per data pixel')
				flatfield_header['NAXIS'] = (2, 'Number of axes')
				flatfield_header['NAXIS1'] = (10560, 'Image width')
				flatfield_header['NAXSI2'] = (10560, 'Image height')
				flatfield_header['BZERO'] = (32768, 'Offset for unsigned short')
				flatfield_header['BSCALE'] = (1, 'Default scaling factor')
				flatfield_header['DATE'] = (date, 'Date of combined flats')
				flatfield_header['EXPTIME'] = (5, 'Exposure time')
				flatfield_header['IMAGETYP'] = ('FLAT', 'Type of image')
				flatfield_header['CFMD'] = (combined_md, 'Combined flat median [ADU]')
				flatfield_header['CFMN'] = (combined_mn, 'Combined flat mean [ADU]')
				flatfield_header['CFSD'] = (combined_sd, 'Combined flat standard deviation [ADU]')
				flatfield_header['FFMD'] = (flatfield_md, 'Flatfield median [ADU]')
				flatfield_header['FFMN'] = (flatfield_mn, 'Flatfield mean [ADU]')
				flatfield_header['FFSD'] = (flatfield_sd, 'Flatfield standard deviation [ADU]')
				flatfield_header['FTCB'] = ('yes', 'Status: flats combined')
				flatfield_header['FTNM'] = (keep_idx, 'Number of combined flat frames')
				flatfield_header['OVSRM'] = ('yes', 'Status: overscan removed')

				# save the flatfield frame
				fits.writeto(flatfield_path, flatfield_data, flatfield_header)

				print(line_preamble, 'Normalized flatfield generated for', date)

			else:
				print(line_preamble, 'No flat frames available')

		else:
			print(line_preamble, 'Normalized flatfield already exists for', date)

	@staticmethod
	def make_mask(frame_data):
		'''This function creates a mask array for a FITS image data array.

		:p frame_data [numpy.ndarray] - the image data array

		:r mask_data [numpy.ndarray] - the mask data array
		:r boxes [photutils.aperture.rectangle.RectangularAperture] - the set of mask boxes
		'''

		# generate mask holder
		frame_shape = frame_data.shape
		dy = frame_shape[0]
		dx = frame_shape[1]
		mask_data = np.zeros(frame_shape, dtype=bool)

		# set the edge region
		edge = 101 // Skopeus.get_binning()

		# create the left edge mask
		mask_data[0:dy, 0:edge] = True
		eml_cen_x = edge // 2
		eml_cen_y = dy // 2
		eml_box = RectangularAperture((eml_cen_x, eml_cen_y), edge, dy, theta=0.)

		# create the top edge mask
		mask_data[dy-edge:dy, 0:dx] = True		
		emt_cen_x = dx // 2
		emt_cen_y = edge // 2
		emt_box = RectangularAperture((emt_cen_x, emt_cen_y + dy - edge), dx, edge, theta=0.)

		# create the right edge mask
		mask_data[0:dy, dx-edge:dx] = True
		emr_cen_x = edge // 2
		emr_cen_y = dy // 2
		emr_box = RectangularAperture((emr_cen_x + dx - edge, emr_cen_y), edge, dy, theta=0.)

		# create the bottom edge mask
		mask_data[0:edge, 0:dx] = True
		emb_cen_x = dx // 2
		emb_cen_y = edge // 2
		emb_box = RectangularAperture((emb_cen_x, emb_cen_y), dx, edge, theta=0.)

		# create the overscan boundary masks
		w = 20
		tx = Skopeus.get_tap_x()
		ty = Skopeus.get_tap_y()

		# set the vertical bands
		mask_data[0:dy, tx-w:tx+w] = True
		mask_data[0:dy, (2*tx)-w:(2*tx)+w] = True
		mask_data[0:dy, (3*tx)-w:(3*tx)+w] = True
		mask_data[0:dy, (4*tx)-w:(4*tx)+w] = True
		mask_data[0:dy, (5*tx)-w:(5*tx)+w] = True
		mask_data[0:dy, (6*tx)-w:(6*tx)+w] = True
		mask_data[0:dy, (7*tx)-w:(7*tx)+w] = True

		# set the horizontal band
		mask_data[ty-w:ty+w, 0:dx] = True

		# create the box dictionary
		boxes = {}
		boxes['eml_box'] = eml_box
		boxes['emt_box'] = emt_box
		boxes['emr_box'] = emr_box
		boxes['emb_box'] = emb_box

		return mask_data, boxes

	@staticmethod
	def make_stack(date, field, frame_table):
		'''This function creates a master stack frame from a time series of individual frames for a particular date and field.

		:p date [str] - the date string [e.g. yyyy-mm-dd]
		:p field [str] - the field ID string [e.g. xx.xxx]
		:p frame_table [astropy.table] - the frame metadata table

		:r stack_data [numpy.ndarray] - the stacked image data array
		:r stack_header [astropy.io.fits.header] - the stacked image header
		'''

		line_preamble = '\t\033[1m[Μετρητής | make_stack]\033[0m :'

		# set the stack directory
		stack_directory = Configuration.OUTPUT_DIRECTORY + 'clean/' + date + '/FIELD_' + field + '/'
		
		# set the stack names
		stack_name = 'STACK_FIELD_' + field + '_' + date
		stack_filename = stack_name + Configuration.DATA_EXTENSION

		# set the stack path
		stack_path = stack_directory + stack_filename

		# set the image paths
		stack_image_path = stack_directory + 'img/' + stack_name + Configuration.IMAGE_EXTENSION
		stack_background_path = stack_directory + 'bkg/' + stack_name + Configuration.IMAGE_EXTENSION

		stack_image_histogram_path = stack_directory + 'hst/img_' + stack_name + Configuration.IMAGE_EXTENSION
		stack_background_histogram_path = stack_directory + 'hst/bkg_' + stack_name + Configuration.IMAGE_EXTENSION

		# make the stack
		if not os.path.exists(stack_path):
			# get the individual frames
			os.chdir(stack_directory)
			align_list = sorted(glob.glob('cln_*' + Configuration.DATA_EXTENSION))
			align_number = len(align_list)

			# align the frames
			print(line_preamble, 'Aligning frames')

			reference_frame = align_list[0]
			reference_frame_path = stack_directory + reference_frame

			align_reference_frame = 'aln_' + reference_frame
			align_reference_frame_path = stack_directory + align_reference_frame

			reference_frame = fits.open(reference_frame_path)
			reference_frame_data = reference_frame[0].data
			reference_frame_header = reference_frame[0].header

			align_reference_frame_hdu = fits.PrimaryHDU(reference_frame_data, header=reference_frame_header)

			if not os.path.exists(align_reference_frame_path):
				align_reference_frame_hdu.writeto(align_reference_frame_path)

			for frm in range(align_number):
				if frm > 0:
					print('\t\tAligning', align_list[frm])
					target_frame_path = stack_directory + align_list[frm]
					align_frame_path = stack_directory + 'aln_' + align_list[frm]

					if not os.path.exists(align_frame_path):
						Metretes.align_frame(target_frame_path, reference_frame_path, align_frame_path)

			# combine the frames
			stack_list = sorted(glob.glob('aln_*' + Configuration.DATA_EXTENSION))
			stack_number = len(stack_list)

			print(line_preamble, 'Generating', stack_number, 'frames into a master stack (FIELD_' + field + ', ' + date + ')')

			stack_ccddata = ccdproc.combine(stack_list, method=Configuration.COMBINE_METHOD, unit=Configuration.UNIT, mem_limit=Configuration.MEMORY_LIMIT, dtype=Configuration.DATA_TYPE)

			stack_data = np.asarray(stack_ccddata)
			stack_data = stack_data.astype(Configuration.DATA_TYPE)

			# draw the stack image
			if not os.path.exists(stack_image_path):
				Grapheus.draw_field_diagram(stack_data, None, None, stack_image_path)

			# draw a stack histogram
			if not os.path.exists(stack_image_histogram_path):
				Grapheus.draw_histogram_field(stack_data, stack_image_histogram_path)

			# calculate the background statistics
			bkg2d_data, bkg2d_md, bkg2d_sd, bkgsc_mn, bkgsc_md, bkgsc_sd = Metretes.measure_background(stack_data)

			# draw a background image
			if not os.path.exists(stack_background_path):
				Grapheus.draw_field_diagram(bkg2d_data, None, None, stack_background_path)

			# draw a background histogram
			if not os.path.exists(stack_background_histogram_path):
				Grapheus.draw_histogram_field(bkg2d_data, stack_background_histogram_path)

			# update the stack header
			stack_hdu = fits.PrimaryHDU(stack_data)
			stack_header = stack_hdu.header
			stack_header['SIMPLE'] = ('T', 'Conform to FITS standard')
			stack_header['BITPIX'] = (32, 'Number of bits per data pixel')
			stack_header['NAXIS'] = (2, 'Number of axes')
			stack_header['NAXIS1'] = (10560, 'Image width')
			stack_header['NAXIS2'] = (10560, 'Image height')
			stack_header['BZERO'] = (32768, 'Offset for unsigned short')
			stack_header['BSCALE'] = (1, 'Default scaling factor')
			stack_header['DATE'] = (date, 'Date of combined lights')
			stack_header['EXPTIME'] = (300., 'Exposure time')
			stack_header['IMAGETYP'] = ('LIGHT', 'Type of image')
			stack_header['STKCB'] = ('yes', 'Status: lights combined')
			stack_header['STKJD'] = (np.median(frame_table['jd'].tolist()), 'Median Julian date of combined light frames')
			stack_header['STKNM'] = (stack_number, 'Number of combined light frames')

			# save the stack frame
			if not os.path.exists(stack_path):
				fits.writeto(stack_path, stack_data, stack_header)

			# solve the stack frame
			print(line_preamble, 'Solving field')
			Metretes.solve_field(stack_filename, stack_directory, field, verbose=Configuration.VERBOSE)

			# open the solved stack frame
			stack_frame = fits.open(stack_path)
			stack_header = stack_frame[0].header

			print(line_preamble, 'Stack generated (FIELD_' + field + ', ' + date + ')')

		# read the existing stack
		else:
			print(line_preamble, 'Reading existing stack (FIELD_' + field + ', ' + date + ')')
			stack_frame = fits.open(stack_path)
			stack_data = stack_frame[0].data
			stack_header = stack_frame[0].header

		# clean up the aligned frames
		if not Configuration.KEEP_ALIGN_FRAMES:
			# change to the stack directory
			os.chdir(stack_directory)

			# get any aligned frames
			frame_list = glob.glob('aln_*' + Configuration.DATA_EXTENSION)

			# remove any aligned frames
			if len(frame_list) > 0:
				for frm in frame_list:
					subprocess.run(['rm', frm])

		return stack_data, stack_header

	@staticmethod
	def measure_background(frame_data):
		'''This function measures the background signal from a FITS image data array.

		:p frame_data [numpy.ndarray] - the image data array

		:r bkg2d_data [numpy.ndarray] - the background data array
		:r bkg2d_md [float] - the background median (mesh) [ADU]
		:r bkg2d_sd [float] - the background error (mesh) [ADU]
		:r bkgsc_mn [float] - the background mean (clip) [ADU]
		:r bkgsc_md [float] - the background median (clip) [ADU]
		:r bkgsc_sd [float] - the background error (clip) [ADU]
		'''

		# clip the sources
		sigma_clip = SigmaClip(sigma=Configuration.SIGMA_SOURCE)
		threshold = detect_threshold(frame_data, nsigma=Configuration.SIGMA_BACKGROUND, sigma_clip=sigma_clip)
		segment_img = detect_sources(frame_data, threshold, npixels=Configuration.MINIMUM_PIXELS)
		footprint = circular_footprint(radius=Configuration.FOOTPRINT_RADIUS)
		mask = segment_img.make_source_mask(footprint=footprint)

		# estimate the background
		bkg_estimator = MedianBackground()
		bkg = Background2D(frame_data, box_size=Configuration.BOX_SIZE, mask=mask, filter_size=Configuration.FILTER_SIZE, sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
		bkg2d_data = bkg.background
		bkg2d_md = bkg.background_median
		bkg2d_sd = bkg.background_rms_median

		# calculate the sigma-clipped statistics
		bkgsc_mn, bkgsc_md, bkgsc_sd = sigma_clipped_stats(frame_data, sigma=Configuration.SIGMA_BACKGROUND)
		bkg2d_md = float(bkg2d_md)
		bkg2d_sd = float(bkg2d_sd)
		bkgsc_mn = float(bkgsc_mn)
		bkgsc_md = float(bkgsc_md)
		bkgsc_sd = float(bkgsc_sd)

		return bkg2d_data, bkg2d_md, bkg2d_sd, bkgsc_mn, bkgsc_md, bkgsc_sd

	@staticmethod
	def remove_overscan(frame_data, frame_header):
		'''This function removes the overscan region from a FITS image and updates the FITS header.

		:p frame_data [numpy.ndarray] - the image data array
		:p frame_header [astropy.io.fits.header] - the image header

		:r clip_data [numpy.ndarray] - the clipped image data array
		:r clip_header [astropy.io.fits.header] - the updated image header
		'''

		# create the clipped data array
		clip_data = np.zeros((Skopeus.get_axis_y(), Skopeus.get_axis_x()))

		# get the overscan size
		ovs_x = Skopeus.get_overscan_x()
		ovs_y = Skopeus.get_overscan_y()

		# get the tap size
		tap_x = Skopeus.get_tap_x()
		tap_y = Skopeus.get_tap_y()

		# calculate the full tape size
		full_tap_x = tap_x + ovs_x
		full_tap_y = tap_y + ovs_y

		# clip the data array
		idx = 0
		for x in range(0, Skopeus.get_axis_x_raw(), full_tap_x):
			idy = 0
			for y in range(0, Skopeus.get_axis_y_raw(), full_tap_y):
				if y == 0:
					clip_data[idy:idy+tap_y, idx:idx+tap_x] = frame_data[y:y+tap_y, x:x+tap_x]
				else:
					clip_data[idy:idy+tap_y, idx:idx+tap_x] = frame_data[y+ovs_y:y+ovs_y+tap_y, x:x+tap_x]
				idy += tap_y
			idx += tap_x

		# change the data type
		clip_data = np.float32(clip_data)

		# update the header
		clip_header = frame_header
		clip_header['OVSRM'] = ('yes', 'Status: overscan removed')

		return clip_data, clip_header

	@staticmethod
	def setup_source_extractor(save_directory):
		'''This function generates the core files (.conv, .nnw, .param, .sex) necessary to run source extractor at the given save directory.

		:p save_directory [str] - the save directory where the source extractor files are generated

		:r - nothing
		'''

		# make the convolution mask file
		default_conv_path = save_directory + 'default.conv'
		if not os.path.isfile(default_conv_path):
			ln01 = 'CONV NORM\n'
			ln02 = '# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.\n'
			ln03 = '1 2 1\n'
			ln04 = '2 4 2\n'
			ln05 = '1 2 1'
			lns = [ln01, ln02, ln03, ln04, ln05]
			default_conv = open(default_conv_path, 'w')
			for ln in lns:
				default_conv.write(ln)
			default_conv.close()

		# make the neural network weights file
		default_nnw_path = save_directory + 'default.nnw'
		if not os.path.isfile(default_nnw_path):
			ln01 = 'NNW\n'
			ln02 = '# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)\n'
			ln03 = '# inputs:  9 for profile parameters + 1 for seeing.\n'
			ln04 = '# outputs: ``Stellarity index'' (0.0 to 1.0)\n'
			ln05 = '# Seeing FWHM range: from 0.025 to 5.5'' (images must have 1.5 < FWHM < 5 pixels)\n'
			ln06 = '# Optimized for Moffat profiles with 2<= beta <= 4.\n\n'
			ln07 = ' 3 10 10 1\n\n'
			ln08 = '-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01 -5.22453e-01  4.61342e-02  8.31957e-01  2.15505e+00  2.64769e-01\n'
			ln09 = ' 3.03477e+00  2.69561e+00  3.16188e+00  3.34497e+00  3.51885e+00  3.65570e+00  3.74856e+00  3.84541e+00  4.22811e+00  3.27734e+00\n\n'
			ln10 = '-3.22480e-01 -2.12804e+00  6.50750e-01 -1.11242e+00 -1.40683e+00 -1.55944e+00 -1.84558e+00 -1.18946e-01  5.52395e-01 -4.36564e-01 -5.30052e+00\n'
			ln11 = ' 4.62594e-01 -3.29127e+00  1.10950e+00 -6.01857e-01  1.29492e-01  1.42290e+00  2.90741e+00  2.44058e+00 -9.19118e-01  8.42851e-01 -4.69824e+00\n'
			ln12 = '-2.57424e+00  8.96469e-01  8.34775e-01  2.18845e+00  2.46526e+00  8.60878e-02 -6.88080e-01 -1.33623e-02  9.30403e-02  1.64942e+00 -1.01231e+00\n'
			ln13 = ' 4.81041e+00  1.53747e+00 -1.12216e+00 -3.16008e+00 -1.67404e+00 -1.75767e+00 -1.29310e+00  5.59549e-01  8.08468e-01 -1.01592e-02 -7.54052e+00\n'
			ln14 = ' 1.01933e+01 -2.09484e+01 -1.07426e+00  9.87912e-01  6.05210e-01 -6.04535e-02 -5.87826e-01 -7.94117e-01 -4.89190e-01 -8.12710e-02 -2.07067e+01\n'
			ln15 = '-5.31793e+00  7.94240e+00 -4.64165e+00 -4.37436e+00 -1.55417e+00  7.54368e-01  1.09608e+00  1.45967e+00  1.62946e+00 -1.01301e+00  1.13514e-01\n'
			ln16 = ' 2.20336e-01  1.70056e+00 -5.20105e-01 -4.28330e-01  1.57258e-03 -3.36502e-01 -8.18568e-02 -7.16163e+00  8.23195e+00 -1.71561e-02 -1.13749e+01\n'
			ln17 = ' 3.75075e+00  7.25399e+00 -1.75325e+00 -2.68814e+00 -3.71128e+00 -4.62933e+00 -2.13747e+00 -1.89186e-01  1.29122e+00 -7.49380e-01  6.71712e-01\n'
			ln18 = '-8.41923e-01  4.64997e+00  5.65808e-01 -3.08277e-01 -1.01687e+00  1.73127e-01 -8.92130e-01  1.89044e+00 -2.75543e-01 -7.72828e-01  5.36745e-01\n'
			ln19 = '-3.65598e+00  7.56997e+00 -3.76373e+00 -1.74542e+00 -1.37540e-01 -5.55400e-01 -1.59195e-01  1.27910e-01  1.91906e+00  1.42119e+00 -4.35502e+00\n\n'
			ln20 = '-1.70059e+00 -3.65695e+00  1.22367e+00 -5.74367e-01 -3.29571e+00  2.46316e+00  5.22353e+00  2.42038e+00  1.22919e+00 -9.22250e-01 -2.32028e+00\n\n\n'
			ln21 = ' 0.00000e+00\n'
			ln22 = ' 1.00000e+00'
			lns = [ln01, ln02, ln03, ln04, ln05, ln06, ln07, ln08, ln09, ln10, ln11, ln12, ln13, ln14, ln15, ln16, ln17, ln18, ln19, ln20, ln21, ln22]
			default_nnw = open(default_nnw_path, 'w')
			for ln in lns:
				default_nnw.write(ln)
			default_nnw.close()

		# make the parameter file
		default_param_path = save_directory + 'default.param'
		if not os.path.isfile(default_param_path):
			number = 'NUMBER\n'
			alphapeak_j2000 = 'ALPHAPEAK_J2000\n'
			deltapeak_j2000 = 'DELTAPEAK_J2000\n'
			xpeak_image = 'XPEAK_IMAGE\n'
			ypeak_image = 'YPEAK_IMAGE\n'
			flux_growth = 'FLUX_GROWTH\n'
			fluxerr_best = 'FLUXERR_BEST\n'
			flux_growthstep = 'FLUX_GROWTHSTEP\n'
			fwhm_image = 'FWHM_IMAGE\n'
			fwhm_world = 'FWHM_WORLD\n'
			lns = [number, alphapeak_j2000, deltapeak_j2000, xpeak_image, ypeak_image, flux_growth, fluxerr_best, flux_growthstep, fwhm_image, fwhm_world]
			default_param = open(default_param_path, 'w')
			for ln in lns:
				default_param.write(ln)
			default_param.close()

		# make the configuration file (default configuration file for SExtractor 2.12.4 EB 2010-10-10)
		default_sex_path = save_directory + 'default.sex'
		if not os.path.isfile(default_sex_path):
			catalog_name = 'CATALOG_NAME test.cat\n'
			catalog_type = 'CATALOG_TYPE ASCII_HEAD\n'
			parameters_name = 'PARAMETERS_NAME default.param\n'
			detect_type = 'DETECT_TYPE CCD\n'
			detect_minarea = 'DETECT_MINAREA 3\n'
			detect_thresh = 'DETECT_THRESH 5.0\n'
			analysis_thresh = 'ANALYSIS_THRESH 5.0\n'
			detection_filter = 'FILTER Y\n'
			detection_filter_name = 'FILTER_NAME default.conv\n'
			deblend_nthresh = 'DEBLEND_NTHRESH 32\n'
			deblend_mincont = 'DEBLEND_MINCONT 0.005\n'
			clean = 'CLEAN Y\n'
			clean_param = 'CLEAN_PARAM 1.0\n'
			weight_type = 'WEIGHT_TYPE NONE\n'
			weight_image = 'WEIGHT_IMAGE weight.fits\n'
			flag_image = 'FLAG_IMAGE flag.fits\n'
			flag_type = 'FLAG_TYPE OR\n'
			phot_apertures = 'PHOT_APERTURES 10\n'
			phot_autoparams = 'PHOT_AUTOPARAMS 2.5 3.5\n'
			phot_petroparams = 'PHOT_PETROPARAMS 2.0 3.5\n'
			phot_autoapers = 'PHOT_AUTOAPERS 0.0 0.0\n'
			satur_level = 'SATUR_LEVEL ' + str(Skopeus.get_peakmax()) + '\n'
			satur_key = 'SATUR_KEY SATURATE\n'
			mag_zeropoint = 'MAG_ZEROPOINT 0.0\n'
			mag_gamma = 'MAG_GAMMA 4.0\n'
			gain = 'GAIN ' + str(Skopeus.get_gain()) + '\n'
			gain_key = 'GAIN_KEY GAIN\n'
			pixel_scale = 'PIXEL_SCALE ' + str(Skopeus.get_pixel_scale()) + '\n'
			seeing_fwhm = 'SEEING_FWHM ' + str(Skopeus.get_seeing()) + '\n'
			starnnw_name = 'STARNNW_NAME default.nnw\n'
			back_type = 'BACK_TYPE AUTO\n'
			back_value = 'BACK_VALUE 0.0\n'
			back_size = 'BACK_SIZE 64\n'
			back_filtersize = 'BACK_FILTERSIZE 3\n'
			checkimage_type = 'CHECKIMAGE_TYPE NONE\n'
			checkimage_name = 'CHECKIMAGE_NAME check.fits\n'
			memory_objstack = 'MEMORY_OBJSTACK 3000\n'
			memory_pixstack = 'MEMORY_PIXSTACK 300000\n'
			memory_bufsize = 'MEMORY_BUFSIZE 1024\n'
			assoc_name = 'ASSOC_NAME sky.list\n'
			assoc_data = 'ASSOC_DATA 2 3 4\n'
			assoc_params = 'ASSOC_PARAMS 2 3 4\n'
			assoc_radius = 'ASSOC_RADIUS 2.0\n'
			assoc_type = 'ASSOC_TYPE NEAREST\n'
			assocselec_type = 'ASSOCSELEC_TYPE MATCHED\n'
			verbose_type = 'VERBOSE_TYPE QUIET\n'
			header_suffix = 'HEADER_SUFFIX .head\n'
			write_xml = 'WRITE_XML N\n'
			xml_name = 'XML_NAME sex.xml\n'
			xsl_url = 'XSL_URL file:///usr/local/share/sextractor/sextractor.xsl\n'
			lns = [catalog_name, catalog_type, parameters_name, detect_type, detect_minarea, detect_thresh, analysis_thresh, detection_filter, detection_filter_name, deblend_nthresh, deblend_mincont, clean, clean_param, weight_type, weight_image, flag_image, flag_type, phot_apertures, phot_autoparams, phot_petroparams, phot_autoapers,satur_level, satur_key, mag_zeropoint, mag_gamma, gain, gain_key, pixel_scale, seeing_fwhm, starnnw_name, back_type, back_value, back_size, back_filtersize, checkimage_type, checkimage_name, memory_objstack, memory_pixstack, memory_bufsize, assoc_name, assoc_data, assoc_params, assoc_radius, assoc_type, assocselec_type, verbose_type, header_suffix, write_xml, xml_name,xsl_url]
			default_sex = open(default_sex_path, 'w')
			for ln in lns:
				default_sex.write(ln)
			default_sex.close()

	@staticmethod
	def solve_field(frame_name, frame_directory, field, verbose=False):
		'''This function obtains and applies the WCS solution to a FITS image via local instance of Astrometry.net code.

		:p frame_name [str] - the FITS image name
		:p frame_directory [str] - the FITS image directory
		:p field [str] - the image survey field
		:p verbose [bool] - a toggle for printing Astrometry.net output to the terminal (default is False)

		:r - nothing
		'''

		# change to the frame directory
		os.chdir(frame_directory)

		# set up the source extractor
		Metretes.setup_source_extractor(frame_directory)

		# define the astrometry output files
		file_str = frame_name[:-5]
		file_axy = file_str + '.axy'
		file_corr = file_str + '.corr'
		file_match = file_str + '.match'
		file_new = file_str + '.new'
		file_rdls = file_str + '.rdls'
		file_solved = file_str + '.solved'
		file_wcs = file_str + '.wcs'
		file_xyls = file_str + '-indx.xyls'

		# define the search parameters
		solve_ra = float(Skopeus.get_field_position(field)[0])
		solve_dec = float(Skopeus.get_field_position(field)[1])
		solve_radius = Configuration.SOLVE_RADIUS

		# run the terminal command
		if verbose:
			subprocess.run(['solve-field', '--use-source-extractor', '--source-extractor-path', '/usr/bin/source-extractor', '--no-plots', frame_name, '--ra', str(solve_ra), '--dec', str(solve_dec), '--radius', str(solve_radius)])
		else:
			subprocess.run(['solve-field', '--use-source-extractor', '--source-extractor-path', '/usr/bin/source-extractor', '--no-plots', frame_name, '--ra', str(solve_ra), '--dec', str(solve_dec), '--radius', str(solve_radius)], stdout=subprocess.DEVNULL)

		# remove unnecessary output files
		subprocess.run(['rm', file_axy])
		subprocess.run(['rm', file_corr])
		subprocess.run(['rm', file_match])
		subprocess.run(['rm', file_rdls])
		subprocess.run(['rm', file_solved])
		subprocess.run(['rm', file_wcs])
		subprocess.run(['rm', file_xyls])
		subprocess.run(['mv', file_new, frame_name])

	@staticmethod
	def subtract_background(raw_frame_data, raw_frame_header, img_ini_path, img_res_path):
		'''This function measures and subtracts the background signal from a FITS image data array.

		:p raw_frame_data [numpy.ndarray] - the raw image data array
		:p raw_frame_header [astropy.io.fits.header] - the raw image header
		:p img_ini_path [str] - the initial background plot path
		:p img_res_path [str] - the residual background plot path

		:r clean_frame_data [numpy.ndarray] - the cleaned image data array
		:r clean_frame_header [astropy.io.fits.header] - the updated image header
		:r bkg2d_ini_md [float] - the initial median (mesh) [ADU]
		:r bkg2d_ini_sd [float] - the initial error (mesh) [ADU]
		:r bkgsc_ini_mn [float] - the initial mean (clip) [ADU]
		:r bkgsc_ini_md [float] - the initial median (clip) [ADU]
		:r bkgsc_ini_sd [float] - the initial error (clip) [ADU]
		:r bkg2d_res_md [float] - the residual median (mesh) [ADU]
		:r bkg2d_res_sd [float] - the residual error (mesh) [ADU]
		:r bkgsc_res_mn [float] - the residual mean (clip) [ADU]
		:r bkgsc_res_md [float] - the residual median (clip) [ADU]
		:r bkgsc_res_sd [float] - the residual error (clip) [ADU]
		'''

		# measure the initial background
		bkg2d_ini_data, bkg2d_ini_md, bkg2d_ini_sd, bkgsc_ini_mn, bkgsc_ini_md, bkgsc_ini_sd = Metretes.measure_background(raw_frame_data)

		# save the initial background plot
		if not os.path.exists(img_ini_path):
			Grapheus.draw_field_diagram(bkg2d_ini_data, None, None, img_ini_path)

		# subtract the initial background
		if Configuration.BACKGROUND_METHOD == 'flat':
			clean_frame_data = raw_frame_data - bkg2d_ini_md
		elif Configuration.BACKGROUND_METHOD == '2d':
			clean_frame_data = raw_frame_data - bkg2d_ini_data

		# measure the residual background
		bkg2d_res_data, bkg2d_res_md, bkg2d_res_sd, bkgsc_res_mn, bkgsc_res_md, bkgsc_res_sd = Metretes.measure_background(clean_frame_data)

		# save the residual background plot
		if not os.path.exists(img_res_path):
			Grapheus.draw_field_diagram(bkg2d_res_data, None, None, img_res_path)

		# update the header
		clean_frame_header = raw_frame_header
		clean_frame_header['BKGSB'] = ('median-' + Configuration.BACKGROUND_METHOD, 'Method')
		clean_frame_header['BKGTH'] = (Configuration.SIGMA_BACKGROUND, 'Threshold [sigma]')
		clean_frame_header['BKGMD'] = (bkg2d_ini_md, 'Median [ADU]')
		clean_frame_header['BKGSD'] = (bkg2d_ini_sd, 'Median RMS [ADU]')
		clean_frame_header['RESMD'] = (bkg2d_res_md, 'Residual median [ADU]')
		clean_frame_header['RESSD'] = (bkg2d_res_sd, 'Residual RMS [ADU]')

		return clean_frame_data, clean_frame_header, bkg2d_ini_md, bkg2d_ini_sd, bkgsc_ini_mn, bkgsc_ini_md, bkgsc_ini_sd, bkg2d_res_md, bkg2d_res_sd, bkgsc_res_mn, bkgsc_res_md, bkgsc_res_sd		

	@staticmethod
	def subtract_dark(raw_frame_data, raw_frame_header, dark_frame_path):
		'''This function corrects an image data array with a master dark data array, saves the file, and returns both the corrected data array and header.

		:p raw_frame_data [numpy.ndarray] - the raw image data array
		:p raw_frame_header [astropy.io.fits.header] - the raw image header
		:p dark_frame_path [str] - the master dark image path

		:r clean_frame_data [numpy.ndarray] - the cleaned image data array
		:r clean_frame_header [astropy.io.fits.header] - the updated image header
		'''

		# open the master dark frame
		dark_frame = fits.open(dark_frame_path)
		dark_frame_data = dark_frame[0].data
		dark_frame_header = dark_frame[0].header

		# subtract the dark
		clean_frame_data = raw_frame_data - dark_frame_data		

		# update the header
		clean_frame_header = raw_frame_header
		clean_frame_header['DRKSB'] = ('yes', 'Status: dark corrected')

		return clean_frame_data, clean_frame_header