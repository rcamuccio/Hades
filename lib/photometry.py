from config import Configuration
from lib.calculator import Calculator
from lib.observatory import Observatory
from lib.plot import Plot
from lib.survey import Survey
from lib.utility import Utility

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

class Photometry:

	@staticmethod
	def align_frame(target_frame_path, reference_frame_path):
		'''This function aligns a target frame at one epoch with a reference frame at another epoch.

		:parameter target_frame_path - The path of the target frame
		:parameter reference_frame_path - The path of the reference frame

		:return - Nothing is returned, but the aligned frame is re-written to the target frame path
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
		align_frame_data, align_frame_footprint = reproject_interp(target_frame_hdu, reference_frame_header)

		# correct bad pixels
		align_frame_data, target_frame_header = Photometry.correct_bad_pixels(align_frame_data, target_frame_header, detect_bad_pixels=False)

		# set the aligned frame data type
		align_frame_data = align_frame_data.astype(Configuration.DATA_TYPE)

		# modify the aligned frame header
		align_frame_header = target_frame_header
		align_frame_header['ALNIM'] = ('yes', 'Status: image aligned')

		# convert the aligned frame data to an HDU object
		align_frame_hdu = fits.PrimaryHDU(align_frame_data, header=align_frame_header)

		# save the aligned frame HDU object
		align_frame_hdu.writeto(target_frame_path, overwrite=True)

	@staticmethod
	def clean_raw_frames(date, field):
		'''This function cleans a series of raw data frames for a given date and field.

		:parameter date - The date of the raw frames
		:parameter field - The field ID of the raw frames

		:return frame_table - A table with the statistics of each raw data frame
		'''

		print('Cleaning available raw frames (' + field + ', ' + date + ')')

		# set the input and output directories
		raw_frame_directory = Configuration.INPUT_DATA_DIRECTORY + 'raw/' + date + '/' + field + '/'
		clean_frame_directory = Configuration.OUTPUT_DATA_DIRECTORY + 'clean/' + date + '/' + field + '/'

		# create the output directory
		if not os.path.exists(clean_frame_directory):
			os.makedirs(clean_frame_directory)

		# set the calibration frame paths
		dark_frame_path = Configuration.OUTPUT_DATA_DIRECTORY + 'dark/' + date + '/mdl_' + date + Configuration.FILE_EXTENSION
		flat_frame_path = Configuration.OUTPUT_DATA_DIRECTORY + 'flat/' + date + '/ff_' + date + Configuration.FILE_EXTENSION

		# generate a frame table
		frame_table_path = clean_frame_directory + 'stat_' + field + '_' + date + '.dat'
		table_columns = ('frm', 'dte', 'jd', 'tex', 'bmd', 'bsd', 'rmn', 'rmd', 'rsd')
		table_datatypes = (str, str, float, float, float, float, float, float, float)
		frame_table = Table(names=table_columns, dtype=table_datatypes)

		# load the raw frames to clean
		os.chdir(raw_frame_directory)
		raw_frame_list = sorted(glob.glob('*.fits'))
		num_raw_frames = len(raw_frame_list)
		frame_idx = 1
		print('\tNumber of raw frames found:', num_raw_frames)

		# define the reference frame for alignment
		reference_frame = 'cln_' + raw_frame_list[0]
		reference_frame_path = clean_frame_directory + reference_frame

		# loop through each frame to clean
		for frm in range(num_raw_frames):
			os.chdir(raw_frame_directory)

			# set the frame name and path
			raw_frame_file = raw_frame_list[frm]
			raw_frame_name = raw_frame_file.split('.')[0] + '.' + raw_frame_file.split('.')[1]

			clean_frame = 'cln_' + raw_frame_file
			clean_frame_name = 'cln_' + raw_frame_name
			clean_frame_path = clean_frame_directory + clean_frame

			raw_bkg_img_path = clean_frame_directory + 'bkg_' + raw_frame_name + '.png'
			clean_bkg_img_path = clean_frame_directory + 'res_' + raw_frame_name + '.png'

			# check if cleaned frame already exists
			if not os.path.exists(clean_frame_path):
				print('\t(' + str(frame_idx) + '/' + str(num_raw_frames) + ') Cleaning frame', raw_frame_file)

				# open the raw frame
				print('\t\tOpening frame')
				raw_frame = fits.open(raw_frame_file)
				raw_frame_data = raw_frame[0].data
				raw_frame_header = raw_frame[0].header

				exptime = raw_frame_header['EXPTIME']
				time = Time(raw_frame_header['DATE'], format='isot', scale='utc')
				jd = time.jd

				# remove the overscan region
				print('\t\tRemoving overscan')
				clean_frame_data, clean_frame_header = Photometry.remove_overscan(raw_frame_data, raw_frame_header)

				# subtract the master dark
				print('\t\tSubtracting dark')
				clean_frame_data, clean_frame_header = Photometry.subtract_dark(clean_frame_data, clean_frame_header, dark_frame_path)

				# divide the flatfield
				print('\t\tDividing flat')
				clean_frame_data, clean_frame_header = Photometry.divide_flat(clean_frame_data, clean_frame_header, flat_frame_path)

				# subtract the background
				print('\t\tSubtracting background')
				clean_frame_data, clean_frame_header, bkg2d_ini_md, bkg2d_ini_sd, bkgsc_ini_mn, bkgsc_ini_md, bkgsc_ini_sd, bkg2d_res_md, bkg2d_res_sd, bkgsc_res_mn, bkgsc_res_md, bkgsc_res_sd = Photometry.subtract_background(clean_frame_data, clean_frame_header, raw_bkg_img_path, clean_bkg_img_path)

				# solve the field
				print('\t\tSolving field')
				field_id = field.split('_')[1]
				fits.writeto(clean_frame_path, clean_frame_data, clean_frame_header)
				Photometry.solve_field(clean_frame, clean_frame_directory, field_id, verbose=Configuration.VERBOSE)

				# align the frames to the first
				if frm > 0:
					print('\t\tAligning to', reference_frame)
					Photometry.align_frame(clean_frame_path, reference_frame_path)

				# add a row to the frame table
				table_row = (raw_frame_file, date, jd, exptime, bkg2d_ini_md, bkg2d_ini_sd, bkgsc_res_mn, bkgsc_res_md, bkgsc_res_sd)
				frame_table.add_row(table_row)

			else:
				print('\t(' + str(frame_idx) + '/' + str(num_raw_frames) + ') Cleaned frame', raw_frame_file, 'already exists')

			frame_idx += 1

		print()

		if os.path.exists(frame_table_path):
			frame_table = ascii.read(frame_table_path, format='fixed_width', delimiter='|')
		else:
			ascii.write(frame_table, frame_table_path, format='fixed_width')

		return frame_table

	@staticmethod
	def correct_bad_pixels(frame_data, frame_header, detect_bad_pixels=True):
		'''This function corrects the pixels (nans and infs) of a data frame, and detects the individual bad pixels if desired.

		:parameter frame_data - The image data array
		:parameter frame_header - The image header
		:parameter detect_bad_pixels - A toggle for determining if indidivual bad pixels (nans and infs) are detected

		:return frame data - The corrected image data array
		:return frame_header - The corrected image header
		'''

		# set the array dimensions
		try:
			overscan_present = frame_header['OVSRM']
			overscan_toggle = True
			x_axis = Observatory.axis_x_raw()
			y_axis = Observatory.axis_y_raw()

		except KeyError:
			overscan_toggle = False
			x_axis = Observatory.axis_x()
			y_axis = Observatory.axis_y()

		# detect the individual nans and infs present
		if detect_bad_pixels:
			nan_idx = 0
			inf_idx = 0
			for y in range(0, y_axis):
				for x in range(0, x_axis):
					if np.isnan(frame_data[y][x]):
						nan_idx += 1
					elif np.isinf(frame_data[y][x]):
						inf_idx += 1

		# calculate the background statistics
		bkg_mn, bkg_md, bkg_sd = sigma_clipped_stats(frame_data, sigma=Configuration.SIG_BKG)

		# create a holder array for the bad pixels
		bad_pixels = np.isnan(frame_data) + np.isinf(frame_data)
		
		# replace the bad pixels with the background mean
		frame_data[bad_pixels] = bkg_mn

		# update the header
		frame_header['BPCOR'] = ('yes', 'Status: pixels corrected')

		return frame_data, frame_header

	@staticmethod
	def difference_frames(field, date):
		'''This function creates a series of differenced images from a reference image and a series of aligned images at different epochs.

		'''

		print('This is where differencing would go... IF I HAD ONE')

		# grab the frames to subtract
		frame_directory = Configuration.OUTPUT_DATA_DIRECTORY + 'clean/' + date + '/' + field + '/'
		os.chdir(frame_directory)
		frame_list = sorted(glob.glob('cln*' + Configuration.FILE_EXTENSION))
		num_frames = len(frame_list)

		# define the reference frame for differencing
		reference_frame = 'stack_' + field + '_' + date + Configuration.FILE_EXTENSION
		reference_frame_path = frame_directory + reference_frame

		# prepare the oisdifference.c file for differencing

		# loop through each frame to subtract
		for frm in range(num_frames):
			os.chdir(frame_directory)

			# set the frame name and path
			raw_frame_file = frame_list[frm]
			dif_frame_file = 'dif_' + raw_frame_file
			dif_frame_path = frame_directory + 'dif_' + raw_frame_file

			# check to see if the difference file already exists

			# run the differencing algorithm


	@staticmethod
	def divide_flat(raw_frame_data, raw_frame_header, flatfield_frame_path):
		'''This function corrects an image data array with a normalized flatfield data array, saves the file, and returns both the corrected data array and header.

		:parameter raw_frame_data - The image data array
		:parameter raw_frame_header - The image header
		:parameter flatfield_frame_path - The path to the flatfield frame

		:return clean_frame_data - The corrected image data array
		:return clean_frame_header - The corrected image header
		'''

		flatfield_frame = fits.open(flatfield_frame_path)
		flatfield_data = flatfield_frame[0].data
		flatfield_header = flatfield_frame[0].header

		clean_frame_data = raw_frame_data / flatfield_data
		clean_frame_header = raw_frame_header
		clean_frame_header['FLTDV'] = ('yes', 'Status: flat divided')

		return clean_frame_data, clean_frame_header

	@staticmethod
	def extract_sources(frame_data, frame_header, table_path):
		'''This function extracts sources from an image data array and returns a source table.

		:parameter frame_data - The image data array
		:parameter frame_header - The image header
		:parameter table_path - The path to the source table

		:return table - The source table
		'''

		if not os.path.exists(table_path):
			print('Extracting sources')
			wcs = WCS(frame_header)
			
			mask, boxes = Photometry.make_mask(frame_data)
			
			frame_mn, frame_md, frame_sd = sigma_clipped_stats(frame_data, mask=mask, sigma=Configuration.THRESHOLD)
			
			daofind = DAOStarFinder(fwhm=Configuration.FWHM, threshold=Configuration.SIG_SRC * frame_sd)
			table = daofind(frame_data, mask=mask)

			ra_list = []
			dec_list = []

			for ln in table:
				x = ln['xcentroid']
				y = ln['ycentroid']

				sky_pos = pixel_to_skycoord(x, y, wcs=wcs)
				
				ra = sky_pos.ra.deg
				dec = sky_pos.dec.deg
				
				ra_list.append(ra)
				dec_list.append(dec)

			table['ra'] = ra_list
			table['dec'] = dec_list

			table.write(table_path, format=Configuration.TABLE_FORMAT)

		else:
			print('Reading existing extracted table')
			table = Table.read(table_path, format=Configuration.TABLE_FORMAT)

		return table

	@staticmethod
	def point_aperture_photometry(date, field, frame_data, frame_header, point_ra, point_dec):

		return 

	@staticmethod
	def frame_aperture_photometry(date, field, frame_data, frame_header, match_table, master_table_path, output_name=None):
		'''This function performs aperture photometry on a FITS image and returns a table of photometric measurements.

		:parameter frame_data - The image data array
		:parameter frame_header - The image header
		:parameter match_table - The matched catalog table
		:parameter master_table_path - The photometry table save path

		:return master_table - The photometry table performed at the matched catalog positions
		'''

		output_directory = Configuration.OUTPUT_DATA_DIRECTORY + 'clean/' + date + '/' + field + '/'

		img_field_path = output_directory + 'plt_' + output_name + Configuration.IMAGE_EXTENSION
		img_colormag_path = output_directory + 'cm_' + output_name + Configuration.IMAGE_EXTENSION

		if not os.path.exists(master_table_path):
			print('Performing aperture photometry')

			# get the mask
			mask, boxes = Photometry.make_mask(frame_data)

			# get the coordinates
			wcs = WCS(frame_header)

			# get the observation time
			time = Time(frame_header['DATE'], format='isot', scale='utc')
			jd = time.jd

			# get the airmass
			latitude = Observatory.latitude()
			longitude = Observatory.longitude()
			elevation = Observatory.elevation()
			location = EarthLocation(lat=Observatory.latitude(), lon=Observatory.longitude(), height=Observatory.elevation())

			field_ra = float(Survey.get_field(field)[0])
			field_dec = float(Survey.get_field(field)[1])

			altitude = Calculator.sky_to_altitude(location, time, field_ra, field_dec)
			airmass = Calculator.airmass(altitude)
			print('TEST:', latitude, longitude, elevation, field_ra, field_dec, altitude, airmass)

			# get the exposure time
			exp_time = frame_header['EXPTIME']

			# calculate the background statistics
			bkg_mn, bkg_md, bkg_sd = sigma_clipped_stats(frame_data, mask=mask, sigma=Configuration.SIG_BKG)

			# extract coordinates from the input table
			sky_coords = []
			for obj in match_table:
				ra = obj['ra']
				dec = obj['dec']
				coord = (ra, dec)
				sky_coords.append(coord)
			sky_pos = SkyCoord(sky_coords, unit='deg')
			pix_pos = skycoord_to_pixel(sky_pos, wcs=wcs)
			pix_coords = []
			for px in range(len(pix_pos[0])):
				xp = pix_pos[0][px]
				yp = pix_pos[1][px]
				coord = (xp, yp)
				pix_coords.append(coord)

			# create apertures and annuli
			apertures = CircularAperture(pix_coords, r=Configuration.APER_SIZE)
			annuli = CircularAnnulus(pix_coords, r_in=Configuration.ANNULI_INNER, r_out=Configuration.ANNULI_OUTER)
			aper_set = [apertures, annuli]
			aperture_area = apertures.area
			annulus_area = annuli.area

			# perform photometry at all positions
			phot_table = aperture_photometry(frame_data, aper_set, mask=mask, wcs=wcs)
			bkg_loc_mean = phot_table['aperture_sum_1'] / annulus_area
			bkg_loc_sum = bkg_loc_mean * aperture_area
			res_loc_sum = phot_table['aperture_sum_0'] - bkg_loc_sum

			phot_table['annulus_mean'] = bkg_loc_mean
			phot_table['aperture_bkg_sum'] = bkg_loc_sum
			phot_table['aperture_res_sum'] = res_loc_sum

			# merge the input and photometry tables
			master_table = hstack([match_table, phot_table])

			# calculate magnitudes and errors
			flux_list = []
			flux_err_list = []
			inst_mag_list = []
			inst_mag_err_list = []
			color_list = []
			delta_mag_list = []
			flux_flag_list = []
			color_flag_list = []

			for src in master_table:
				flux = src['aperture_res_sum']

				if flux < 0:
					flux_flag_list.append(True)
					flux = abs(bkg_md)
				else:
					flux_flag_list.append(False)

				flux_err = np.sqrt(flux + aperture_area * bkg_sd**2)
				inst_mag = - 2.5 * np.log10(flux) + 2.5 * np.log10(exp_time)
				inst_mag_err = (2.5 * flux_err) / (np.log(10.) * flux)

				if (src['phot_bp_n_obs'] == 0) or (src['phot_rp_n_obs'] == 0):
					color = -9.999999
					cat_mag = -9.999999
					delta_mag = -9.999999
					color_flag_list.append(True)
				else:
					b = src['phot_bp_mean_mag']
					r = src['phot_rp_mean_mag']
					color = b - r
					cat_mag = src['phot_g_mean_mag']
					delta_mag = cat_mag - inst_mag
					color_flag_list.append(False)

				flux_list.append(flux)
				flux_err_list.append(flux_err)
				inst_mag_list.append(inst_mag)
				inst_mag_err_list.append(inst_mag_err)
				color_list.append(color)
				delta_mag_list.append(delta_mag)

			master_table['flux'] = flux_list
			master_table['flux_flag'] = flux_flag_list
			master_table['flux_err'] = flux_err_list
			master_table['inst_mag'] = inst_mag_list
			master_table['inst_mag_err'] = inst_mag_err_list
			master_table['color'] = color_list
			master_table['color_flag'] = color_flag_list
			master_table['delta_mag'] = delta_mag_list

			# calculate transform and zero point
			color_list = []
			delta_list = []
			for src in master_table:
				if (src['flux_flag'] == True) or (src['color_flag'] == True):
					pass
				else:
					color = src['color']
					delta = src['delta_mag']
					color_list.append(color)
					delta_list.append(delta)
			yfit, slope, intercept, delta_slope, delta_intercept = Calculator.unweighted_fit(color_list, delta_list)

			# plot an annotated frame
			Plot.field(frame_data, img_field_path, apertures, boxes)

			# plot a color-magnitude diagram
			Plot.colormag(color_list, delta_list, yfit, img_colormag_path)

			# save the table
			if not os.path.exists(master_table_path):
				master_table.write(master_table_path, format=Configuration.TABLE_FORMAT)

		else:
			print('Reading existing master photometry table')
			master_table = Table.read(master_table_path, format=Configuration.TABLE_FORMAT)

		return master_table

	@staticmethod
	def get_overscan(frame_data):

		# horizontal overscan regions
		ovs_h_x = Observatory.tap_x()
		ovs_h_y = 2 * Observatory.overscan_y()
		ovs_h_size = u.Quantity([ovs_h_y, ovs_h_x], u.pixel)

		ovs_h_center_01 = (660, 5300)
		ovs_h_center_02 = (2160, 5300)
		ovs_h_center_03 = (3660, 5300)
		ovs_h_center_04 = (5160, 5300)
		ovs_h_center_05 = (6660, 5300)
		ovs_h_center_06 = (8160, 5300)
		ovs_h_center_07 = (9660, 5300)
		ovs_h_center_08 = (11160, 5300)
		ovs_h_center_list = [ovs_h_center_01, ovs_h_center_02, ovs_h_center_03, ovs_h_center_04, ovs_h_center_05, ovs_h_center_06, ovs_h_center_07, ovs_h_center_08]

		ovs_h_mn_list = []
		ovs_h_sd_list = []
		n_ovs_h = len(ovs_h_center_list)

		for i in range(n_ovs_h):
			ovs = i+1
			ovs_center = ovs_h_center_list[i]

			cutout = Cutout2D(frame_data, ovs_center, ovs_h_size)
			bbox = cutout.bbox_original
			x1 = bbox[1][0]
			x2 = bbox[1][1]
			y1 = bbox[0][0]
			y2 = bbox[0][1]

			box = RectangularAperture(ovs_center, ovs_h_x, ovs_h_y, theta=0.)

			cutout_data = cutout.data
			cutout_min = np.min(cutout_data)
			cutout_max = np.max(cutout_data)
			cutout_mn, cutout_md, cutout_sd = sigma_clipped_stats(cutout_data, sigma=Configuration.SIG_BKG)

			ovs_h_mn_list.append(float(np.around(cutout_mn, decimals=2)))
			ovs_h_sd_list.append(float(np.around(cutout_sd, decimals=2)))

		# vertical overscan regions
		ovs_v_x = Observatory.overscan_x()
		ovs_v_y = Observatory.axis_y() + ovs_h_y
		ovs_v_size = u.Quantity([ovs_v_y, ovs_v_x], u.pixel)

		ovs_v_center_01 = (1410, 5300)
		ovs_v_center_02 = (2910, 5300)
		ovs_v_center_03 = (4410, 5300)
		ovs_v_center_04 = (5910, 5300)
		ovs_v_center_05 = (7410, 5300)
		ovs_v_center_06 = (8910, 5300)
		ovs_v_center_07 = (10410, 5300)
		ovs_v_center_08 = (11910, 5300)
		ovs_v_center_list = [ovs_v_center_01, ovs_v_center_02, ovs_v_center_03, ovs_v_center_04, ovs_v_center_05, ovs_v_center_06, ovs_v_center_07, ovs_v_center_08]

		ovs_v_mn_list = []
		ovs_v_sd_list = []
		n_ovs_v = len(ovs_v_center_list)

		for i in range(n_ovs_v):
			ovs = i+1
			ovs_center = ovs_v_center_list[i]

			cutout = Cutout2D(frame_data, ovs_center, ovs_v_size)
			bbox = cutout.bbox_original
			x1 = bbox[1][0]
			x2 = bbox[1][1]
			y1 = bbox[0][0]
			y2 = bbox[0][1]

			box = RectangularAperture(ovs_center, ovs_v_x, ovs_v_y, theta=0.)

			cutout_data = cutout.data
			cutout_min = np.min(cutout_data)
			cutout_max = np.max(cutout_data)
			cutout_mn, cutout_md, cutout_sd = sigma_clipped_stats(cutout_data, sigma=Configuration.SIG_BKG)

			ovs_v_mn_list.append(float(np.around(cutout_mn, decimals=2)))
			ovs_v_sd_list.append(float(np.around(cutout_sd, decimals=2)))

		return ovs_h_mn_list, ovs_h_sd_list, ovs_v_mn_list, ovs_v_sd_list

	@staticmethod
	def get_taps(frame_data, frame_type):

		tap_x = Observatory.tap_x()
		tap_y = Observatory.tap_y()
		tap_size = u.Quantity([tap_y, tap_x], u.pixel)

		if frame_type == 'raw':
			tap_center_01 = (660, 2640)
			tap_center_02 = (660, 7960)
			tap_center_03 = (2160, 2640)
			tap_center_04 = (2160, 7960)
			tap_center_05 = (3660, 2640)
			tap_center_06 = (3660, 7960)
			tap_center_07 = (5160, 2640)
			tap_center_08 = (5160, 7960)
			tap_center_09 = (6660, 2640)
			tap_center_10 = (6660, 7960)
			tap_center_11 = (8160, 2640)
			tap_center_12 = (8160, 7960)
			tap_center_13 = (9660, 2640)
			tap_center_14 = (9660, 7960)
			tap_center_15 = (11160, 2640)
			tap_center_16 = (11160, 7960)

		elif frame_type == 'clean':
			tap_center_01 = (660, 2640)
			tap_center_02 = (660, 7920)
			tap_center_03 = (1980, 2640)
			tap_center_04 = (1980, 7920)
			tap_center_05 = (3300, 2640)
			tap_center_06 = (3300, 7920)
			tap_center_07 = (4620, 2640)
			tap_center_08 = (4620, 7920)
			tap_center_09 = (5940, 2640)
			tap_center_10 = (5940, 7920)
			tap_center_11 = (7260, 2640)
			tap_center_12 = (7260, 7920)
			tap_center_13 = (8580, 2640)
			tap_center_14 = (8580, 7920)
			tap_center_15 = (9900, 2640)
			tap_center_16 = (9900, 7920)

		tap_center_list = [tap_center_01, tap_center_02, tap_center_03, tap_center_04, tap_center_05, tap_center_06, tap_center_07, tap_center_08, tap_center_09, tap_center_10, tap_center_11, tap_center_12, tap_center_13, tap_center_14, tap_center_15, tap_center_16]

		ntaps = len(tap_center_list)

		tap_mn_list = []
		tap_sd_list = []
		tap_box_list = []

		for i in range(ntaps):
			tap = i+1
			tap_center = tap_center_list[i]

			cutout = Cutout2D(frame_data, tap_center, tap_size)
			bbox = cutout.bbox_original
			x1 = bbox[1][0]
			x2 = bbox[1][1]
			y1 = bbox[0][0]
			y2 = bbox[0][1]
			box = RectangularAperture(tap_center, tap_x, tap_y, theta=0.)
			tap_box_list.append(box)

			cutout_data = cutout.data
			cutout_min = np.min(cutout_data)
			cutout_max = np.max(cutout_data)
			cutout_mn, cutout_md, cutout_sd = sigma_clipped_stats(cutout_data, sigma=Configuration.SIG_BKG)

			tap_mn_list.append(float(np.around(cutout_mn, decimals=2)))
			tap_sd_list.append(float(np.around(cutout_sd, decimals=2)))

			# calculate masked statistics
			#sigma_clip = SigmaClip(sigma=Configuration.SIG_BKG, maxiters=10)
			#threshold = detect_threshold(cutout_data, nsigma=2.0, sigma_clip=sigma_clip)
			#segment_img = detect_sources(cutout_data, threshold, npixels=10)
			#footprint = circular_footprint(radius=10)
			#source_mask = segment_img.make_source_mask(footprint=footprint)
			#cutout_mask_mean, cutout_mask_median, cutout_mask_std = sigma_clipped_stats(cutout_data, sigma=Configuration.SIG_BKG, mask=source_mask)

		return tap_mn_list, tap_sd_list, tap_box_list	

	@staticmethod
	def make_dark(date, dark_type):

		# set the required directories
		dark_frame_dir = Configuration.INPUT_DATA_DIRECTORY + 'dark/' + date + '/'
		master_dark_frame_dir = Configuration.OUTPUT_DATA_DIRECTORY + 'dark/' + date + '/'

		# set the required file paths
		if dark_type == 'flat':
			master_dark_frame_path = master_dark_frame_dir + 'mdf_' + date + Configuration.FILE_EXTENSION
			dark_table_path = master_dark_frame_dir + 'stat_mdf_' + date + '.dat'
		elif dark_type == 'light':
			master_dark_frame_path = master_dark_frame_dir + 'mdl_' + date + Configuration.FILE_EXTENSION
			dark_table_path = master_dark_frame_dir + 'stat_mdl_' + date + '.dat'

		if not os.path.exists(master_dark_frame_path):
			if not os.path.exists(master_dark_frame_dir):
				os.mkdir(master_dark_frame_dir)

			dark_frame_list = []
			dark_ccddata_list = []
			dark_keep = 0

			# load the dark frames to combine
			os.chdir(dark_frame_dir)

			print('Generating a master dark (' + dark_type + ') for', date)

			if dark_type == 'flat':
				search_str = 'dark_' + str(int(Observatory.exposure_time('dark-f'))) + 's*' + Configuration.FILE_EXTENSION
			elif dark_type == 'light':
				search_str = 'dark_' + str(int(Observatory.exposure_time('dark-l'))) + 's*' + Configuration.FILE_EXTENSION
			
			dark_frame_list = sorted(glob.glob(search_str))

			num_dark_frames = len(dark_frame_list)
			print('    Number of dark frames found:', num_dark_frames)

			# generate a dark table
			tbl_cols = ('frm', 'dte', 'jd', 'tex', 'mn', 'md', 'sd', 'kp')
			tbl_dtype = (str, str, float, float, float, float, float, str)
			dark_table = Table(names=tbl_cols, dtype=tbl_dtype)

			for frm in range(num_dark_frames):
				# open the dark frame
				dark_frame = fits.open(dark_frame_list[frm])
				dark_frame_data = dark_frame[0].data
				dark_frame_header = dark_frame[0].header

				# grab the time information
				time = Time(dark_frame_header['DATE'], format='isot', scale='utc')
				jd = time.jd

				# calculate background statistics
				dark_mn, dark_md, dark_sd = sigma_clipped_stats(dark_frame_data, sigma=Configuration.SIG_BKG)

				# remove the overscan from the dark
				dark_frame_data, dark_frame_header = Photometry.remove_overscan(dark_frame_data, dark_frame_header)

				# add ccddata object to list
				dark_ccddata = ccdproc.CCDData(dark_frame_data, unit='adu')
				dark_ccddata_list.append(dark_ccddata)

				# count the kept dark
				dark_keep += 1
				dark_keep_flag = 'Y'
				print('    (++) Frame:', dark_frame_list[frm])

				# add a row to the dark table
				tbl_row = (dark_frame_list[frm], date, jd, 300., dark_mn, dark_md, dark_sd, dark_keep_flag)
				dark_table.add_row(tbl_row)

			print('    Number of darks to combine:', dark_keep)

			if num_dark_frames > 0:
				# save the table
				if os.path.exists(dark_table_path):
					dark_table = ascii.read(dark_table_path, format='fixed_width', delimiter='|')
				else:
					ascii.write(dark_table, dark_table_path, format='fixed_width')

				# create a master dark
				print('    Generating a master dark for', date)
				master_dark_ccddata = ccdproc.combine(dark_ccddata_list, method=Configuration.COMBINE_METHOD, unit=Configuration.UNIT, mem_limit=Configuration.MEMORY_LIMIT, dtype=Configuration.DATA_TYPE)

				# convert the ccddata object to a numpy array
				master_dark_data = np.asarray(master_dark_ccddata)

				# read the master dark statistics
				master_dark_mn, master_dark_md, master_dark_sd = sigma_clipped_stats(master_dark_data, sigma=Configuration.SIG_BKG)

				# modify the master dark header
				master_dark_hdu = fits.PrimaryHDU(master_dark_data)
				master_dark_header = master_dark_hdu.header

				master_dark_header['SIMPLE'] = ('T', 'Conform to FITS standard')
				master_dark_header['BITPIX'] = (32, 'Number of bits per data pixel')
				master_dark_header['NAXIS'] = (2, 'Number of axes')
				master_dark_header['NAXIS1'] = (10560, 'Image width')
				master_dark_header['NAXIS2'] = (10560, 'Image height')
				master_dark_header['BZERO'] = (32768, 'Offset for unsigned short')
				master_dark_header['BSCALE'] = (1, 'Default scaling factor')
				master_dark_header['DATE'] = (date, 'Date of combined darks')
				master_dark_header['EXPTIME'] = (300, 'Exposure time')
				master_dark_header['IMAGETYP'] = ('DARK', 'Type of image')
				master_dark_header['DRKCB'] = ('yes', 'Status: darks combined')
				master_dark_header['DRKNM'] = (dark_keep, 'Number of combined dark frames')
				master_dark_header['MDMD'] = (master_dark_md, 'Master dark median [ADU]')
				master_dark_header['MDMN'] = (master_dark_mn, 'Master dark mean [ADU]')
				master_dark_header['MDSD'] = (master_dark_sd, 'Master dark standard deviation [ADU]')
				master_dark_header['OVSRM'] = ('yes', 'Status: overscan removed')

				fits.writeto(master_dark_frame_path, master_dark_data, master_dark_header)

				print('    Master dark generated for', date)
				print()

			else:
				print('    No dark frames available')
				print()

		else:
			print('Master dark already exists for', date)
			print()

	@staticmethod
	def make_flat(date):

		# set the required directories
		flat_frame_dir = Configuration.INPUT_DATA_DIRECTORY + 'flat/' + date + '/'
		flatfield_frame_dir = Configuration.OUTPUT_DATA_DIRECTORY + 'flat/' + date + '/'
		master_dark_frame_dir = Configuration.OUTPUT_DATA_DIRECTORY + 'dark/' + date + '/'

		# set the required file paths
		flatfield_frame_path = flatfield_frame_dir + 'ff_' + date + Configuration.FILE_EXTENSION
		master_dark_frame_path = master_dark_frame_dir + 'mdf_' + date + Configuration.FILE_EXTENSION
		flat_table_path = flatfield_frame_dir + 'stat_ff_' + date + '.dat'

		if not os.path.exists(flatfield_frame_path):
			print('Generating a normalized flatfield for', date)

			if not os.path.exists(flatfield_frame_dir):
				os.mkdir(flatfield_frame_dir)

			# load the master dark
			master_dark_ccddata = ccdproc.fits_ccddata_reader(master_dark_frame_path, unit=Configuration.UNIT)

			dark_exp = Observatory.exposure_time('dark-f') * u.second
			flat_exp = Observatory.exposure_time('flat') * u.second

			flat_frame_list = []
			flat_ccddata_list = []
			flat_keep = 0

			# load the flat frames to combine
			os.chdir(flat_frame_dir)

			search_str = 'flat_' + str(int(Observatory.exposure_time('flat'))) + 's*' + Configuration.FILE_EXTENSION

			flat_frame_list = sorted(glob.glob(search_str))

			num_flat_frames = len(flat_frame_list)
			print('    Number of flat frames found:', num_flat_frames)

			# generate a flat table
			tbl_cols = ('frm', 'dte', 'jd', 'tex', 'mn', 'md', 'sd', 'kp')
			tbl_dtype = (str, str, float, float, float, float, float, str)
			flat_table = Table(names=tbl_cols, dtype=tbl_dtype)

			for frm in range(num_flat_frames):
				# open the flat frame
				flat_frame = fits.open(flat_frame_list[frm])
				flat_frame_data = flat_frame[0].data
				flat_frame_header = flat_frame[0].header

				# grab the time information
				time = Time(flat_frame_header['DATE'], format='isot', scale='utc')
				jd = time.jd

				# calculate background statistics
				flat_mn, flat_md, flat_sd = sigma_clipped_stats(flat_frame_data, sigma=Configuration.SIG_BKG)

				# filter out flats with backgrounds outside the limits
				if (flat_mn > 3000.) and (flat_mn < 40000.):
					# remove the overscan from the flat
					flat_frame_data, flat_frame_header = Photometry.remove_overscan(flat_frame_data, flat_frame_header)

					# subtract the master dark from the flat
					flat_ccddata = ccdproc.CCDData(flat_frame_data, unit=Configuration.UNIT)
					flat_ccddata = ccdproc.subtract_dark(flat_ccddata, master_dark_ccddata, dark_exposure=dark_exp, data_exposure=flat_exp)

					# add ccddata object to list
					flat_ccddata_list.append(flat_ccddata)

					# count the kept flat
					flat_keep += 1
					flat_keep_flag = 'Y'
					print('    (+) Frame:', flat_frame_list[frm])

				else:
					flat_keep_flag = 'N'
					print('    (-) Frame:', flat_frame_list[frm])

				# add a row to the flat table
				tbl_row = (flat_frame_list[frm], date, jd, 5., flat_mn, flat_md, flat_sd, flat_keep_flag)
				flat_table.add_row(tbl_row)

			print('    Number of flat frames to combine:', flat_keep)

			if num_flat_frames > 0:
				# save the table
				if os.path.exists(flat_table_path):
					flat_table = ascii.read(flat_table_path, format='fixed_width', delimiter='|')
				else:
					ascii.write(flat_table, flat_table_path, format='fixed_width')

				# create a combined flat
				print('    Creating a combined flat')
				combined_flat_ccddata = ccdproc.combine(flat_ccddata_list, method=Configuration.COMBINE_METHOD, unit=Configuration.UNIT, mem_limit=Configuration.MEMORY_LIMIT, dtype=Configuration.DATA_TYPE)
				combined_flat_data = np.asarray(combined_flat_ccddata)

				# read the combined flat statistics
				combined_flat_mn, combined_flat_md, combined_flat_sd = sigma_clipped_stats(combined_flat_data, sigma=Configuration.SIG_BKG)

				# create a normalized flatfield
				print('    Creating a normalized flatfield')
				flatfield_data = combined_flat_data / combined_flat_mn
				flatfield_data = flatfield_data.astype(Configuration.DATA_TYPE)

				# read the normalized flatfield statistics
				flatfield_mn, flatfield_md, flatfield_sd = sigma_clipped_stats(flatfield_data, sigma=Configuration.SIG_BKG)

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
				flatfield_header['DATE'] = (date, 'Date of combined darks')
				flatfield_header['EXPTIME'] = (5, 'Exposure time')
				flatfield_header['IMAGETYP'] = ('FLAT', 'Type of image')
				flatfield_header['CFMD'] = (combined_flat_md, 'Combined flat median [ADU]')
				flatfield_header['CFMN'] = (combined_flat_mn, 'Combined flat mean [ADU]')
				flatfield_header['CFSD'] = (combined_flat_sd, 'Combined flat standard deviation [ADU]')
				flatfield_header['FFMD'] = (flatfield_md, 'Flatfield median [ADU]')
				flatfield_header['FFMN'] = (flatfield_mn, 'Flatfield mean [ADU]')
				flatfield_header['FFSD'] = (flatfield_sd, 'Flatfield standard deviation [ADU]')
				flatfield_header['FTCB'] = ('yes', 'Status: flats combined')
				flatfield_header['FTNM'] = (flat_keep, 'Number of combined flat frames')
				flatfield_header['OVSRM'] = ('yes', 'Status: overscan removed')

				fits.writeto(flatfield_frame_path, flatfield_data, flatfield_header)

				print('    Normalized flatfield generated for', date)
				print()

			else:
				print('    No flat frames available')
				print()

		else:
			print('Normalized flatfield already exists for', date)
			print()

	@staticmethod
	def make_mask(frame_data):

		# generate mask holder
		frame_shape = frame_data.shape
		dy = frame_shape[0]
		dx = frame_shape[1]
		mask = np.zeros(frame_shape, dtype=bool)

		# generate edge region masks
		edge = 101 // Observatory.binning()

		# left edge mask
		mask[0:dy, 0:edge] = True
		eml_cen_x = edge // 2
		eml_cen_y = dy // 2
		eml_box = RectangularAperture((eml_cen_x, eml_cen_y), edge, dy,  theta=0.)

		# top edge mask
		mask[dy-edge:dy, 0:dx] = True		
		emt_cen_x = dx // 2
		emt_cen_y = edge // 2
		emt_box = RectangularAperture((emt_cen_x, emt_cen_y + dy - edge), dx, edge, theta=0.)

		# right edge mask
		mask[0:dy, dx-edge:dx] = True
		emr_cen_x = edge // 2
		emr_cen_y = dy // 2
		emr_box = RectangularAperture((emr_cen_x + dx - edge, emr_cen_y), edge, dy, theta=0.)

		# bottom edge mask
		mask[0:edge, 0:dx] = True
		emb_cen_x = dx // 2
		emb_cen_y = edge // 2
		emb_box = RectangularAperture((emb_cen_x, emb_cen_y), dx, edge, theta=0.)

		# generate overscan boundary masks
		w = 20
		tx = Observatory.tap_x()
		ty = Observatory.tap_y()

		# vertical bands
		mask[0:dy, tx-w:tx+w] = True
		mask[0:dy, (2*tx)-w:(2*tx)+w] = True
		mask[0:dy, (3*tx)-w:(3*tx)+w] = True
		mask[0:dy, (4*tx)-w:(4*tx)+w] = True
		mask[0:dy, (5*tx)-w:(5*tx)+w] = True
		mask[0:dy, (6*tx)-w:(6*tx)+w] = True
		mask[0:dy, (7*tx)-w:(7*tx)+w] = True

		# horizontal band
		mask[ty-w:ty+w, 0:dx] = True

		#box_list = [eml_box, emt_box, emr_box, emb_box]

		boxes = {}
		boxes['eml_box'] = eml_box
		boxes['emt_box'] = emt_box
		boxes['emr_box'] = emr_box
		boxes['emb_box'] = emb_box

		return mask, boxes

	@staticmethod
	def make_stack(date, field, frame_table):
		'''This function creates a master stack for a particular field and night. The frame table is used for additional information in the stack header.

		:parameter date - The date of the stack
		:parameter field - The field of the stack
		:parameter frame_table - The table of information for the individual frames used in the stack

		:return stack_data - The stacked image data array
		:return stack_header - The stacked image header
		'''

		stack_directory = Configuration.OUTPUT_DATA_DIRECTORY + 'clean/' + date + '/' + field + '/'
		stack_subname = 'stack_' + field + '_' + date
		stack_name = stack_subname + Configuration.FILE_EXTENSION
		stack_path = stack_directory + stack_name

		# image paths
		img_stack_path = stack_directory + 'plt_stk_' + field + '_' + date + Configuration.IMAGE_EXTENSION
		img_bkg_path = stack_directory + 'plt_bkg_' + field + '_' + date + Configuration.IMAGE_EXTENSION

		# histogram paths
		hist_stack_path = stack_directory + 'hst_stk_' + field + '_' + date + Configuration.IMAGE_EXTENSION
		hist_bkg_path = stack_directory + 'hst_bkg_' + field + '_' + date + Configuration.IMAGE_EXTENSION

		if not os.path.exists(stack_path):
			# grab the individual frames
			os.chdir(stack_directory)
			stack_list = glob.glob('cln_*.fits')
			num_stack_frames = len(stack_list)
			print('Generating', num_stack_frames, 'frames into master stack (' + field + ', ' + date + ')')

			# combine the frames
			stack_ccddata = ccdproc.combine(stack_list, method=Configuration.COMBINE_METHOD, unit=Configuration.UNIT, mem_limit=Configuration.MEMORY_LIMIT, dtype=Configuration.DATA_TYPE)
			stack_data = np.asarray(stack_ccddata)
			stack_data = stack_data.astype(Configuration.DATA_TYPE)

			if not os.path.exists(img_stack_path):
				Plot.field(stack_data, img_stack_path)

			if not os.path.exists(hist_stack_path):
				Plot.histogram(stack_data, hist_stack_path)

			# calculate the background statistics
			print('\tCalculating background statistics')
			bkg2d_data, bkg2d_md, bkg2d_sd, bkgsc_mn, bkgsc_md, bkgsc_sd = Photometry.measure_background(stack_data)

			if not os.path.exists(img_bkg_path):
				Plot.field(bkg2d_data, img_bkg_path)

			if not os.path.exists(hist_bkg_path):
				Plot.histogram(bkg2d_data, hist_bkg_path)

			# edit the stack header
			stk_hdu = fits.PrimaryHDU(stack_data)
			stack_header = stk_hdu.header
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
			stack_header['STKNM'] = (num_stack_frames, 'Number of combined light frames')
			fits.writeto(stack_path, stack_data, stack_header)

			# solve the stack frame
			print('\tSolving field')
			field_id = field.split('_')[1]
			Photometry.solve_field(stack_name, stack_directory, field_id, verbose=Configuration.VERBOSE)

			print('\tStack generated (' + field + ', ' + date + ')\n')
			stack_frame = fits.open(stack_path)
			stack_header = stack_frame[0].header

		else:
			print('Reading existing stack (' + field + ', ' + date + ')\n')
			stack_frame = fits.open(stack_path)
			stack_data = stack_frame[0].data
			stack_header = stack_frame[0].header

		return stack_data, stack_header

	@staticmethod
	def match_catalogs(source_table, query_table, match_table_path):
		'''This function matches an extracted catalog of sources with a queried catalog of sources.

		:parameter source_table - An Astropy table of extracted sources
		:parameter query_table - An Astropy table of queried sources
		:parameter match_table_path - The path of the matched table

		:return match_table - An Astropy table of matched sources
		'''

		if not os.path.exists(match_table_path):
			print('Matching catalogs')
			source_coord = SkyCoord(source_table['ra']*u.deg, source_table['dec']*u.deg)
			query_coord = SkyCoord(query_table['ra'], query_table['dec'], unit='deg')

			idx, d2d, d3d = query_coord.match_to_catalog_sky(source_coord)
			t = d2d < 5 * u.arcsec

			query_table['sep_flag'] = t
			query_table['idx'] = idx

			ct = query_table.dtype
			names, dtypes = zip(*ct.descr)
			match_table = Table(names=names, dtype=dtypes)

			for row, sep_flag in enumerate(query_table['sep_flag']):
				if sep_flag:
					match_table.add_row(query_table[row].as_void())

			match_table.write(match_table_path, format=Configuration.TABLE_FORMAT)

		else:
			print('Reading matched catalog')
			match_table = Table.read(match_table_path, format=Configuration.TABLE_FORMAT)

		return match_table

	@staticmethod
	def measure_background(frame_data):

		# calculate the two-dimensional background
		sigma_clip = SigmaClip(sigma=Configuration.SIG_SRC)
		threshold = detect_threshold(frame_data, nsigma=Configuration.SIG_BKG, sigma_clip=sigma_clip)
		segment_img = detect_sources(frame_data, threshold, npixels=Configuration.NPIXELS)
		footprint = circular_footprint(radius=Configuration.FOOTPRINT_RADIUS)
		mask = segment_img.make_source_mask(footprint=footprint)
		bkg_estimator = MedianBackground()
		bkg = Background2D(frame_data, box_size=Configuration.BOX_SIZE, mask=mask, filter_size=Configuration.FILTER_SIZE, sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
		bkg2d_data = bkg.background
		bkg2d_md = bkg.background_median
		bkg2d_sd = bkg.background_rms_median

		# calculate the sigma-clipped statistics
		bkgsc_mn, bkgsc_md, bkgsc_sd = sigma_clipped_stats(frame_data, sigma=Configuration.SIG_BKG)

		bkg2d_md = float(bkg2d_md)
		bkg2d_sd = float(bkg2d_sd)
		bkgsc_mn = float(bkgsc_mn)
		bkgsc_md = float(bkgsc_md)
		bkgsc_sd = float(bkgsc_sd)

		return bkg2d_data, bkg2d_md, bkg2d_sd, bkgsc_mn, bkgsc_md, bkgsc_sd

	@staticmethod
	def read_raw_frame(path, display=True, show_mask=True):

		# load frame
		frame = fits.open(path)
		frame_data = frame[0].data
		frame_header = frame[0].header

		# read header
		frame_x = frame_header['NAXIS1']
		frame_y = frame_header['NAXIS2']
		date = frame_header['DATE']
		time = Time(frame_header['DATE'], format='isot', scale='utc')
		jd = time.jd
		exp_time = frame_header['EXPTIME']
		imagetype = frame_header['IMAGETYP']

		# read coordinates
		wcs = WCS(frame_header)

		# read overscan
		ovs_h_mn_list, ovs_h_sd_list, ovs_v_mn_list, ovs_v_sd_list = Photometry.get_overscan(frame_data)

		# read taps
		tap_mn_list, tap_sd_list, tap_box_list = Photometry.get_taps(frame_data)

		if display:
			interval = MinMaxInterval()
			vmin, vmax = interval.get_limits(frame_data)
			vmin = 2200
			vmax = 3300
			norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LinearStretch())

			plt.figure(figsize=(10, 10))
			plt.imshow(frame_data, cmap='gray', origin='lower', norm=norm, interpolation='nearest')
			plt.xlim([0, 12000])
			plt.ylim([0, 10600])
			plt.colorbar()

			#if show_mask:
				#for ovs_box in ovs_box_list:
					#ovs_box.plot(color='magenta', ls='dashed')

				#for ovs in ovs_list:
					#plt.plot(ovs[0], ovs[1], color='magenta', marker='+', markersize=10)

				#for tap_box in tap_box_list:
					#tap_box.plot(color='red', ls='dashed')

				#for tap in tap_list:
					#plt.plot(tap[0], tap[1], color='red', marker='+', markersize=10)

				#for box in box_list:
					#box.plot(color='cyan', ls='dashed')

			save = False
			if save:
				plt.savefig('test.png', dpi=400)
			else:
				plt.show()
			plt.close()

		return frame_data, frame_header, jd, tap_mn_list, tap_sd_list, ovs_h_mn_list, ovs_h_sd_list, ovs_v_mn_list, ovs_v_sd_list

	@staticmethod
	def remove_overscan(frame_data, frame_header):

		# make the clipped data frame
		clipped_data = np.zeros((Observatory.axis_y(), Observatory.axis_x()))

		# overscan size
		ovs_x = Observatory.overscan_x()
		ovs_y = Observatory.overscan_y()

		# tap size
		tap_x = Observatory.tap_x()
		tap_y = Observatory.tap_y()

		# full tap size
		full_tap_x = tap_x + ovs_x
		full_tap_y = tap_y + ovs_y

		# move through x and y across the frame
		idx = 0
		for x in range(0, Observatory.axis_x_raw(), full_tap_x):
			idy = 0
			for y in range(0, Observatory.axis_y_raw(), full_tap_y):
				if y == 0:
					clipped_data[idy:idy+tap_y, idx:idx+tap_x] = frame_data[y:y+tap_y, x:x+tap_x]
				else:
					clipped_data[idy:idy+tap_y, idx:idx+tap_x] = frame_data[y+ovs_y:y+ovs_y+tap_y, x:x+tap_x]
				idy += tap_y
			idx += tap_x

		# change the data type
		clipped_data = np.float32(clipped_data)

		# update the header
		frame_header['OVSRM'] = ('yes', 'Status: overscan removed')

		return clipped_data, frame_header

	@staticmethod
	def setup_source_extractor(save_directory):

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
			satur_level = 'SATUR_LEVEL ' + str(Observatory.peakmax()) + '\n'
			satur_key = 'SATUR_KEY SATURATE\n'
			mag_zeropoint = 'MAG_ZEROPOINT 0.0\n'
			mag_gamma = 'MAG_GAMMA 4.0\n'
			gain = 'GAIN ' + str(Observatory.gain()) + '\n'
			gain_key = 'GAIN_KEY GAIN\n'
			pixel_scale = 'PIXEL_SCALE ' + str(Observatory.pixel_scale()) + '\n'
			seeing_fwhm = 'SEEING_FWHM ' + str(Observatory.seeing()) + '\n'
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
	def solve_field(frame_name, frame_dir, field, verbose=False):

		# change to the frame directory
		os.chdir(frame_dir)

		# set up source extractor
		Photometry.setup_source_extractor(frame_dir)

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
		ra = float(Survey.get_field(field)[0])
		dec = float(Survey.get_field(field)[1])
		radius = Configuration.SOLVE_RADIUS

		# run the command
		if verbose:
			subprocess.run(['solve-field', '--use-source-extractor', '--source-extractor-path', '/usr/bin/source-extractor', '--no-plots', frame_name, '--ra', str(ra), '--dec', str(dec), '--radius', str(radius)])
		else:
			subprocess.run(['solve-field', '--use-source-extractor', '--source-extractor-path', '/usr/bin/source-extractor', '--no-plots', frame_name, '--ra', str(ra), '--dec', str(dec), '--radius', str(radius)], stdout=subprocess.DEVNULL)

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
	def source_aperture_photometry(frame_data, frame_header, ra, dec):

		# load directories
		dirs = Utility.get_directories(location='output', create=False)

		# load mask
		mask, boxes = Photometry.mk_mask(frame_data, frame_header)

		# gather exposure time information
		exp_time = frame_header['EXPTIME']

		# gather coordinate information
		wcs = WCS(frame_header)

		# gather time information
		time = Time(frame_header['DATE'], format='isot', scale='utc')
		jd = time.jd

		# calculate global background statistics
		bkg_mn, bkg_md, bkg_sd = sigma_clipped_stats(frame_data, mask=mask, sigma=Configuration.SIG_BKG)

		# calculate airmass
		location = EarthLocation(lat=Observatory.latitude(), lon=Observatory.longitude(), height=Observatory.elevation())
		altitude = Calculator.altitude(location, time, ra, dec)
		airmass = Calculator.airmass(altitude)

		sky_pos = SkyCoord(ra, dec, unit='deg')
		pix_pos = skycoord_to_pixel(sky_pos, wcs=wcs)
		pix_coord = (pix_pos[0], pix_pos[1])

		aperture = CircularAperture(pix_coord, r=Configuration.APER_SIZE)
		aperture_area = aperture.area
		annulus = CircularAnnulus(pix_coord, r_in=Configuration.ANNULI_INNER, r_out=Configuration.ANNULI_OUTER)

		aperstats = ApertureStats(frame_data, annulus)
		bkg_mean = aperstats.mean
		bkg_sum = np.abs(frame_header['SKY'] + bkg_mean) * aperture_area

		source_table = aperture_photometry(frame_data, aperture, mask=mask, wcs=wcs)

		tot_sum = source_table['aperture_sum'][0]
		res_sum = (tot_sum - bkg_mean * aperture_area) * Observatory.gain()
		
		flux = res_sum
		if flux < 0:
			flux = abs(bkg_md)

		flux_err = np.sqrt(flux + aperture_area * bkg_sd ** 2)
		inst_mag = 25. -2.5*np.log10(flux) + 2.5*np.log10(exp_time)
		inst_mag_err = (2.5 * flux_err) / (np.log(10.) * flux)

		print('RA:', ra, 'deg')
		print('Dec:', dec, 'deg')
		print('x:', pix_pos[0])
		print('y:', pix_pos[1])
		print('Bkg Mean:', bkg_mean)
		print('Bkg Sum:', bkg_sum)
		print('Total Flux:', tot_sum)
		print('Residual Flux:', res_sum)
		print('Flux Error:', flux_err)
		print('Inst Mag:', inst_mag)
		print('Inst Mag Error:', inst_mag_err)

		return flux, flux_err, inst_mag, inst_mag_err, airmass, jd

	@staticmethod
	def subtract_background(raw_frame_data, raw_frame_header, ini_bkg_img_path, res_bkg_img_path):
		'''This function measures and subtracts the background signal from an image data array.

		:parameter raw_frame_data - 
		:parameter raw_frame_header - 
		:parameter ini_bkg_img_path - 
		:parameter res_bkg_img_path - 

		:return clean_frame_data - 
		:return cf_header - 
		:return bkg2d_ini_md - 
		:return bkg2d_ini_sd - 
		:return bkgsc_ini_mn - 
		:return bkgsc_ini_md - 
		:return bkgsc_ini_sd - 
		:return bkg2d_res_md - 
		:return bkg2d_res_sd - 
		:return bkgsc_res_mn - 
		:return bkgsc_res_md - 
		:return bkgsc_res_sd - 
		'''

		# measure the initial background
		bkg2d_ini_data, bkg2d_ini_md, bkg2d_ini_sd, bkgsc_ini_mn, bkgsc_ini_md, bkgsc_ini_sd = Photometry.measure_background(raw_frame_data)

		# save the initial background
		if not os.path.exists(ini_bkg_img_path):
			Plot.field(bkg2d_ini_data, ini_bkg_img_path)

		# subtract the initial background
		bkg_method = Configuration.BKG_METHOD
		if bkg_method == 'flat':
			clean_frame_data = raw_frame_data - bkg2d_ini_md
		elif bkg_method == '2d':
			clean_frame_data = raw_frame_data - bkg2d_ini_data

		# measure the residual background
		bkg2d_res_data, bkg2d_res_md, bkg2d_res_sd, bkgsc_res_mn, bkgsc_res_md, bkgsc_res_sd = Photometry.measure_background(clean_frame_data)

		# save the residual background
		if not os.path.exists(res_bkg_img_path):
			Plot.field(bkg2d_res_data, res_bkg_img_path)
			
		# edit the header
		cf_header = raw_frame_header
		cf_header['BKGSB'] = ('median-' + bkg_method, 'Background method')
		cf_header['BKGTH'] = (Configuration.SIG_BKG, 'Background threshold')
		cf_header['BKGMD'] = (bkg2d_ini_md, 'Median [ADU]')
		cf_header['BKGSD'] = (bkg2d_ini_sd, 'Median RMS [ADU]')
		cf_header['RESMD'] = (bkg2d_res_md, 'Residual median [ADU]')
		cf_header['RESSD'] = (bkg2d_res_sd, 'Residual RMS [ADU]')

		return clean_frame_data, cf_header, bkg2d_ini_md, bkg2d_ini_sd, bkgsc_ini_mn, bkgsc_ini_md, bkgsc_ini_sd, bkg2d_res_md, bkg2d_res_sd, bkgsc_res_mn, bkgsc_res_md, bkgsc_res_sd

	@staticmethod
	def subtract_dark(raw_frame_data, raw_frame_header, master_dark_frame_path):

		master_dark_frame = fits.open(master_dark_frame_path)
		master_dark_data = master_dark_frame[0].data
		master_dark_header = master_dark_frame[0].header

		clean_frame_data = raw_frame_data - master_dark_data
		clean_frame_header = raw_frame_header
		clean_frame_header['DRKSB'] = ('yes', 'Master dark subtracted')

		return clean_frame_data, clean_frame_header