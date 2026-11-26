'''
Σκοπεύς - The Observer
'''

import numpy as np
import os
import pandas as pd
import time
from astroplan import Observer
from astropy import units as u
from astropy.coordinates import EarthLocation, get_sun, SkyCoord
from astropy.time import Time
from datetime import timedelta
from ligo.skymap.postprocess import crossmatch

from conf import Configuration
from lib.logistes import Logistes

class Skopeus:

	@staticmethod
	def find_survey_field(survey_fields, ra, dec, moon_phase=2, ang_ext=0, program='main_survey'):
		'''This function selects the appropriate survey field (or set of fields) based on the right ascension, declination, and angular extent of the target of interest.

		:parameter survey_fields [] - the set of survey fields
		:parameter ra [float] - the right ascension of the target [deg]
		:parameter dec [float] - the declination of the target [deg]
		:parameter moon_phase [int] - the moon phase (0 = dark time, 1 = bright time, 2 = any time)
		:parameter ang_ext [float] - the angular radius of the object
		:parameter program [str] - the survey mode

		:return fields [] - the field (or set of fields)
		:return ang_dists [] - 

		return fields, ang_dists
		'''

		# get the angular distance between the survey fields and the target
		ang_dist = survey_fields.apply(lambda x: Logistes.get_angular_distance(float(x.ra), float(x.dec), ra, dec), axis=1)

		# if the target has angular extent find the fields within that extent
		if ang_ext > 0:
			if len(survey_fields[ang_dist < ang_ext]) > 0:
				fields = survey_fields[ang_dist < ang_ext].copy()
				ang_dists = ang_dist[ang_dist < ang_ext]

				if program != 'main_survey':
					fields['program'] = program
					fields['moon_phase'] = moon_phase

			else:
				fields = 'No fields found.'

		# otherwise find the closest field within a field radius
		else:
			if np.min(ang_dist) <= Skopeus.get_fov() / 2:
				fields = survey_fields.loc[np.argmin(ang_dist)].copy()
				ang_dists = np.argmin(ang_dist)

				if program != 'main_survey':
					fields['program'] = program
					fields['moon_phase'] = moon_phase

			else:
				fields = 'No field found.'

		return fields, ang_dists

	@staticmethod
	def generate_survey_fields(field_sep):
		'''This function generates a set of survey fields based on a given field separation.

		:p field_sep [float] - the separation between each field center [deg]

		:r field_list [Pandas.DataFrame] - the set of survey fields
		'''

		# set the survey field path
		field_path = Configuration.PLOUTON_DIRECTORY + 'fields/' + Configuration.OBSERVATORY + '_fields' + Configuration.TABLE_EXTENSION

		# initialize the parameters
		deg_to_rad = np.pi / 180.
		gal_idx = 0
		tot_idx = 1

		if not os.path.isfile(field_path):
			print('Generating survey fields for', Configuration.OBSERVATORY.upper())
			st = time.time()

			if Skopeus.get_latitude() < 0e0:
				# create a sky coordinate object for the first field
				c = SkyCoord(ra=0e0, dec=-90e0, frame='icrs', unit='deg')

				# initialize the data frame with the first field
				field_list = pd.DataFrame(data=[['00.000', 0e0, -90e0, c.galactic.l.to_value(), c.galactic.b.to_value(), 'main_survey', Skopeus.get_exposure_time('light'), 1, 0, 0, 0, 0]], columns=['field_id', 'ra', 'dec', 'l', 'b', 'program', 'exposure_time', 'cadence', 'ephemeris', 'period', 'observations', 'moon_phase'])

				# calculate the number of declination strips
				field_number = int(np.ceil((90 + Skopeus.get_declination_limit()) / field_sep))

				# create the array of declination strips
				declination_strips = -90 + np.arange(0, field_number) * field_sep

			else:
				# create a sky coordinate object for the first field
				c = SkyCoord(ra=0e0, dec=90e0, frame='icrs', unit='deg')

				# initialize the data frame with the first field
				field_list = pd.DataFrame(data=[['00.000', 0e0, 90e0, c.galactic.l.to_value(), c.galactic.b.to_value(), 'main_survey', Skopeus.get_exposure_time('light'), 1, 0, 0, 0, 0]], columns=['field_id', 'ra', 'dec', 'l', 'b', 'program', 'exposure_time', 'cadence', 'ephemeris', 'period', 'observations', 'moon_phase'])

				# calculate the number of declination strips
				field_number = int(np.ceil((90 - Skopeus.get_declination_limit()) / field_sep))

				# create the array of declination strips
				declination_strips = 90 - np.arange(0, field_number) * field_sep

			print('Number of declination strips:', field_number)

			# generate the fields
			eo = 0
			for idx in range(1, field_number):
				nfr = np.ceil(360. * np.cos(declination_strips[idx] * deg_to_rad) / field_sep)
				ra_sep = 360. / nfr

				if eo == 1:
					ra_off = 0.5 * ra_sep
					eo = 0
				else:
					ra_off = 0.
					eo = 1

				for idy in range(0, int(nfr)):
					# set up the first hex field
					if len(hex(idx).split('x')[1]) == 1:
						field_1 = '0' + hex(idx).split('x')[1]
					else:
						field_1 = hex(idx).split('x')[1]

					# set up the second hex field
					if len(hex(idy).split('x')[1]) == 1:
						field_2 = '00' + hex(idy).split('x')[1]
					elif len(hex(idy).split('x')[1]) == 2:
						field_2 = '0' + hex(idy).split('x')[1]
					else:
						field_2 = hex(idy).split('x')[1]

					# set up the field name and position
					field_name = field_1 + '.' + field_2
					field_ra = idy * ra_sep + ra_off
					field_dec = declination_strips[idx]

					# get the galactic coordinates of the field
					c_idy = SkyCoord(ra=field_ra, dec=field_dec, frame='icrs', unit='deg')
					if (c_idy.galactic.b.to_value() < Configuration.GALACTIC_PLANE) & (c_idy.galactic.b.to_value() > -1 * Configuration.GALACTIC_PLANE):
						moon_phase = 1
						gal_idx += 1
					else:
						moon_phase = 0

					# set up the series for appending
					field = pd.Series(data=[field_name, field_ra, field_dec, c_idy.galactic.l.to_value(), c_idy.galactic.b.to_value(), 'main_survey', Skopeus.get_exposure_time('light'), 1, 0, 0, 0, moon_phase], index=['field_id', 'ra', 'dec', 'l', 'b', 'program', 'exposure_time', 'cadence', 'ephemeris', 'period', 'observations', 'moon_phase'])

					# append the series
					field_list.loc[tot_idx] = field
					tot_idx += 1

			# write the field list to file
			field_list.to_csv(field_path, sep=' ', header=True, index=False, float_format='%.3f')

			print('Field generation complete')
			print('The main survey for', Configuration.OBSERVATORY.upper(), 'consists of:')
			print('\tTotal survey fields:', tot_idx)
			print('\tGalactic fields:', gal_idx)
			print('\tExtragalactic fields:', tot_idx - gal_idx)

			fn = time.time()
			print('Survey fields generated in', np.around((fn - st), decimals=2), 's')

		else:
			print('Reading existing survey fields for', Configuration.OBSERVATORY.upper())
			field_list = pd.read_csv(field_path, header=0, sep=' ')

		return field_list

	@staticmethod
	def generate_targets():

		if Configuration.GALAXY_CATALOG == 'glade24':
			print('Reading GLADE 2.4 galaxy catalog')
			path = Configuration.PLOUTON_DIRECTORY + 'catalogs/GLADE_2.4.txt'
			usecols = [1, 6, 7, 8]

		elif Configuration.GALAXY_CATALOG == 'glade+':
			print('Reading GLADE+ galaxy catalog')
			path = Configuration.PLOUTON_DIRECTORY + 'catalogs/GLADE+.txt'
			usecols = [2, 8, 9, 32]

		delimiter = ' '
		header = None
		names = ['GWGC', 'ra', 'dec', 'dist']
		low_memory = False

		df = pd.read_csv(catalog, delimiter=delimiter, usecols=usecols, header=header, names=names, low_memory=low_memory)

		print('Slicing catalog')

		# drop NaNs from the catalog
		df = df.dropna()

		# make declination cuts to the catalog
		lim_dec = df.dec > Skopeus.get_declination_limit()
		new_df = df[lim_dec]

		# make right ascension cuts to the catalog
		right_now = Time.now()
		right_now = Time(right_now)

		sun = get_sun(right_now)
		horizon = -15 * u.degree
		min_height = 30 * u.degree

		alpha_obs_min = (sun.ra - horizon + min_height).degree
		alpha_obs_max = (sun.ra + horizon - min_height).degree

		circum = 90. - abs(new_df.dec) < abs(Skopeus.get_latitude())
		alfa_min = new_df.ra > float(alpha_obs_min)
		alfa_max = new_df.ra <= float(alpha_obs_max)

		case_1 = (alfa_min & alfa_max) | circum
		case_2 = (alfa_min | alfa_max) | circum

		if alpha_obs_max > alpha_obs_min:
			final_df = new_df[case_1]
		else:
			final_df = new_df[case_2]

		print('Galaxy targets generated')

		return final_df

	@staticmethod
	def get_axis_x():

		if Configuration.OBSERVATORY == 'toros':
			axis_x = 10560
		else:
			axis_x = None

		return axis_x

	@staticmethod
	def get_axis_x_raw():

		if Configuration.OBSERVATORY == 'toros':
			axis_x_raw = 12000
		else:
			axis_x_raw = None

		return axis_x_raw

	@staticmethod
	def get_axis_y():

		if Configuration.OBSERVATORY == 'toros':
			axis_y = 10560
		else:
			axis_y = None

		return axis_y

	@staticmethod
	def get_axis_y_raw():

		if Configuration.OBSERVATORY == 'toros':
			axis_y_raw = 10600
		else:
			axis_y_raw = None

		return axis_y_raw

	@staticmethod
	def get_bandpass():

		if Configuration.OBSERVATORY == 'toros':
			bandpass = [147, 141, 147, 147]
		else:
			bandpass = None

		return bandpass

	@staticmethod
	def get_binning():

		if Configuration.OBSERVATORY == 'toros':
			binning = 1
		else:
			binning = None

		return binning

	@staticmethod
	def get_center_wavelength():

		if Configuration.OBSERVATORY == 'toros':
			center_wavelength = [473.5, 638.5, 775.5, 922.5]
		else:
			center_wavelength = None

		return center_wavelength

	@staticmethod
	def get_declination_limit():

		if Configuration.OBSERVATORY == 'toros':
			declination_limit = 26.66
		else:
			declination_limit = None

		return declination_limit

	@staticmethod
	def get_elevation():

		if Configuration.OBSERVATORY == 'toros':
			elevation = 2420
		else:
			elevation = None

		return elevation

	@staticmethod
	def get_exposure_number():

		if Configuration.OBSERVATORY == 'toros':
			exposure_number = 1
		else:
			exposure_number = None

		return exposure_number

	@staticmethod
	def get_exposure_time(type='light', scale='s'):

		if Configuration.OBSERVATORY == 'toros':
			if type == 'dark-f':
				exposure_time = 5.
			elif type == 'dark-l':
				exposure_time = 300.
			elif type == 'flat':
				exposure_time = 5.
			elif type == 'light':
				exposure_time = 300.

			if scale == 'd':
				exposure_time = exposure_time / 60. / 60. / 24.

		return exposure_time

	@staticmethod
	def get_exposure_total():

		if Configuration.OBSERVATORY == 'toros':
			exposure_number = Skopeus.get_exposure_number()
			exposure_time = Skopeus.get_exposure_time()
			overhead = Skopeus.get_overhead()
			read_time = Skopeus.get_read_time()

			exposure_total = (exposure_time + read_time) * exposure_number + overhead

		else:
			exposure_total = None

		return exposure_total

	@staticmethod
	def get_field_dates(field_id):
		'''This function retrieves a list of observation dates for a given survey field.

		:parameter field_id [str] - the survey field hexadecimal code
		
		:return date_list [list] - the list of dates
		'''

		# read the field list
		field_path = Configuration.PLOUTON_DIRECTORY + 'fields/comm_fields' + Configuration.TABLE_EXTENSION

		field_file = open(field_path, 'r')

		# find each of the field dates
		for ln in field_file:
			ln = ln.split()
			if ln[0] == field_id:
				date_list = ln[-1].split('/')

		return date_list

	@staticmethod
	def get_commissioning_fields():

		field_path = Configuration.PLOUTON_DIRECTORY + 'field/comm_fields' + Configuration.TABLE_EXTENSION

		field_file = open(field_path, 'r')

		com_field_ra_list = []
		com_field_dec_list = []
		grb_field_ra_list = []
		grb_field_dec_list = []
		gw_field_ra_list = []
		gw_field_dec_list = []
		ot_field_ra_list = []
		ot_field_dec_list = []
		sn_field_ra_list = []
		sn_field_dec_list = []
		xrb_field_ra_list = []
		xrb_field_dec_list = []

		for ln in field_file:
			ln = ln.split()

			field_id = ln[0]

			field_ra = float(ln[1])
			field_ra = np.remainder(field_ra + 180, 360) - 180
			field_ra = float(np.radians(field_ra))

			field_dec = float(ln[2])
			field_dec = float(np.radians(field_dec))

			field_type = ln[3]

			event_type = ln[4]

			if field_type == 'C':
				com_field_ra_list.append(field_ra)
				com_field_dec_list.append(field_dec)

			else:
				if event_type == 'GW':
					gw_field_ra_list.append(field_ra)
					gw_field_dec_list.append(field_dec)

				elif (event_type == 'LGRB') or (event_type == 'SGRB'):
					grb_field_ra_list.append(field_ra)
					grb_field_dec_list.append(field_dec)

				elif event_type == 'OT':
					ot_field_ra_list.append(field_ra)
					ot_field_dec_list.append(field_dec)

				elif (event_type == 'SNIA') or (event_type == 'SNII'):
					sn_field_ra_list.append(field_ra)
					sn_field_dec_list.append(field_dec)

				elif event_type == 'XRB':
					xrb_field_ra_list.append(field_ra)
					xrb_field_dec_list.append(field_dec)

		return com_field_ra_list, com_field_dec_list, gw_field_ra_list, gw_field_dec_list, grb_field_ra_list, grb_field_dec_list, ot_field_ra_list, ot_field_dec_list, sn_field_ra_list, sn_field_dec_list, xrb_field_ra_list, xrb_field_dec_list

	@staticmethod
	def get_survey_fields():

		field_path = Configuration.PLOUTON_DIRECTORY + 'field/toros_fields' + Configuration.TABLE_EXTENSION

		field_file = open(field_path, 'r')

		gal_ra_list = []
		gal_dec_list = []
		egal_ra_list = []
		egal_dec_list = []

		for ln in field_file:
			ln = ln.split()
			if ln[0] != 'field_id':
				ra = float(ln[1])
				ra = np.remainder(ra + 180, 360) - 180
				ra = np.radians(ra)

				dec = float(ln[2])
				dec = np.radians(dec)

				l = float(ln[3])
				b = float(ln[4])

				if abs(b) <= Configuration.GALACTIC_PLANE:
					gal_ra_list.append(ra)
					gal_dec_list.append(dec)
				else:
					egal_ra_list.append(ra)
					egal_dec_list.append(dec)

		return gal_ra_list, gal_dec_list, egal_ra_list, egal_dec_list

	@staticmethod
	def get_field_position(field_id):
		'''This function returns the right ascension and declination for a given survey field.

		:parameter field_id [str] - the survey field hexadecimal code

		:return ra [float] - the right ascension of the field [deg]
		:return dec [float] - the declination of the field [deg]
		'''

		# read the field list
		field_path = Configuration.PLOUTON_DIRECTORY + 'fields/toros_fields' + Configuration.TABLE_EXTENSION
		field_list = pd.read_csv(field_path, header=0, sep=' ')

		# find the field position
		for i in range(len(field_list.field_id)):
			if field_id == field_list.field_id[i]:
				ra = float(field_list.iloc[i].ra)
				dec = float(field_list.iloc[i].dec)

		return ra, dec

	@staticmethod
	def get_fov():

		if Configuration.OBSERVATORY == 'toros':
			pixel_scale = Skopeus.get_pixel_scale()
			axis_x = Skopeus.get_axis_x()
			fov = (pixel_scale * axis_x) / 3600.
		else:
			fov = None

		return fov

	@staticmethod
	def get_gain():

		if Configuration.OBSERVATORY == 'toros':
			gain = 0.380
		else:
			gain = None

		return gain

	@staticmethod
	def get_latitude():

		if Configuration.OBSERVATORY == 'toros':
			latitude = -31.8023
		else:
			latitude = None

		return latitude

	@staticmethod
	def get_longitude():

		if Configuration.OBSERVATORY == 'toros':
			longitude = -69.3265
		else:
			longitude = None

		return longitude

	@staticmethod
	def get_mirror_diameter():

		if Configuration.OBSERVATORY == 'toros':
			mirror_diameter = 0.610
		else:
			mirror_diameter = None

		return mirror_diameter

	@staticmethod
	def get_mirror_radius():

		if Configuration.OBSERVATORY == 'toros':
			mirror_diameter = Skopeus.get_mirror_diameter()
			mirror_radius = mirror_diameter / 2
		else:
			mirror_radius = None

		return mirror_radius

	@staticmethod
	def get_overhead():

		if Configuration.OBSERVATORY == 'toros':
			overhead = 30.
		else:
			overhead = None

		return overhead

	@staticmethod
	def get_overscan_x():

		if Configuration.OBSERVATORY == 'toros':
			overscan_x = 180
		else:
			overscan_x = None

		return overscan_x

	@staticmethod
	def get_overscan_y():

		if Configuration.OBSERVATORY == 'toros':
			overscan_y = 20
		else:
			overscan_y = None

		return overscan_y

	@staticmethod
	def get_peakmax():

		if Configuration.OBSERVATORY == 'toros':
			peakmax = 45000.
		else:
			peakmax = None

		return peakmax

	@staticmethod
	def get_pixel_scale():

		if Configuration.OBSERVATORY == 'toros':
			pixel_scale = 0.4959
		else:
			pixel_scale = None

		return pixel_scale

	@staticmethod
	def get_qe_atmosphere():

		if Configuration.OBSERVATORY == 'toros':
			qe_atmosphere = [0.8, 0.9, 0.9, 0.9]
		else:
			qe_atmosphere = None

		return qe_atmosphere

	@staticmethod
	def get_qe_camera():

		if Configuration.OBSERVATORY == 'toros':
			qe_camera = 0.85
		else:
			qe_camera = None

		return qe_camera

	@staticmethod
	def get_qe_filter():

		if Configuration.OBSERVATORY == 'toros':
			qe_filter = 0.9
		else:
			qe_filter = None

		return qe_filter

	@staticmethod
	def get_qe_primary():

		if Configuration.OBSERVATORY == 'toros':
			qe_primary = 0.96
		else:
			qe_primary = None

		return qe_primary

	@staticmethod
	def get_qe_secondary():

		if Configuration.OBSERVATORY == 'toros':
			qe_secondary = 0.96
		else:
			qe_secondary = None

		return qe_secondary

	@staticmethod
	def get_read_noise():

		if Configuration.OBSERVATORY == 'toros':
			read_noise = 5.
		else:
			read_noise = None

		return read_noise

	@staticmethod
	def get_read_time():

		if Configuration.OBSERVATORY == 'toros':
			read_time = 90.
		else:
			read_time = None

		return read_time

	@staticmethod
	def get_seeing():

		if Configuration.OBSERVATORY == 'toros':
			seeing = 0.93
		else:
			seeing = None

		return seeing

	@staticmethod
	def get_sky():

		if Configuration.OBSERVATORY == 'toros':
			sky = [22.1, 21.1, 20.1, 18.7]
		else:
			sky = None

		return sky

	@staticmethod
	def get_tap_x():

		if Configuration.OBSERVATORY == 'toros':
			tap_x = 1320
		else:
			tap_x = None

		return tap_x

	@staticmethod
	def get_tap_y():

		if Configuration.OBSERVATORY == 'toros':
			tap_y = 5280
		else:
			tap_y = None

		return tap_y

	@staticmethod
	def get_throughput():

		if Configuration.OBSERVATORY == 'toros':
			qe_atmosphere = float(np.mean(Skopeus.get_qe_atmosphere()))
			qe_camera = Skopeus.get_qe_camera()
			qe_filter = Skopeus.get_qe_filter()
			qe_primary = Skopeus.get_qe_primary()
			qe_secondary = Skopeus.get_qe_secondary()
			vignetting = Skopeus.get_vignetting()

			throughput = qe_atmosphere * qe_camera * qe_filter * qe_primary * qe_secondary * vignetting

		else:
			throughput = None

		return throughput

	@staticmethod
	def get_utc():

		if Configuration.OBSERVATORY == 'toros':
			utc = -3
		else:
			utc = None

		return utc

	@staticmethod
	def get_vignetting():

		if Configuration.OBSERVATORY == 'toros':
			vignetting = 0.756
		else:
			vignetting = None

		return vignetting

	@staticmethod
	def sort_field_skymap(survey_fields, skymap, event_name):

		# copy the data frame so we don't ruin the original data
		selected_fields = survey_fields.copy().reset_index(drop=True)

		# get the coordinate values for each survey field in appropriate units
		field_coords = SkyCoord(selected_fields.ra, selected_fields.dec, unit=u.deg)

		# crossmatch the survey fields with the sky map probabilities
		cross_match = crossmatch(skymap, field_coords)

		# move through each survey field and get the integrated area within field area
		selected_fields['prob'] = cross_match.probdensity

		# sort all survey fields based on the probability strip
		selected_fields = selected_fields.sort_values(by='prob', ascending=False).copy().reset_index(drop=True)
		selected_fields['program'] = 'lvc_alert_' + event_name

		return selected_fields

	@staticmethod
	def sort_galaxy_skymap(catalog, skymap, event_name):

		selected_galaxies = catalog.copy().reset_index(drop=True)

		galaxy_coords = SkyCoord(selected_galaxies.ra, selected_galaxies.dec, unit=u.deg)

		cross_match = crossmatch(skymap, galaxy_coords)

		selected_galaxies['prob'] = cross_match.probdensity

		selected_galaxies = selected_galaxies.sort_values(by='prob', ascending=False).copy().reset_index(drop=True)

		selected_galaxies = selected_galaxies.head(Configuration.GALAXY_LIST)

		return selected_galaxies