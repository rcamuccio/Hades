from config import Configuration
from lib.utilities import Utils
from astroplan import Observer
from astropy import units as u
from astropy.coordinates import EarthLocation, get_sun, SkyCoord
from astropy.time import Time
from datetime import timedelta
from ligo.skymap.postprocess import crossmatch
import numpy as np
import os
import pandas as pd
import time

class Priority:

	@staticmethod
	def field_generator(field_sep):
		'''This function will generate survey fields based on the provided separation of each field.

		:parameter field_sep - The desired field separation required between each field center [deg]

		:return field_list - The field list is both saved to file and returned as a Pandas DataFrame
		'''

		# start clock
		st = time.time()

		# pull in the observatory
		observatory = Utils.config_observatory(Configuration.OBSERVATORY)
		output_dirs = Utils.config_output_dir(create=False)

		# set the survey field catalog path
		field_path = output_dirs['analysis'] + Configuration.OBSERVATORY + '_fields.cat'

		# initialize parameters
		deg_to_rad = np.pi / 180.
		gal_idx = 0
		tot_idx = 1

		if not os.path.isfile(field_path):
			Utils.log('No existing survey fields for ' + str(Configuration.OBSERVATORY).upper(), 'info')

			# set up the field data frame and start with the first field
			if observatory['latitude'] < 0e0:
				# create sky coordinate object of first field
				c = SkyCoord(ra=0e0, dec=-90e0, frame='icrs', unit='deg')

				# initialize the data frame with the first field
				field_list = pd.DataFrame(data=[['00.000', 0e0, -90e0, c.galactic.l.to_value(), c.galactic.b.to_value(), 'main_survey', observatory['exp_light'], 1, 0, 0, 0, 0]], columns=['field_id', 'ra', 'dec', 'l', 'b', 'program', 'exposure_time', 'cadence', 'ephemeris', 'period', 'observations', 'moon_phase'])

				# calculate number of declination strips
				field_number = int(np.ceil((90 + observatory['dec_limit']) / field_sep))

				# create array of declination strips
				declination_strips = -90 + np.arange(0, field_number) * field_sep

			else:
				# create sky coordinate object of first field
				c = SkyCoord(ra=0e0, dec=90e0, frame='icrs', unit='deg')
				
				# initialize the data frame with the first field
				field_list = pd.DataFrame(data=[['00.000', 0e0, 90e0, c.galactic.l.to_value(), c.galactic.b.to_value(), 'main_survey', observatory['exp_light'], 1, 0, 0, 0, 0]], columns=['field_id', 'ra', 'dec', 'l', 'b', 'program', 'exposure_time', 'cadence', 'ephemeris', 'period', 'observations', 'moon_phase'])

				# calculate number of declination strips
				field_number = int(np.ceil((90 - observatory['dec_limit']) / field_sep))

				# create array of declination strips
				declination_strips = 90 - np.arange(0, field_number) * field_sep

			Utils.log('Number of declination strips: ' + str(field_number), 'info')
			Utils.log('Generating survey fields for ' + str(Configuration.OBSERVATORY).upper(), 'info')

			# loop through and generate the fields
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
					# set up the field name with the first hex field
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

					field_name = field_1 + '.' + field_2
					field_ra = idy * ra_sep + ra_off
					field_de = declination_strips[idx]

					# get the galactic coordinates of the field
					c_idy = SkyCoord(ra=field_ra, dec=field_de, frame='icrs', unit='deg')
					if (c_idy.galactic.b.to_value() < Configuration.GALACTIC_PLANE) & (c_idy.galactic.b.to_value() > -1 * Configuration.GALACTIC_PLANE):
						moon_phase = 1
						gal_idx += 1
					else:
						moon_phase = 0

					# set up the series for appending
					field = pd.Series(data=[field_name, field_ra, field_de, c_idy.galactic.l.to_value(), c_idy.galactic.b.to_value(), 'main_survey', observatory['exp_light'], 1, 0, 0, 0, moon_phase], index=['field_id', 'ra', 'dec', 'l', 'b', 'program', 'exposure_time', 'cadence', 'ephemeris', 'period', 'observations', 'moon_phase'])

					# append the series
					field_list.loc[tot_idx] = field
					tot_idx += 1

			field_list.to_csv(field_path, sep=' ', header=True, index=False, float_format='%.3f')

			Utils.log('Field generation complete!', 'info')
			Utils.log('The main survey for ' + str(Configuration.OBSERVATORY).upper() + ' consists of:', 'info')
			Utils.log('    Total survey fields: ' + str(tot_idx), 'info')
			Utils.log('    Galactic fields: ' + str(gal_idx), 'info')
			Utils.log('    Extragalactic fields: ' + str(tot_idx - gal_idx), 'info')

		else:
			# if the file exists already, then just read the field list in
			Utils.log('Reading extant survey fields.', 'info')
			field_list = pd.read_csv(field_path, header=0, sep=' ')

		fn = time.time()
		Utils.log('Survey fields generated in ' + str(np.around((fn - st), decimals=2)) + ' s.', 'info')

		return field_list

	@staticmethod
	def find_field(survey_fields, ra, de, moon_phase=2, ang_extent=0, program='main_survey'):
		'''This function will select the appropriate survey field based on the RA/DE of the target you are interested in observing.

		:parameter survey_fields - The set of survey fields
		:parameter ra - The right ascension of the target [deg]
		:parameter de - The declination of the target [deg]
		:parameter moon_phase - Update this to be 0 for dark time, 1 for bright time, or 2 for any time
		:parameter ang_extent - If this is an extended object, provide its angular radius to grab moer than one survey field
		:parameter program - If the input is something other than 'main_survey', the program is overwritten

		:return fields - 
		:return ang_dists - 
		'''

		# load the observatory
		observatory = Utils.config_observatory(Configuration.OBSERVATORY)

		# get the angular distance between the given position and the survey fields
		ang_dist = survey_fields.apply(lambda x: Utils.angular_distance(float(x.ra), float(x.de), ra, de), axis=1)

		# if angular extent is given, find the fields within that extent
		if ang_extent > 0:
			if len(survey_fields[ang_dist < ang_extent]) > 0:
				fields = survey_fields[ang_dist < ang_extent].copy()
				ang_dists = ang_dist[ang_dist < ang_extent]

				if program != 'main_survey':
					fields['program'] = program
					fields['moon_phase'] = moon_phase

			else:
				fields = 'No fields found.'

		# if angular extent is not given, find the closest field within a field radius
		else:
			if np.min(ang_dist) <= observatory['fov'] / 2:
				fields = survey_fields.loc[np.argmin(ang_dist)].copy()
				ang_dists = np.argmin(ang_dist)

				if program != 'main_survey':
					fields['program'] = program
					fields['moon_phase'] = moon_phase

			else:
				fields = 'No field found.'

		return fields, ang_dists

	@staticmethod
	def generate_targets():

		observatory = Utils.config_observatory(Configuration.OBSERVATORY)
		output_dirs = Utils.config_output_dir(create=False)

		if Configuration.CATALOG == 'glade24':
			path =	output_dirs['catalogs'] + 'GLADE_2.4.txt'
			usecols = [1, 6, 7, 8]
			Utils.log('Reading GLADE 2.4 galaxy catalog.', 'info')

		elif Configuration.CATALOG == 'glade+':
			path =	output_dirs['catalogs'] + 'GLADE+.txt'
			usecols = [2, 8, 9, 32]
			Utils.log('Reading GLADE+ galaxy catalog.', 'info')

		delimiter = ' '
		header = None
		names = ['GWGC', 'ra', 'dec', 'dist']
		low_memory = False

		df = pd.read_csv(catalog, delimiter=delimiter, usecols=usecols, header=header, names=names, low_memory=low_memory)

		Utils.log('Slicing catalog.', 'info')

		# drop NaNs from catalog
		df = df.dropna()

		# make declination cuts
		lim_dec = df.dec > observatory['dec_limit']
		new_df = df[lim_dec]

		# make right ascension cuts
		right_now = Time.now()
		right_now = Time(right_now)

		sun = get_sun(right_now)
		horizon = -15 * u.degree
		min_height = 30 * u.degree

		alpha_obs_min = (sun.ra - horizon + min_height).degree
		alpha_obs_max = (sun.ra + horizon - min_height).degree

		circum = 90. - abs(new_df.dec) < abs(observatory['latitude'])
		alfa_min = new_df.ra > float(alpha_obs_min)
		alfa_max = new_df.ra <= float(alpha_obs_max)

		case_1 = (alfa_min & alfa_max) | circum
		case_2 = (alfa_min | alfa_max) | circum

		if alpha_obs_max > alpha_obs_min:
			final_df = new_df[case_1]
		else:
			final_df = new_df[case_2]

		Utils.log('Catalog slicing complete.', 'info')

		return final_df

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