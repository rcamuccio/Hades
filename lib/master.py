from config import Configuration
from lib.photometry import Photometry
from lib.preprocessing import Preproc
from lib.utilities import Utils
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astroquery.mast import Catalogs
from photutils.aperture import aperture_photometry, ApertureStats, CircularAnnulus, CircularAperture
from photutils.centroids import centroid_sources
import numpy as np
import os
import pandas as pd

class Master:

	@staticmethod
	def master_phot(master, master_header):
		'''This script will generate the star list for the master frame.

		:parameter master - 
		:parameter master_header - 

		:return star_list - 
		'''

		observatory = Utils.config_observatory(Configuration.OBSERVATORY)
		data_dirs = Utils.config_data_dir(create=False)

		if not os.path.isfile(data_dirs['master'] + Configuration.FIELD + '_star_list.csv'):

			if not os.path.isfile(data_dirs['master'] + Configuration.FIELD + '_gaia_dump.csv'):

				# create the string useful for the query region
				field_str = Configuration.FIELD
				field_hex = field_str.split('_')[1]
				field_ra, field_de = Utils.config_field(field_hex)
				field = str(field_ra) + ' ' + str(field_de)

				# select the columns we want to import into the data table
				columns = ['field_id', 'source_id', 'ra', 'dec', 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'teff_val', 'parallax', 'parallax_error', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error']

				# run the query
				Utils.log('Querying MAST for ' + Configuration.FIELD, 'info')
				Utils.log('    RA [deg]: ' + str(field_ra), 'info')
				Utils.log('    DE [deg]: ' + str(field_de), 'info')
				Utils.log('    Search radius [deg]: ' + str(observatory['search_dist']), 'info')
				catalog_data = Catalogs.query_region(field, radius=observatory['search_dist'], catalog='Gaia').to_pandas()
				Utils.log('MAST query complete.', 'info')
				Utils.log('Found ' + str(len(catalog_data)) + ' stars.', 'info')

				# add the field to the catalog data
				catalog_data['field_id'] = Configuration.FIELD
				catalog_data['object_type'] = 'star'

				# pull out the necessary columns
				star_list = catalog_data[columns]

				# get the header file and convert to x/y pixel positions
				w = WCS(master_header)
				ra = star_list.ra.to_numpy()
				de = star_list.dec.to_numpy()

				# convert to x, y
				x, y = w.all_world2pix(ra, de, 0)

				# add the x/y to the star data frame
				star_list['x'] = x
				star_list['y'] = y

				# dump list to file
				star_list.to_csv(data_dirs['master'] + Configuration.FIELD + '_gaia_dump.csv', sep=' ')

			else:
				Utils.log('Reading saved query file.', 'info')
				star_list = pd.read_csv(data_dirs['master'] + Configuration.FIELD + '_gaia_dump.csv', sep=' ', index_col=0)

			# check for any known transient and variable star files
			#if Configuration.KNOWN_VARIABLES == 'Y':
			#	star_list = Photometry.add_variable_list(star_list, master_header)

			# remove stars outside the frame
			star_list = star_list[(star_list.x >= 530) & (star_list.x < 10465) & (star_list.y >= 490) & (star_list.y < 10045)].copy().reset_index(drop=True)

			# centroid the star list
			star_list['xcen'], star_list['ycen'] = centroid_sources(master, star_list.x.to_numpy(), star_list.y.to_numpy(), box_size=5)

			# identify bad indices
			bd_idxs = np.where(np.isnan(star_list.xcen) | np.isnan(star_list.ycen))
			if len(bd_idxs[0]) > 0:
				for bd_idx in bd_idxs[0]:
					star_list.loc[bd_idx, 'xcen'] = star_list.loc[bd_idx, 'x']
					star_list.loc[bd_idx, 'ycen'] = star_list.loc[bd_idx, 'y']

			# add 'tap' to the star list
			star_list['tap'] = 1
			kk = 1
			for idx in range(0, observatory['axs_x'], observatory['tap_x']):
				for idy in range(0, observatory['axs_y'], observatory['tap_y']):
					star_list['tap'] = np.where((star_list.xcen > idx) & (star_list.xcen < idx + observatory['tap_x']) & (star_list.ycen > idy) & (star_list.ycen < idy + observatory['tap_y']), kk, star_list.tap)
					kk += 1

			# centroid the positions (x, y)
			positions = star_list[['xcen', 'ycen']].copy().reset_index(drop=True)

			# set up stellar apertures and annuli
			aperture = CircularAperture(positions, r=Configuration.APER_SIZE)
			aperture_area = aperture.area
			annulus = CircularAnnulus(positions, r_in=Configuration.ANNULI_INNER, r_out=Configuration.ANNULI_OUTER)

			# get background statistics
			aperstats = ApertureStats(master, annulus)
			bkg_mean = aperstats.mean
			bkg_total = bkg_mean * aperture_area

			# run photometry to get the data table
			phot_table = aperture_photometry(master, aperture, method='exact')

			# extract the flux from the table
			star_flux = np.array(phot_table['aperture_sum']) * observatory['gain']

			# calculate the expected photometric error
			star_error = star_flux
			bkg_error = master_header['SKY'] * aperture_area * observatory['gain']

			# combine sky and signal errors in quadrature
			star_flux_err = np.sqrt(star_error + bkg_error)

			# convert to magnitude
			mag = 25. - 2.5 * np.log10(star_flux / observatory['exp_light'])
			mag_err = (np.log(10.) / 2.5) * (star_flux_err / star_flux)

			# initialize the light curve data frame
			star_list['master_mag'] = mag
			star_list['master_mag_err'] = mag_err
			star_list['master_flux'] = star_flux
			star_list['master_flux_err'] = star_flux_err
			star_list['master_sky'] = bkg_total

			# only keep stars with reasonable photometry in the master list
			star_list = star_list[star_list['master_flux'] > 0]

			# index is reset twice to ensure star ID matches brightness on the master frame
			star_list = star_list.sort_values(by='master_mag').reset_index(drop=True).reset_index()
			star_list = star_list.rename(columns={'index':'star_id'})

			star_list.to_csv(data_dirs['master'] + Configuration.FIELD + '_star_list.csv', sep=' ', index=False)

		else:
			star_list = pd.read_csv(data_dirs['master'] + Configuration.FIELD + '_star_list.csv', delimiter=' ', header=0)

		return star_list

	@staticmethod
	def mk_master():
		'''This function will make the master frame that will be used for differencing.

		:return master - 
		:return master_header -
		'''

		observatory = Utils.config_observatory(Configuration.OBSERVATORY)
		data_dirs = Utils.config_data_dir(create=False)

		file_name = Configuration.FIELD + '_master' + Configuration.FILE_EXTENSION

		if not os.path.isfile(data_dirs['master'] + file_name):
			chk_tmp_files = Utils.get_file_list(data_dirs['master_tmp'], Configuration.FILE_EXTENSION)

			# get the image list
			full_image_list, dates = Utils.get_all_files_per_field(data_dirs['clean'], Configuration.FIELD, 'clean', Configuration.FILE_EXTENSION)

			# determine number of loops needed to move through for each image
			full_nfiles = len(full_image_list)

			# get the nights the images were observed
			img_by_night = np.array([line.split('/')[7] for line in full_image_list]).reshape(-1)

			# loop through the image headers and get the sky background
			sky_values = np.zeros(full_nfiles)
			bd_wcs = np.zeros(full_nfiles)
			for idx, file in enumerate(full_image_list):
				# pull in the header file
				h_chk = fits.getheader(file)

				# get the sky background
				sky_values[idx] = h_chk['sky']

			# get the statistics on the images
			img_mn, img_mdn, img_std = sigma_clipped_stats(sky_values, sigma=2)

			# make a list of image to use for the master frame
			image_list = []
			for idx, file in enumerate(full_image_list):
				# remove nights with high sky background, bad wcs, or are on hand-picked bad nights
				if (sky_values[idx] <= img_mdn + 2 * img_std) & (bd_wcs[idx] == 0) & (img_by_night[idx] != '2024-10-03'):
					image_list.append(file)

			nfiles = len(image_list)

			if len(chk_tmp_files) == 0:
				Utils.log('No temporary master files found. Generating new ones.', 'info')

				# number of loops needed to move through for each image
				nbulk = 20

				# get integer and remainder for the combination
				full_bulk = nfiles // nbulk
				part_bulk = nfiles % nbulk

				if part_bulk > 0:
					hold_bulk = full_bulk + 1
				else:
					hold_bulk = full_bulk

				# here is the holder
				hold_data = np.ndarray(shape=(hold_bulk, observatory['axs_y'], observatory['axs_x']))

				# update log
				Utils.log('Generating master frame in ' + str(nbulk) + ' image bulks.', 'info')
				Utils.log('There are ' + str(nfiles) + ' images to combine.', 'info')
				Utils.log('There are ' + str(hold_bulk) + ' mini-files to median combine.', 'info')

				cnt_img = 0
				for kk in range(0, hold_bulk):

					# loop through the images in sets of nbulk
					if kk < full_bulk:
						# generate image holder
						block_hold = np.ndarray(shape=(nbulk, observatory['axs_y'], observatory['axs_x']))
						# generate max index
						mx_index = nbulk

					else:
						# generate image holder
						block_hold = np.ndarray(shape=(part_bulk, observatory['axs_y'], observatory['axs_x']))
						# generate max index
						mx_index = part_bulk

					# make the starting index
					loop_start = kk * nbulk
					idx_cnt = 0

					Utils.log('Making mini-file ' + str(kk) + '.', 'info')

					# now loop through the images
					for jj in range(loop_start, mx_index + loop_start):

						# read in the image directly into the block hold
						master_tmp, master_tmp_head = fits.getdata(image_list[jj], header=True)

						if (kk == 0) & (jj == 0):
							block_hold[idx_cnt] = master_tmp - master_tmp_head['sky']
							master_header = master_tmp_head
							del master_tmp
							Utils.log('Mini-file ' + str(kk) + ' image ' + str(jj) + ' aligned. ' + str(nfiles - cnt_img) + ' remain.', 'info')

						else:
							tmp = Preproc.align_img(master_tmp, master_tmp_head, master_header)
							Utils.log('Mini-file ' + str(kk) + ' image ' + str(jj) + ' aligned. ' + str(nfiles - cnt_img) + ' remain.', 'info')
							block_hold[idx_cnt] = tmp - master_tmp_head['sky']
							del tmp
							del master_tmp

						# increase the iteration
						cnt_img += 1
						idx_cnt += 1

					# median the data into a single file
					hold_data[kk] = np.median(block_hold, axis=0)
					del block_hold
					if kk < 10:
						fits.writeto(data_dirs['master_tmp'] + '0' + str(kk) + '_tmp_master' + Configuration.FILE_EXTENSION, hold_data[kk], master_tmp_head, overwrite=True)
					else:
						fits.writeto(data_dirs['master_tmp'] + str(kk) + '_tmp_master' + Configuration.FILE_EXTENSION, hold_data[kk], master_tmp_head, overwrite=True)

			else:
				Utils.log('Legacy files found. Creating master frame from these files. Delete if you do not want this!', 'info')
				hold_bulk = len(chk_tmp_files)

				# here is the holder
				hold_data = np.ndarray(shape=(hold_bulk, observatory['axs_y'], observatory['axs_x']))

				for kk, tmp_file in enumerate(chk_tmp_files):
					master_tmp, master_tmp_head = fits.getdata(data_dirs['master_tmp'] + tmp_file, header=True)
					hold_data[kk] = master_tmp
					master_header = master_tmp_head

			# median combine the mini-images into one large image
			master = np.median(hold_data, axis=0)

			# update the header
			master_header['MAST_COMB'] = 'median'
			master_header['NUM_MAST'] = nfiles

			# write the image out to the master directory
			fits.writeto(data_dirs['master'] + file_name, master, master_header, overwrite=True)

		else:
			master, master_header = fits.getdata(data_dirs['master'] + file_name, header=True)

		return master, master_header

	@staticmethod
	def pull_master():
		'''This script will generate the master file and photometry file for data reduction.

		:return master - The master frame
		:return star_list - The list of kernel stars
		'''

		#cen_ra, cen_de, ref_img = Preproc.find_alignment_position()
		#Preproc.align_images(cen_ra, cen_de, ref_img)

		master, master_header = Master.mk_master()

		star_list = Master.master_phot(master, master_header)

		return master, star_list