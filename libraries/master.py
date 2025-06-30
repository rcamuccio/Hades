from config import Configuration
from libraries.utils import Utils
import numpy as np
import os
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.mast import Catalogs
from photutils.aperture import aperture_photometry, CircularAnnulus, CircularAperture
from photutils.centroids import centroid_sources

import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=Warning)

class Master:

	@staticmethod
	def master_phot(master, master_header):
		''' This function will generate the star list for the master frame and provide a photometry file.

		:parameter master - 
		:parameter master_header - 

		:return star_list - 
		'''

		if os.path.isfile(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_star_list.txt') == 0:

			field = str(Configuration.RA) + ' ' + str(Configuration.DEC)

			columns = ['field_id', 'source_id', 'ra', 'dec', 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'teff_val', 'parallax', 'parallax_error', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error']

			Utils.log('Querying MAST for all stars within the field ' + str(Configuration.FIELD), 'info')
			catalog_data = Catalogs.query_region(field, radius=Configuration.SEARCH_DIST/1.5, catalog='Gaia').to_pandas()
			Utils.log('Query finished. ' + str(len(catalog_data)) + ' stars found.', 'info')

			catalog_data['field_id'] = Configuration.FIELD

			star_list = catalog_data[columns]

			w = WCS(master_header)
			ra = star_list.ra.to_numpy()
			dec = star_list.dec.to_numpy()

			x, y = w.all_world2pix(ra, dec, 0)

			star_list['x'] = x
			star_list['y'] = y
			star_list = star_list[(star_list.x >= 10) & (star_list.x < (Configuration.AXS_X - 10)) & (star_list.y >= 10) & (star_list.y < (Configuration.AXS_X - 10))].copy().reset_index(drop=True)

			star_list['xcen'], star_list['ycen'] = centroid_sources(master, star_list.x.to_numpy(), star_list.y.to_numpy(), box_size=5)
			bd_idxs = np.where(np.isnan(star_list.xcen) | np.isnan(star_list.ycen))
			if len(bd_idxs[0]) > 0:
				for bd_idx in bd_idxs[0]:
					star_list.loc[bd_idx, 'xcen'] = star_list.loc[bd_idx, 'x']
					star_list.loc[bd_idx, 'ycen'] = star_list.loc[bd_idx, 'y']

			positions = star_list[['xcen', 'ycen']].copy().reset_index(drop=True)

			aperture = CircularAperture(positions, r=Configuration.APER_SIZE)
			annulus_aperture = CircularAnnulus(positions, r_in=Configuration.ANNULI_INNER, r_out=Configuration.ANNULI_OUTER)
			apers = [aperture, annulus_aperture]

			phot_table = aperture_photometry(master, apers, method='exact')

			sky = phot_table['aperture_sum_1'] / annulus_aperture.area

			flux = np.array(phot_table['aperture_sum_0'] - sky * aperture.area)

			flux_er = np.sqrt((phot_table['aperture_sum_0']))

			mag = 25. - 2.5*np.log10(flux)
			mag_er = (np.log(10.) / 2.5) * (flux_er / flux)

			star_list['master_mag'] = mag
			star_list['master_mag_er'] = mag_er
			star_list['master_flux'] = flux
			star_list['master_flux_er'] = flux_er
			star_list['sky'] = sky
			star_list = star_list.reset_index()
			star_list = star_list.rename(columns={'index':'star_id'})

			star_list.to_csv(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_star_list.txt', sep=' ', index=False)

		else:
			star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_star_list.txt', delimiter=' ', header=0)

		return star_list

	@staticmethod
	def mk_master(combine_type='median'):
		''' This function will make the master frame that will be used for differencing.

		:parameter combine_type - 

		:return master -
		:return master_header -
		'''

		file_name = Configuration.FIELD + '_master' + Configuration.FILE_EXTENSION

		if os.path.isfile(Configuration.MASTER_DIRECTORY + file_name) == 0:

			image_list, dates = Utils.get_all_files_per_field(Configuration.CLEAN_DIRECTORY, Configuration.FIELD, Configuration.FILE_EXTENSION)

			nfiles = len(image_list)

			if combine_type == 'mean':

				Utils.log('Generating the master frame from multiple files using a mean combination. There are ' + str(nfiles) + ' images to combine.', 'info')

				for kk in range(0, nfiles):

					img_tmp = fits.getdata(image_list[kk])

					if kk == 0:
						master = img_tmp
					else:
						master += img_tmp

				master /= nfiles

				master_header = fits.getheader(image_list[0])

				master_header['COMB'] = 'mean'
				master_header['NUM_COMB'] = nfiles

				fits.writeto(Configuration.MASTER_DIRECTORY + file_name, master, master_header, overwrite=True)

			elif combine_type == 'median':

				nbulk = 20

				full_bulk = nfiles // nbulk
				part_bulk = nfiles % nbulk

				if part_bulk > 0:
					hold_bulk = full_bulk + 1
				else:
					hold_bulk = full_bulk

				hold_data = np.ndarray(shape=(hold_bulk, Configuration.AXS_Y, Configuration.AXS_X))

				Utils.log('Generating a master frame from multiple files in bulks of ' + str(nbulk) + ' images. There are ' + str(nfiles) + ' images to combine, which means there should be ' + str(hold_bulk) + ' mini-files to median combine.', 'info')

				for kk in range(0, hold_bulk):

					if kk < full_bulk:
						block_hold = np.ndarray(shape=(nbulk, Configuration.AXS_Y, Configuration.AXS_X))
						mx_index = nbulk
					else:
						block_hold = np.ndarray(shape=(part_bulk, Configuration.AXS_Y, Configuration.AXS_X))
						mx_index = part_bulk

					loop_start = kk * nbulk
					idx_cnt = 0

					Utils.log('Making mini-file ' + str(kk) + '.', 'info')

					for jj in range(loop_start, mx_index + loop_start):

						master_tmp, master_head = fits.getdata(image_list[jj], header=True)

						if (kk == 0) & (jj == 0):
							block_hold[idx_cnt] = master_tmp
							master_hold = master_tmp
							master_header = master_head
						else:
							tmp = hcongrid(master_tmp, master_head, master_header)
							block_hold[idx_cnt] = tmp

						idx_cnt += 1

					hold_data[kk] = np.median(block_hold, axis=0)

				master = np.median(hold_data, axis=0)

				master_header['COMB'] = 'median'
				master_header['NUM_COMB'] = nfiles

				fits.writeto(Configuration.MASTER_DIRECTORY + file_name, master, master_header, overwrite=True)

		else:
			master, master_header = fits.getdata(Configuration.MASTER_DIRECTORY + file_name, header=True)

		return master, master_header

	@staticmethod
	def pull_master():
		''' This script will generate the master file and photometry file for the data reduction.

		:return master - 
		:return star_list - 
		'''

		master, master_header = Master.mk_master(combine_type='median')

		star_list = Master.master_phot(master, master_header)

		return master, star_list

