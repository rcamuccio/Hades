from config import Configuration
from libraries.photometry import Photometry
from libraries.utils import Utils
import numpy as np
import os
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from photutils.aperture import aperture_photometry, CircularAperture, CircularAnnulus
from photutils.centroids import centroid_sources

import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=Warning)

class Lightcurves:

	@staticmethod
	def generate_flux_files(star_list):
		''' This function will generate light curves for all of the stars in a given star list.

		:parameter star_list - The star list generated from the master frame

		:return - Nothing is returned, but light curve files are generated for all stars
		'''

		files, dates = Utils.get_all_files_per_field(Configuration.DIFFERENCED_DIRECTORY, Configuration.FIELD, Configuration.FILE_EXTENSION)

		if Configuration.TRANSIENT_LC == 'Y':
			trans = pd.Series(index=star_list.columns.to_list())
			master, master_head = fits.getdata(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_master.fits', header=True)
			w = WCS(master_head)
			ra = Configuration.TRANSIENT_RA
			dec = Configuration.TRANSIENT_DEC

			x, y = w.all_world2pix(ra, dec, 0)

			trans['x'] = x
			trans['y'] = y

			xcen, ycen = centroid_sources(master, trans.x.item(), trans.y.item(), box_size=5)
			if np.isnan(xcen) | np.isnan(ycen):
				trans['xcen'] = x
				trans['ycen'] = y
			else:
				trans['xcen'] = xcen[0]
				trans['ycen'] = ycen[0]

			trans['star_id'] = Configuration.TRANSIENT_NAME
			trans['toros_field_id'] = Configuration.FIELD
			trans['source_id'] = '-999'
			trans['ra'] = Configuration.TRANSIENT_RA
			trans['dec'] = Configuration.TRANSIENT_DEC

			positions = np.transpose((trans['x'], trans['y']))

			aperture = CircularAperture(positions, r=Configuration.APER_SIZE)
			aperture_annulus = CircularAnnulus(positions, r_in=Configuration.ANNULI_INNER, r_out=Configuration.ANNULI_OUTER)
			apers = [aperture, aperture_annulus]
			phot_table = aperture_photometry(master, apers, method='exact')

			sky = np.median(master)
			flux = np.array(phot_table['aperture_sum_0'] - (sky * aperture.area))
			trans['master_flux'] = np.around(flux, decimals=6)
			star_error = np.sqrt(np.abs(phot_table['aperture_sum_0']))
			flux_er = np.array(np.sqrt(star_error**2))
			trans['master_flux_er'] = np.around(flux_er, decimals=6)

			trans['master_mag'] = np.around(25. - 2.5*np.log10(flux), decimals=6)
			trans['master_mag_er'] = np.around((np.log(10.) / 2.5) * (flux_er / flux), decimals=6)

			star_list.loc[len(star_list)] = trans

		for idx, file in enumerate(files):

			fin_nme = file.split('.fits')[0] + '.flux'

			if os.path.isfile(fin_nme) == 0:
				Utils.log('Flux file ' + fin_nme + ' found. Skipping...', 'info')

			if os.path.isfile(fin_nme) == 1:
				Utils.log('Working to extract flux from ' + file + '.', 'info')
				Photometry.single_frame_aperture_photometry(star_list, file, fin_nme)

		Utils.log('Flux extraction complete for ' + Configuration.FIELD + '.', 'info')

	@staticmethod
	def mk_raw_lightcurves(star_list):
		''' This function will create the individual raw light curve files for each star in the specific star list.

		:parameter star_list - 

		:return - Nothing is returned, but each light curve is output
		'''

		Photometry.combine_flux_files(star_list)