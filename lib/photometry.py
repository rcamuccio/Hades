from config import Configuration
from lib.preprocessing import Preproc
from lib.utilities import Utils
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from photutils.aperture import aperture_photometry, ApertureStats, CircularAperture, CircularAnnulus
from photutils.centroids import centroid_sources
from photutils.detection import DAOStarFinder
import numpy as np
import os
import pandas as pd
import subprocess
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=Warning)

class Photometry:

	@staticmethod
	def extract_sources(path):

		Utils.log('Extracting sources for ' + str(path) + '.', 'info')
		img = fits.open(path)
		img_data = img[0].data
		img_header = img[0].header

		wcs = WCS(img_header)

		#mask, boxes = Preproc.mk_mask(img_data)

		mn, md, sd = sigma_clipped_stats(img_data, sigma=Configuration.THRESHOLD)
		#mn, md, sd = sigma_clipped_stats(img_data, mask=mask, sigma=Configuration.THRESHOLD)

		daofind = DAOStarFinder(fwhm=Configuration.FWHM, threshold=Configuration.THRESHOLD * sd)

		table = daofind(img_data)
		#table = daofind(img_data, mask=mask)

		ra_list = []
		de_list = []

		for ln in table:
			x = ln['xcentroid']
			y = ln['ycentroid']

			sky_pos = pixel_to_skycoord(x, y, wcs=wcs)

			ra = sky_pos.ra.deg
			de = sky_pos.dec.deg

			ra_list.append(ra)
			de_list.append(de)

		table['ra'] = ra_list
		table['dec'] = de_list

		Utils.log('Writing extracted catalog to file.', 'info')
		table.write('table-source.cat', format='ascii.fixed_width', overwrite=True)

		return table

	@staticmethod
	def match_catalogs(source_table, query_table):

		coo_source = SkyCoord(source_table['ra'] * u.deg, source_table['dec'] * u.deg)

		if (Configuration.QUERY == 'gaia-cone') or (Configuration.QUERY == 'gaia-square'):
			coo_query = SkyCoord(query_table['ra'], query_table['dec'])

		idx, d2d, d3d = coo_query.match_to_catalog_sky(coo_source)
		t = d2d < 1 * u.arcsec

		query_table['sep_flag'] = t
		query_table['idx'] = idx

		ct = query_table.dtype
		names, dtypes = zip(*ct.descr)
		match_table = Table(names=names, dtype=dtypes)

		for row, sep_flag in enumerate(query_table['sep_flag']):
			if sep_flag:
				match_table.add_row(query_table[row].as_void())

		match_table.write('table-match.cat', format='ascii.fixed_width', overwrite=True)

		return match_table

	@staticmethod
	def add_variable_list(star_list, master_header):
		'''This function reads in the known variable/transient list for the star field. It will get master frame photometry and x/y pixel positions so the stars can be written to file.

		:parameter star_list - The data frame with the original star list
		:parameter master_header - The header of the master frame

		:return star_list - The star list with the variables/transient included
		'''

		observatory = Utils.config_observatory(Configuration.OBSERVATORY)
		output_dirs = Utils.config_output_dir()

		# add a column to link in the star list
		star_list['var_id'] = '--'
		star_list['var_type'] = '--'
		star_list['var_period'] = 0
		star_list['object_type'] = 'star'

		# get the variable list
		known = pd.read_csv(output_dirs['sources'] + Configuration.FIELD + '_known_objects.csv', sep=',')

		# get the variable/transient list
		#known = pd.read_csv(Configuration.MASTER_DIRECTORY + '/known_objects/' + Configuration.FIELD + '_known_objects.csv', sep=',')

		# update the ra and de in the known data frame to be in degrees
		#known['ra'] = known.apply(lambda x: ((float(x['coords'].split(' ')[0]) / 24) + (float(x['coords'].split(' ')[1]) / 60 / 24) + (float(x['coords'].split(' ')[2]) / 60 / 60 / 24)) * 360, axis=1)

		#known['dec'] = known.apply(lambda x: float(x['coords'].split(' ')[3]) - (float(x['coords'].split(' ')[4]) / 60) - (float(x['coords'].split(' ')[5]) / 60 / 60) if float(x['coords'].split(' ')[3]) < 0 else float(x['coords'].split(' ')[3]) + (float(x['coords'].split(' ')[4]) / 60) + (float(x['coords'].split(' ')[5]) / 60 / 60), axis=1)

		# get the header file and convert to x/y pixel positions
		w = WCS(master_header)
		ra = known.ra.to_numpy()
		de = known.dec.to_numpy()
		x, y = w.all_world2pix(ra, de, 0)

		# add the x/y to the star data frame
		known['x'] = x
		known['y'] = y

		# add the variable ID to the star list so it can be linked to the table
		for idx, row in known.iterrows():
			dist = np.min(np.sqrt((star_list.x - row.x) ** 2 - (star_list.y - row.y) ** 2))

			try:
				nme_chk = star_list[star_list.source_id == int(row.source_id)].index.values[0]
			except:
				nme_chk = -99

			if nme_chk > 0:
				star_list.loc[nme_chk, 'var_id'] = row.source_id
				star_list.loc[nme_chk, 'var_type'] = row.var_type
				star_list.loc[nme_chk, 'var_period'] = row.var_period
				star_list.loc[nme_chk, 'object_type'] = row.object_type

			elif (dist < 5. / observatory['pixel_scale']) & (nme_chk < 0):
				min_pos = np.argmin(np.sqrt((star_list.x - row.x) ** 2 - (star_list.y - row.y) ** 2))

		return star_list

	@staticmethod
	def combine_flux_files(star_list):
		'''This function deconstructs the flux files for each light curve.

		:parameter star_list - The star list to be used for photometry

		:return - Nothing is returned, but the raw files are saved
		'''

		# pull in the star list for photometry
		#star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_star_list.txt', delimiter=' ', header=0)

		observatory = Utils.config_observatory(Configuration.OBSERVATORY)
		data_dirs = Utils.config_data_dir(create=False)

		# get the flux files to read in
		files, dates = Utils.get_all_files_per_field(data_dirs['flux'], Configuration.FIELD, 'diff', '.flux')

		nfiles = len(files)
		nstars = len(star_list)
		Utils.log(str(nfiles) + ' flux files found for ' + Configuration.FIELD + '.', 'info')

		# make holders for the light curves
		jd = np.zeros(nfiles)
		mag = np.zeros((nstars, nfiles))
		err = np.zeros((nstars, nfiles))
		err_scl = np.zeros((nstars, nfiles))
		trd = np.zeros((nstars, nfiles))
		zpt = np.zeros((nstars, nfiles))
		sky = np.zeros((nstars, nfiles))
		bkg = np.zeros((nstars, nfiles))

		for idy, file in enumerate(files):
			# read in the data frame with the flux information
			img_flux = pd.read_csv(file, header=0)

			if idy == 0:
				star_list['tap'] = img_flux['tap'].to_numpy()
				star_list['object_type'] = img_flux['object_type'].to_numpy()

			# set the data to the numpy array
			jd[idy] = img_flux.loc[0, 'jd']
			mag[:, idy] = img_flux['mag'].to_numpy()
			err[:, idy] = img_flux['mag_err'].to_numpy()
			err_scl[:, idy] = img_flux['mag_err'].to_numpy()
			zpt[:, idy] = img_flux['zpt'].to_numpy()
			sky[:, idy] = img_flux['sky'].to_numpy()
			bkg[:, idy] = img_flux['bkg'].to_numpy()

			if (idy % 100 == 0) & (idy > 0):
				Utils.log('100 flux files read. ' + str(nfiles - idy - 1) + ' files remain.', 'info')

		# if there is an object (like 47 Tuc) in the field then block it
		tuc = False
		if tuc:
			star_list['bd_stars'] = np.where((star_list['xcen'] > 4300) & (star_list['xcen'] < 9300) & (star_list['ycen'] > 3600) & (star_list['ycen'] < 8200), 1, 0)
		else:
			star_list['bd_stars'] = np.where((star_list['xcen'] > 12000) & (star_list['xcen'] < 0) & (star_list['ycen'] > 12000) & (star_list['ycen'] < 0), 1, 0)

		# convert the master frame magnitude to e/s
		star_list['master_mag'] = star_list['master_mag'] + 2.5 * np.log10(observatory['exp_light'])

		src_id = star_list.source_id.to_numpy()

		for idy, row in star_list.iterrows():
			# get the distance to all stars
			dd = np.sqrt((row.xcen - star_list.xcen.to_numpy()) ** 2 + (row.ycen - star_list.ycen.to_numpy()) ** 2)

			# get the difference in magnitude
			dmag = np.abs(row.master_mag - star_list.master_mag.to_numpy())

			# only get nearby stars of similar magnitudes
			#vv = np.argwhere((dd < 500) & (dd > 0) & (dmag > 0) & (dmag < .5)).reshape(-1)

			# only get nearby stars of similar magnitude, on the same chip, aren't variables, and aren't in a bad area
			vv = np.argwhere((dmag < 2) & (dd > 0) & (star_list['tap'] == row.tap) & (star_list['bd_stars'] == 0) & (star_list['object_type'] == 'star')).reshape(-1)
			vv_all = np.argwhere((dmag < 0.5) & (dd > 0) & (star_list['bd_stars'] == 0) & (star_list['object_type'] == 'star')).reshape(-1)

			# generate the trend for the trend stars
			if len(vv) > 0:
				for idz in range(len(jd)):
					_, trd[idy, idz], _ = sigma_clipped_stats(mag[vv, idz] - star_list.loc[vv].master_mag.to_numpy(), sigma=3)
			else:
				for idz in range(len(jd)):
					_, trd[idy, idz], _ = sigma_clipped_stats(mag[vv_all, idz] - star_list.loc[vv_all].master_mag.to_numpy(), sigma=3)

			# get the updated error based on similar stars
			_, _, ss = sigma_clipped_stats(mag[vv_all] - trd[idy], sigma=3, axis=1)
			m_er, _, _ = sigma_clipped_stats(err[idy], sigma=3)
			if len(ss) > 0:
				sigma_scale = m_er / np.quantile(ss, 0.01)
			else:
				sigma_scale = 1.

			# rescale the errors
			err_scl[idy] = err[idy] / sigma_scale

			#if len(vv) > 0:
				# make the trend holder and standard deviation holder
				#holder = np.zeros((len(vv), len(jd)))

				# loop through all OK stars removing outliers
				#for idz in range(len(vv)):
					#_, mdn, _ = sigma_clipped_stats(mag[vv[idz], :], sigma=2.5)
					#holder[idz, :] = mag[vv[idz], :] - mdn

				#for idz in range(len(jd)):
					#_, trd[idy, idz], _ = sigma_clipped_stats(holder[:, idz], sigma=2)

			# Update the log
			if (idy % 1000 == 0) & (idy > 0):
				Utils.log('1000 stars had trends found. ' + str(nstars - idy - 1) + ' stars remain.', 'info')

		# write out the light curve data
		Photometry.write_light_curves(nstars, jd, mag, err, err_scl, trd, zpt, sky, bkg, src_id)

	@staticmethod
	def generate_flux_files():
		'''This function will generate flux files for all stars in a given star list.

		:return - Nothing is returned, but light curve files are generated for all stars
		'''

		data_dirs = Utils.config_data_dir(create=False)

		# pull in the star list for photometry
		star_list = pd.read_csv(data_dirs['master'] + Configuration.FIELD + '_star_list.csv', delimiter=' ', header=0)

		# get the image list for flux extraction
		files, dates = Utils.get_all_files_per_field(data_dirs['differenced'], Configuration.FIELD, 'diff', Configuration.FILE_EXTENSION)

		# make the output directories for the flux files
		output_dirs = []
		for dte in dates:
			output_dirs.append(data_dirs['flux'] + dte)
			output_dirs.append(data_dirs['flux'] + dte + '/' + Configuration.FIELD)
		Utils.create_directories(output_dirs)

		# begin the algorithm to run photometry
		for idx, file in enumerate(files):
			file_nme = file.split(Configuration.FILE_EXTENSION)[0] + '.flux'
			file_nme = file_nme.replace('/diff/', '/flux/')
			tmp_nme = file_nme.split('/')[-1]

			# check if flux file already exists
			if os.path.isfile(file_nme):
				Utils.log('Flux file ' + tmp_nme + ' found. Skipping.', 'info')

			# run aperture photometry on the star list
			else:
				Utils.log('Extracting flux from ' + tmp_nme + '.', 'info')
				Photometry.single_frame_aperture_photometry(star_list, file, file_nme)

		Utils.log('Flux extraction complete for ' + Configuration.FIELD + '.', 'info')

	@staticmethod
	def mk_raw_lightcurves():
		'''This function will create individual raw light curve files for each star in the specific star list.

		:return - Nothing is returned, but each light curve is generated
		'''

		data_dirs = Utils.config_data_dir(create=False)

		# pull in the star list for the photometry
		star_list = pd.read_csv(data_dirs['master'] + Configuration.FIELD + '_star_list.csv', delimiter=' ', header=0)

		# combine the flux from the flux files and write the raw light curves
		Photometry.combine_flux_files(star_list)

	@staticmethod
	def mk_sextractor_conv(path):

		ln1 = 'CONV NORM\n'
		ln2 = '# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.\n'
		ln3 = '1 2 1\n'
		ln4 = '2 4 2\n'
		ln5 = '1 2 1'

		lns = [ln1, ln2, ln3, ln4, ln5]

		default_conv = open(path, 'w')
		for ln in lns:
			default_conv.write(ln)
		default_conv.close()

	@staticmethod
	def mk_sextractor_nnw(path):

		ln1 = 'NNW\n'
		ln2 = '# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)\n'
		ln3 = '# inputs:  9 for profile parameters + 1 for seeing.\n'
		ln4 = '# outputs: ``Stellarity index'' (0.0 to 1.0)\n'
		ln5 = '# Seeing FWHM range: from 0.025 to 5.5'' (images must have 1.5 < FWHM < 5 pixels)\n'
		ln6 = '# Optimized for Moffat profiles with 2<= beta <= 4.\n\n'
		ln7 = ' 3 10 10 1\n\n'
		ln8 = '-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01 -5.22453e-01  4.61342e-02  8.31957e-01  2.15505e+00  2.64769e-01\n'
		ln9 = ' 3.03477e+00  2.69561e+00  3.16188e+00  3.34497e+00  3.51885e+00  3.65570e+00  3.74856e+00  3.84541e+00  4.22811e+00  3.27734e+00\n\n'
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

		lns = [ln1, ln2, ln3, ln4, ln5, ln6, ln7, ln8, ln9, ln10, ln11, ln12, ln13, ln14, ln15, ln16, ln17, ln18, ln19, ln20, ln21, ln22]

		default_nnw = open(path, 'w')
		for ln in lns:
			default_nnw.write(ln)
		default_nnw.close()

	@staticmethod
	def mk_sextractor_param(path):

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

		default_param = open(path, 'w')
		for ln in lns:
			default_param.write(ln)
		default_param.close()

	@staticmethod
	def mk_sextractor_sex(path):

		observatory = Utils.config_observatory(Configuration.OBSERVATORY)

		# Default configuration file for SExtractor 2.12.4
		# EB 2010-10-10

		catalog_name = 'CATALOG_NAME ' + str(Configuration.FIELD) + '.cat\n'
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
		satur_level = 'SATUR_LEVEL ' + str(observatory['peakmax']) + '\n'
		satur_key = 'SATUR_KEY SATURATE\n'
		mag_zeropoint = 'MAG_ZEROPOINT 0.0\n'
		mag_gamma = 'MAG_GAMMA 4.0\n'
		gain = 'GAIN ' + str(Configuration.GAIN) + '\n'
		gain_key = 'GAIN_KEY GAIN\n'
		pixel_scale = 'PIXEL_SCALE ' + str(observatory['pixel_scale']) + '\n'
		seeing_fwhm = 'SEEING_FWHM ' + str(observatory['seeing']) + '\n'
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
		verbose_type = 'VERBOSE_TYPE NORMAL\n'
		header_suffix = 'HEADER_SUFFIX .head\n'
		write_xml = 'WRITE_XML N\n'
		xml_name = 'XML_NAME sex.xml\n'
		xsl_url = 'XSL_URL file:///usr/local/share/sextractor/sextractor.xsl\n'

		lns = [catalog_name, catalog_type, parameters_name, detect_type, detect_minarea, detect_thresh, analysis_thresh, detection_filter, detection_filter_name, deblend_nthresh, deblend_mincont, clean, clean_param, weight_type, weight_image, flag_image, flag_type, phot_apertures, phot_autoparams, phot_petroparams, phot_autoapers,satur_level, satur_key, mag_zeropoint, mag_gamma, gain, gain_key, pixel_scale, seeing_fwhm, starnnw_name, back_type, back_value, back_size, back_filtersize, checkimage_type, checkimage_name, memory_objstack, memory_pixstack, memory_bufsize, assoc_name, assoc_data, assoc_params, assoc_radius, assoc_type, assocselec_type, verbose_type, header_suffix, write_xml, xml_name,xsl_url]

		default_sex = open(path, 'w')
		for ln in lns:
			default_sex.write(ln)
		default_sex.close()

	@staticmethod
	def run_sextractor():

		data_dirs = Utils.config_data_dir(create=False)

		conv_path = data_dirs['master'] + '/default.conv'
		Photometry.mk_sextractor_conv(conv_path)

		nnw_path = data_dirs['master'] + '/default.nnw'
		Photometry.mk_sextractor_nnw(nnw_path)

		param_path = data_dirs['master'] + '/default.param'
		Photometry.mk_sextractor_param(param_path)

		sex_path = data_dirs['master'] + '/default.sex'
		Photometry.mk_sextractor_sex(sex_path)

		file_path = data_dirs['master'] + Configuration.FIELD + '_master' + Configuration.FILE_EXTENSION

		os.chdir(data_dirs['master'])
		subprocess.run(['source-extractor', file_path])

	@staticmethod
	def single_frame_aperture_photometry(star_list, img_name, fin_name):
		'''This function will find the subtraction stars to use for the differencing. They will be the same stars for every frame. This will help in detrending later.

		:parameter star_list - The data frame with the list of stars to use for subtraction
		:parameter img_name - The file from which to extract flux
		:parameter fin_name - The name of the output flux file

		:return - Nothing is returned, but the flux file is written
		'''

		observatory = Utils.config_observatory(Configuration.OBSERVATORY)

		# get the image for photometry
		img, header = fits.getdata(img_name, header=True)

		# get the various important header information
		time = Time(header['DATE'], format='isot', scale='utc')
		jd = time.jd

		# get the stellar positions from the master frame
		positions = np.transpose((star_list['xcen'], star_list['ycen']))

		# set the apertures to each stellar position
		aperture = CircularAperture(positions, r=Configuration.APER_SIZE)
		aperture_area = aperture.area
		annulus = CircularAnnulus(positions, r_in=Configuration.ANNULI_INNER, r_out=Configuration.ANNULI_OUTER)

		# get the background statistics
		aperstats = ApertureStats(img, annulus)
		bkg_mean = aperstats.mean
		bkg_total = bkg_mean * aperture_area

		# run photometry to get the data table
		phot_table = aperture_photometry(img, aperture, method='exact')

		# extract flux from the table
		star_flux = np.array(phot_table['aperture_sum']) * observatory['gain']

		# calculate the expected photometric error
		star_err = np.abs(star_flux.astype(float) + star_list['master_flux'].to_numpy().astype(float))
		bkg_err = np.abs(header['sky'] + bkg_mean) * aperture_area * observatory['gain']
		#star_error = np.abs(star_flux)
		#bkg_error = np.sqrt((header['sky'] * aperture_area) ** 2 + bkg_total ** 2) * Configuration.GAIN

		# combine sky and signal errors in quadrature
		star_flux_err = np.sqrt(star_err + bkg_err)
		#star_flux_snr = np.abs(star_flux / star_flux_err)

		# combine the fluxes
		flux = star_flux.astype(float) + star_list['master_flux'].to_numpy().astype(float)
		flux_err = np.sqrt(star_flux_err.astype(float) ** 2 + star_list['master_flux_err'].to_numpy().astype(float) ** 2)

		# convert to magnitude
		mag = 25. - 2.5 * np.log10(flux)
		mag_nbkg = 25. - 2.5 * np.log10(flux - bkg_total)
		mag_err = (np.log(10.) / 2.5) * (flux_err / flux)

		# initialize the master magnitude of all frames
		m_mag = star_list['master_mag'].to_numpy()

		# get the delta magnitude between the science and master frames
		dmag = mag[~np.isnan(mag)] - m_mag[~np.isnan(mag)]

		# initialize the offset vector
		off = np.zeros(len(mag))

		# set up the holders for interpolation
		f_mags = np.arange(6, 16) + 0.5
		n_mags = np.zeros(len(f_mags))

		# loop through each tap to find the zero point offset
		for tap in range(1, 17):
			for mag_idx, m_lw in enumerate(f_mags):
				# get the zero point using non-nan stars in the tap between the magnitude range
				dmag_bin = dmag[(mag[~np.isnan(mag)] > m_lw) & (mag[~np.isnan(mag)] < m_lw + 1) & (star_list[~np.isnan(mag)].tap.to_numpy() == tap) & (star_list[~np.isnan(mag)].object_type.to_numpy() == 'star')]

				#dmag_bin = dmag[(mag[~np.isnan(mag)] > m_lw) & (mag[~np.isnan(mag)] < m_lw + 1) & (star_list[~np.isnan(mag)].tap.to_numpy() == tap)]

				# sigma clip outliers
				dmag_mn, dmag_md, dmag_sg = sigma_clipped_stats(np.array(dmag_bin, dtype=float), sigma=2)
				n_mags[mag_idx] = dmag_md

			# interpolate to all magnitudes
			off[star_list.tap.to_numpy() == tap] = np.interp(mag[star_list.tap.to_numpy() == tap], f_mags, n_mags)

		# correct the magnitudes for exposure time
		mag = mag + 2.5 * np.log10(observatory['exp_light'])

		# replace nans with -9.999999
		off = np.where(np.isnan(off), -9.999999, off)
		mag = np.where(np.isnan(mag), -9.999999, mag)

		# generate the final flux file
		flux_file = star_list.copy().reset_index(drop=True)
		flux_file['flux'] = flux
		flux_file['flux_err'] = flux_err
		flux_file['mag'] = mag
		flux_file['mag_nbkg'] = mag_nbkg
		flux_file['mag_err'] = mag_err
		flux_file['sky'] = header['SKY']
		flux_file['bkg'] = bkg_mean
		flux_file['jd'] = jd
		flux_file['zpt'] = off
		flux_file['exp_time'] = observatory['exp_light']

		flux_file.to_csv(fin_name, header=True, index=False)

	@staticmethod
	def write_light_curves(nstars, jd, mag, err, err_scl, trd, zpt, sky, bkg, src_id):
		'''This function writes the flux columns to light curves for each source ID.

		:parameter nstars - The number of stars to make light curves
		:parameter jd - The numpy array of julian dates (one per file)
		:parameter mag - The magnitudes for each stars at each time
		:parameter err - The photometric error for each star at each time
		:parameter err_scl - The photometric error for each star with scaling
		:parameter trd - The trend from nearby similar magnitude stars
		:parameter zpt - The zeropoint using all stars in the frame
		:parameter sky - The median sky background
		:parameter bkg - The local background for the source
		:parameter src_id - An nstars-length array of source IDs

		:return - Nothing is returned, but the light curve files are written
		'''

		data_dirs = Utils.config_data_dir(create=False)

		Utils.log('Writing light curves for ' + str(nstars) + ' stars.', 'info')

		# initialize the light curve data frame
		lc = pd.DataFrame(columns=['jd', 'mag', 'err', 'trd', 'zpt', 'sky', 'bkg'])

		for idx in range(0, nstars):
			star_id = str(src_id[idx])

			# add the time, magnitude, and error to the data frame
			lc['jd'] = np.around(jd, decimals=6)
			lc['mag'] = np.around(mag[idx, :], decimals=6)
			lc['err'] = np.around(err_scl[idx, :], decimals=6)
			lc['trd'] = np.around(trd[idx, :], decimals=6)
			lc['org_err'] = np.around(err[idx, :], decimals=6)
			lc['zpt'] = np.around(zpt[idx, :], decimals=6)
			lc['sky'] = np.around(sky[idx, :], decimals=6)
			lc['bkg'] = np.around(bkg[idx, :], decimals=6)

			# make sure the data are in order
			lc = lc.sort_values(by = 'jd').reset_index(drop=True)
			lc['err'] = np.where(lc['mag'] < 0, -9.999999, lc['err'])
			lc['org_err'] = np.where(lc['mag'] < 0, -9.999999, lc['org_err'])
			lc['trd'] = np.where(lc['mag'] < 0, -9.999999, lc['trd'])
			lc['zpt'] = np.where(lc['mag'] < 0, -9.999999, lc['zpt'])
			lc['sky'] = np.where(lc['mag'] < 0, -9.999999, lc['sky'])
			lc['bkg'] = np.where(lc['mag'] < 0, -9.999999, lc['bkg'])
			lc['mag'] = np.where(lc['mag'] < 0, -9.999999, lc['mag'])

			# write the new file
			lc[['jd', 'mag', 'err', 'org_err', 'zpt', 'sky', 'bkg']].to_csv(data_dirs['lightcurve_field'] + Configuration.FIELD + '_' + star_id + '.lc', sep=' ', index=False, na_rep='-9.999999')

			if (idx > 0) & (idx / 10000 % 1 == 0):
				Utils.log('10000 stars have had their light curves written. ' + str(nstars - idx - 1) + ' stars remain.', 'info')

		Utils.log('All light curves written.', 'info')