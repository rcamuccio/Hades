from config import Configuration
from libraries.preprocessing import Preprocessing
from libraries.utils import Utils
import numpy as np
import os
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.apertures import aperture_photometry, CircularAperture, CircularAnnulus

class Difference:

	@staticmethod
	def difference_images(star_list):
		''' This function will generate a master frame and position files, or if a single frame is chosen, then only the position file is generated.

		:parameter star_list - A data frame with the aperture photometry from the master image

		:return - Nothing is returned, but the images are differenced
		'''

		files, dates = Utils.get_all_files_per_field(Configuration.CLEAN_DIRECTORY, Configuration.FIELD, Configuration.FILE_EXTENSION)

		nfiles = len(files)

		master, master_header = fits.getdata(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_master' + Configuration.FILE_EXTENSION, header=True)

		Difference.prep_ois(master, master_header)

		for ii in range(0, nfiles):

			fin_nme = Preprocessing.mk_nme(files[ii], 'Y', 'N', 'N', 'N', 'N')

			if os.path.isfile(fin_nme) == 1:
				Utils.log('File ' + files[ii] + ' found. Skipping...', 'info')

			if os.path.isfile(fin_nme) == 0:
				Utils.log('Working to difference file ' + files[ii] + '.', 'info')
				Difference.diff_img(star_list, files[ii], fin_nme)

		Utils.log('Differencing complete for ' + Configuration.FIELD + '.', 'info')

	@staticmethod
	def diff_img(star_list, file, out_name):
		''' This function will check for and determine the reference stars. It will then difference the image.

		:parameter file - The file name to difference
		:parameter out_name - The final file name

		:return - Nothing is returned, but the image is differenced
		'''

		org_img, org_header = fits.getdata(file, header=True)
		img_sky_mean, img_sky_median, img_sky_std = sigma_clipped_stats(org_img, sigma=3.0)

		master_header = fits.getheader(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_master.fits')

		img_sbkg = org_img - img_sky_median
		img_align = hcongrid(img_sbkg, org_header, master_header)
		org_header['WCSAXES'] = master_header['WCSAXES']
		org_header['CRPIX1'] = master_header['CRPIX1']
		org_header['CRPIX2'] = master_header['CRPIX2']
		org_header['PC1_1'] = master_header['PC1_1']
		org_header['PC1_2'] = master_header['PC1_2']
		org_header['PC2_1'] = master_header['PC2_1']
		org_header['PC2_2'] = master_header['PC2_2']
		org_header['CDELT1'] = master_header['CDELT1']
		org_header['CDELT2'] = master_header['CDELT2']
		org_header['CUNIT1'] = master_header['CUNIT1']
		org_header['CUNIT2'] = master_header['CUNIT2']
		org_header['CTYPE1'] = master_header['CTYPE1']
		org_header['CTYPE2'] = master_header['CTYPE2']
		org_header['CRVAL1'] = master_header['CRVAL1']
		org_header['CRVAL2'] = master_header['CRVAL2']
		org_header['LONPOLE'] = master_header['LONPOLE']
		org_header['LATPOLE'] = master_header['LATPOLE']
		org_header['MJDREF'] = master_header['MJDREF']
		org_header['RADESYS'] = master_header['RADESYS']
		org_header['ALIGNED'] = 'Y'

		fits.writeto(Configuration.CODE_DIFFERENCE_DIRECTORY + 'img.fits', img_align, org_header, overwrite=True)

		kernel_stars = Difference.find_subtraction_stars_img(org_img, star_list)

		Difference.ois_difference(out_name, org_header)

	@staticmethod
	def find_subtraction_stars_img(img, star_list):
		''' This function will find the subtraction stars to use for the differencing. They will be the same stars for every frame. This will help in detrending later.

		:parameter star_list - The data frame with the list of stars to use for the subtraction

		:return diff_list - 
		'''

		Utils.log('Finding stars for kernel from the star list.', 'info')

		positions = np.transpose((star_list['xcen'], star_list['ycen']))

		aperture = CircularAperture(positions, r=Configuration.APER_SIZE)
		aperture_annulus = CircularAnnulus(positions, r_in=Configuration.ANNULI_INNER, r_out=Configuration.ANNULI_OUTER)
		apers = [aperture, aperture_annulus]

		phot_table = aperture_photometry(img, apers, method='exact')

		sky = phot_table['aperture_sum_1'] / aperture_annulus.area

		flux = np.array(phot_table['aperture_sum_0'] - (sky * aperture.area)) / Configuration.EXP_TIME

		star_error = np.sqrt((phot_table['aperture_sum_0'] - (sky * aperture.area)) * Configuration.GAIN)
		sky_error = np.sqrt(aperture.area * sky * Configuration.GAIN)

		flux_er = np.array(np.sqrt(star_error**2 + sky_error**2))

		mag = 25 - 2.5*np.log10(flux)
		mag_er = (np.log(10.) / 2.5) * (flux_er / flux)

		diff_list = star_list.copy().reset_index(drop=True)
		diff_list['min_dist'] = diff_list.apply(lambda x: np.sort(np.sqrt((x['x'] - diff_list['x'])**2 + x['y'] - diff_list['y'])**2)[1], axis=1)
		dist_cut = 2*Configuration.STMP + 1
		diff_list['dmag'] = np.abs(diff_list['master_mag'].to_numpy() - mag)
		mn, md, sg = sigma_clipped_stats(diff_list.dmag, sigma=2)
		mag_plus = md + sg
		mag_minus = md - sg

		diff_list = star_list[(star_list['xcen'] > Configuration.AXS_LIMIT)
							& (star_list['xcen'] < Configuration.AXS_X - Configuration.AXS_LIMIT)
							& (star_list['ycen'] > Configuration.AXS_LIMIT)
							& (star_list['ycen'] < Configuration.AXS_Y - Configuration.AXS_LIMIT)
							& (diff_list['min_dist'] > dist_cut)
							& ((diff_list['dmag'] < mag_plus) | (diff_list['dmag'] > mag_minus))].copy().reset_index(drop=True)

		if len(diff_list) > Configuration.NRSTARS:
			diff_list = diff_list.sample(n=Configuration.NRSTARS)
		else:
			Utils.log('There are not enough stars on the frame to use ' + str(Configuration.NRSTARS) + ' in the subtraction. Using all available stars.', 'info')

		diff_list['x'] = np.around(diff_list['xcen'] + 1, decimals=0)
		diff_list['y'] = np.around(diff_list['ycen'] + 1, decimals=0)

		Utils.write_txt(Configuration.CODE_DIFFERENCE_DIRECTORY + 'params.txt', 'w', '%1d %1d %1d %4d\n' % (Configuration.STMP, Configuration.KRNL, Configuration.ORDR, len(diff_list)))

		diff_list[['x', 'y']].astype(int).to_csv(Configuration.CODE_DIFFERENCE_DIRECTORY + 'refstars.txt', index=0, header=0, sep=' ')

		return diff_list

	@staticmethod
	def ois_difference(out_name, header):
		''' This function will run the c code oisdifference.

		:parameter out_name - The file name for the difference file
		:parameter header - The header of the image

		:return Nothing is returned, but the image is differenced
		'''

		Utils.log('Now starting image subtraction.', 'info')

		Utils.log('The kernel size is: ' + str(Configuration.KRNL * 2 + 1) + 'x' + str(Configuration.KRNL * 2 + 1) + '; the stamp size is: ' + str(Configuration.STMP * 2 + 1) + 'x' + str(Configuration.STMP * 2 + 1) + '; the polynomial is: ' + str(Configuration.ORDR) + ' order; and ' + str(Configuration.NRSTARS) + ' stars were used in the subtraction.', 'info')

		os.chdir(Configuration.CODE_DIFFERENCE_DIRECTORY)

		shh = os.system('./a.out')

		dimg, diff_header = fits.getdata('dimg.fits', header=True)
		header['diffed'] = 'Y'

		fits.writeto('dimg.fits', dimg, header, overwrite=True)

		os.system('mv dimg.fits ' + out_name)

		os.chdir(Configuration.WORKING_DIRECTORY)

		Utils.log('Image subtraction complete.', 'info')

	@staticmethod
	def prep_ois(master, master_header):
		''' This function will prepare the files necessary for the ois difference.

		:parameter master - The master image for differencing
		:parameter master_header - The header file for the master image

		:return - Nothing is returned, but the necessary text files are written, and the code is compiled for differencing
		'''

		os.system('cp oisdifference.c ' + Configuration.CODE_DIFFERENCE_DIRECTORY)
		os.chdir(Configuration.CODE_DIFFERENCE_DIRECTORY)
		os.system('gcc oisdifference.c -L/Users/yuw816/Development/cfitsio-4.6.2/lib -I/Users/yuw816/Development/cfitsio-4.6.2/include -lcfitsio -lm')
		os.chdir(Configuration.WORKING_DIRECTORY)

		master_sky_mean, master_sky_median, master_sky_std = sigma_clipped_stats(master, sigma=3.0)

		master_sbkg = master - master_sky_median

		fits.writeto(Configuration.CODE_DIFFERENCE_DIRECTORY + 'ref.fits', master_sbkg, master_header, overwrite=True)

		Utils.write_txt(Configuration.CODE_DIFFERENCE_DIRECTORY + 'ref.txt', 'w', 'ref.fits')
		Utils.write_txt(Configuration.CODE_DIFFERENCE_DIRECTORY + 'img.txt', 'w', 'img.fits')