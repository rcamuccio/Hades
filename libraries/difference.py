from config import Configuration
from libraries.preprocessing import Preprocessing
from libraries.utils import Utils
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from photutils.aperture import aperture_photometry, CircularAperture
import numpy as np
import os
import pandas as pd
import shutil
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=Warning)

class BigDiff:

    @staticmethod
    def difference_images():
        '''This function will generate the master frame and position files, or if a single frame is chosen, then only a position file is generated.

        :return - Nothing is returned, but the images are differenced
        '''

        # get the image list to difference
        files, dates = Utils.get_all_files_per_field(Configuration.CLEAN_DIRECTORY, Configuration.FIELD, 'clean', Configuration.FILE_EXTENSION)
        nfiles = len(files)

        # read in the master frame information
        master, master_header = fits.getdata(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_master' + Configuration.FILE_EXTENSION, header=True)

        # read in the star list for processing
        star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_star_list.txt', delimiter=' ', header=0)

        # prepare the oisdifference.c file for differencing
        BigDiff.prep_ois(master, master_header)

        # begin with the algorithm to difference the images
        for file_path in files:
            file_nme = Preprocessing.mk_nme(file_path, 'Y', 'N', 'N', 'N', 'N')

            mod_file_nme = file_nme.split('/')[-1]

            # check to see if the difference file already exists
            if os.path.isfile(file_nme):
                Utils.log('File ' + mod_file_nme + ' found. Skipping.', 'info')

            else:
                Utils.log('Differencing file ' + mod_file_nme + '.', 'info')
                BigDiff.diff_img(star_list, file_path, file_nme)

        Utils.log('Differencing complete for ' + Configuration.FIELD + '.', 'info')

        return

    @staticmethod
    def diff_img(star_list, file, out_name):
        '''This function will check for and determine reference stars. It will then difference the image.

        :parameter star_list - 
        :parameter file - The file name to difference
        :parameter out_name - The final file name

        :return - Nothing is returned, but the image is differenced
        '''

        # read in the image
        org_img, org_header = fits.getdata(file, header=True)

        # read in the master frame header to align the images
        master_header = fits.getheader(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_master' + Configuration.FILE_EXTENSION)

        # write the new image file
        img_sbkg = org_img - org_header['SKY']
        img_align = Preprocessing.align_img(img_sbkg, org_header, master_header)

        org_header['WCSAXES'] = master_header['WCSAXES']
        org_header['CRPIX1'] = master_header['CRPIX1']
        org_header['CRPIX2'] = master_header['CRPIX2']

        try:
            org_header['PC1_1'] = master_header['PC1_1']
            org_header['PC1_2'] = master_header['PC1_2']
            org_header['PC2_1'] = master_header['PC2_1']
            org_header['PC2_2'] = master_header['PC2_2']
            org_header['CDELT1'] = master_header['CDELT1']
            org_header['CDELT2'] = master_header['CDELT2']

        except KeyError as e:
            Utils.log(f"{e}... Trying with CDi_j instead.", "info")
            org_header['CD1_1'] = master_header['CD1_1']
            org_header['CD1_2'] = master_header['CD1_2']
            org_header['CD2_1'] = master_header['CD2_1']
            org_header['CD2_2'] = master_header['CD2_2']

        if WCS(master_header).has_distortion:
            org_header['A_ORDER'] = master_header['A_ORDER']
            org_header['A_0_0'] = master_header['A_0_0']
            org_header['A_0_1'] = master_header['A_0_1']
            org_header['A_0_2'] = master_header['A_0_2']
            org_header['A_1_0'] = master_header['A_1_0']
            org_header['A_1_1'] = master_header['A_1_1']
            org_header['A_2_0'] = master_header['A_2_0']
            org_header['B_ORDER'] = master_header['B_ORDER']
            org_header['B_0_0'] = master_header['B_0_0']
            org_header['B_0_1'] = master_header['B_0_1']
            org_header['B_0_2'] = master_header['B_0_2']
            org_header['B_1_0'] = master_header['B_1_0']
            org_header['B_1_1'] = master_header['B_1_1']
            org_header['B_2_0'] = master_header['B_2_0']

            org_header['AP_ORDER'] = master_header['AP_ORDER']
            org_header['AP_0_0'] = master_header['AP_0_0']
            org_header['AP_0_1'] = master_header['AP_0_1']
            org_header['AP_0_2'] = master_header['AP_0_2']
            org_header['AP_1_0'] = master_header['AP_1_0']
            org_header['AP_1_1'] = master_header['AP_1_1']
            org_header['AP_2_0'] = master_header['AP_2_0']
            org_header['BP_ORDER'] = master_header['BP_ORDER']
            org_header['BP_0_0'] = master_header['BP_0_0']
            org_header['BP_0_1'] = master_header['BP_0_1']
            org_header['BP_0_2'] = master_header['BP_0_2']
            org_header['BP_1_0'] = master_header['BP_1_0']
            org_header['BP_1_1'] = master_header['BP_1_1']
            org_header['BP_2_0'] = master_header['BP_2_0']

        org_header['CUNIT1'] = master_header['CUNIT1']
        org_header['CUNIT2'] = master_header['CUNIT2']
        org_header['CTYPE1'] = master_header['CTYPE1']
        org_header['CTYPE2'] = master_header['CTYPE2']
        org_header['CRVAL1'] = master_header['CRVAL1']
        org_header['CRVAL2'] = master_header['CRVAL2']
        org_header['LONPOLE'] = master_header['LONPOLE']
        org_header['LATPOLE'] = master_header['LATPOLE']

        try:
            org_header['MJDREF'] = master_header['MJDREF']
            org_header['RADESYS'] = master_header['RADESYS']
        except KeyError as e:
            Utils.log(f'{e}... Skipping keywords MJDREF, RADESYS.', 'info')

        org_header['ALIGNED'] = 'Y'

        img_name = 'img' + Configuration.FILE_EXTENSION
        fits.writeto(Configuration.CODE_DIFFERENCE_DIRECTORY + img_name, img_align, org_header, overwrite=True)

        # get the kernel stars for the subtraction
        nstars = BigDiff.find_subtraction_stars_img(img_align, star_list)

        BigDiff.ois_difference(file, out_name, org_header, nstars)

        return

    @staticmethod
    def find_subtraction_stars_img(img, star_list):
        '''This function will find the subtraction stars to use for differencing. They will be the same stars for every frame. This will help in detrending later.

        :parameter img - 
        :parameter star_list - 

        :return nstars - 
        '''

        Utils.log('Finding stars for kernel from the star list.', 'info')
        diff_list = star_list.copy().reset_index(drop=True)

        # check for stars based on their magnitude differences
        positions = np.transpose((diff_list['xcen'], diff_list['ycen']))
        aperture = CircularAperture(positions, r=Configuration.APER_SIZE)

        # run the photometry to get the data table
        phot_table = aperture_photometry(img, aperture, method='exact')

        # the global background was subtracted, and the master magnitude does not have exposure time corrected
        flux = np.array(phot_table['aperture_sum']) * Configuration.GAIN

        # convert to magnitude
        mag = 25. - 2.5 * np.log10(flux)

        # clip likely variables - these objects have large magnitude changes relative to the master frame
        diff_list['dmag'] = diff_list['master_mag'].to_numpy() - mag

        # the clip will be 2-sigma-clipped above and below the mean offset
        dmn, dmd, dsg = sigma_clipped_stats(diff_list.dmag, sigma=2)
        dmag_plus = dmd + dsg
        dmag_minus = dmd - dsg
        diff_list = diff_list[(diff_list['dmag'] < dmag_plus) & (diff_list['dmag'] > dmag_minus)].copy().reset_index(drop=True)

        # now check for non-crowded stars
        diff_list['prox'] = diff_list.apply(lambda x: np.sort(np.sqrt((x.xcen - diff_list.xcen) ** 2 + (x.ycen - diff_list.ycen) ** 2))[1], axis=1)
        diff_list = diff_list[diff_list.prox > (2 * Configuration.STMP + 1)].copy().reset_index(drop=True)

        # make a magnitude cut
        diff_list = diff_list[(diff_list.phot_g_mean_mag < 15) & (diff_list.master_mag > 10)].copy().reset_index(drop=True)

        if len(diff_list) > Configuration.NRSTARS:
            diff_list = diff_list.sample(n=Configuration.NRSTARS)
            nstars = Configuration.NRSTARS
        else:
            nstars = len(diff_list)
            Utils.log('Insufficient stars on frame for subtraction (' + str(Configuration.NRSTARS) + '). All available stars will be used (' + str(nstars) + '.', 'info')

        # add 1 for indexing in c vs indexing in python
        diff_list['x'] = np.around(diff_list['xcen'] + 1, decimals=0)
        diff_list['y'] = np.around(diff_list['ycen'] + 1, decimals=0)

        # write the parameter file for the C code to use
        Utils.write_txt(Configuration.CODE_DIFFERENCE_DIRECTORY + 'parms.txt', 'w', '%1d %1d %1d %4d\n' % (Configuration.STMP, Configuration.KRNL, Configuration.ORDR, nstars))

        refstars = 'refstars.txt'

        # export the differencing stars
        diff_list[['x', 'y']].astype(int).to_csv(Configuration.CODE_DIFFERENCE_DIRECTORY + refstars, index=0, header=0, sep=' ')

        return nstars

    @staticmethod
    def ois_difference(filepath, out_name, header, nstars):
        '''This function will run the c code oisdifference.

        :parameter filepath - 
        :parameter out_name - The file name for the difference file
        :parameter header - The header of the image
        :parameter nstars - The number of stars in the subtraction

        :return - Nothing is returned, but the image is differenced
        '''

        Utils.log('Starting image subtraction.', 'info')
        '''
        Utils.log('The kernel size is: ' + str(Configuration.KRNL * 2 + 1) + 'x' + str(Configuration.KRNL * 2 + 1) + '; the stamp size is: ' + str(Configuration.STMP * 2 + 1) + 'x' + str(Configuration.STMP * 2 + 1) + '; the polynomial is: ' + str(Configuration.ORDR) + ' order; and ' + str(nstars) + ' stars are being used in the subtraction.', 'info')
        '''

        # change to the directory
        os.chdir(Configuration.CODE_DIFFERENCE_DIRECTORY)

        differenced_img_name = 'dimg' + Configuration.FILE_EXTENSION
        command = './a.out'

        # run the c code
        shh = os.system(command)

        # update the header file
        dimg, diff_header = fits.getdata(differenced_img_name, header=True)
        header['diffed'] = 'Y'
        header['nstars'] = nstars

        # update the image with the new file header
        fits.writeto(differenced_img_name, dimg, header, overwrite=True)

        # move the differenced file to the difference directory
        shutil.move(differenced_img_name, out_name)

        # change back to the working directory
        os.chdir(Configuration.OUTPUT_DIRECTORY)

        # get the photometry from the differenced image
        Utils.log('Image subtraction complete.', 'info')

        return

    @staticmethod
    def prep_ois(master, master_header):
        '''This function will prepare the files necessary for ois difference.

        :parameter master - The master image for differencing
        :parameter master_header - The header file for the master image

        :return - Nothing is returned, but the necessary text files are written and the code is compiled for differencing
        '''

        codebase = 'oisdifference.c'
        exec_name = 'a.out'

        # compile the oisdifference.c code
        shutil.copy2(codebase,  Configuration.CODE_DIFFERENCE_DIRECTORY)
        os.chdir(Configuration.CODE_DIFFERENCE_DIRECTORY)

        system_command = 'gcc ' + codebase + ' -lcfitsio -lm -o ' + exec_name

        os.system(system_command)
        os.chdir(Configuration.OUTPUT_DIRECTORY)

        # write the new master file
        fits.writeto(Configuration.CODE_DIFFERENCE_DIRECTORY + 'ref.fits', master, master_header, overwrite=True)

        # prepare the text files
        Utils.write_txt(Configuration.CODE_DIFFERENCE_DIRECTORY + 'ref.txt', 'w', 'ref' + Configuration.FILE_EXTENSION)
        Utils.write_txt(Configuration.CODE_DIFFERENCE_DIRECTORY + 'img.txt', 'w', 'img' + Configuration.FILE_EXTENSION)

        return