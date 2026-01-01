from config import Configuration
from libraries.preprocessing import Preprocessing
from libraries.utils import Utils
from astropy.io import fits
import numpy as np
import os
import sys
import time

class Clean:

    @staticmethod
    def clean_images():
        '''This is the main script to clean multiple images.

        return - no value is returned, the values images from in_path are cleaned and deposited in out_path
        '''

        # start clock
        st = time.time()

        # get file list for all dates field was observed
        Utils.log('Getting file list.', 'info')
        files, date_dirs = Utils.get_all_files_per_field(Configuration.RAW_DIRECTORY, Configuration.FIELD, 'raw', Configuration.FILE_EXTENSION)

        # make output directories for clean/difference/flux files
        output_dirs = []

        for dte in date_dirs:
            output_dirs.append(Configuration.DATA_DIRECTORY + 'clean/' + dte)
            output_dirs.append(Configuration.DATA_DIRECTORY + 'clean/' + dte + '/' + Configuration.FIELD)
            output_dirs.append(Configuration.DATA_DIRECTORY + 'diff/' + dte)
            output_dirs.append(Configuration.DATA_DIRECTORY + 'diff/' + dte + '/' + Configuration.FIELD)
            output_dirs.append(Configuration.DATA_DIRECTORY + 'flux/' + dte)
            output_dirs.append(Configuration.DATA_DIRECTORY + 'flux/' + dte + '/' + Configuration.FIELD)
            output_dirs.append(Configuration.DATA_DIRECTORY + 'review/' + dte)
            output_dirs.append(Configuration.DATA_DIRECTORY + 'review/' + dte + '/' + Configuration.FIELD)

        Utils.create_directories(output_dirs)

        # break if there are no files
        if len(files) == 0:
            Utils.log('No .fits files found for ' + Configuration.FIELD + '!' + ' Breaking...', 'debug')
            sys.exit()

        Utils.log('Starting to clean ' + str(len(files)) + ' images.', 'info')
        
        for idx, file in enumerate(files):
            # make a new name for the file based on which actions are taken
            file_name = Preprocessing.mk_nme(file, 'N', Configuration.SUBTRACT_BIAS, Configuration.SUBTRACT_DARK, Configuration.DIVIDE_FLAT, Configuration.CLIP_IMAGE, Configuration.SUBTRACT_SKY, Configuration.PLATE_SOLVE)

            # check if file exists
            if os.path.isfile(file_name):
                Utils.log('Image ' + file_name + ' already exists. Skipping.', 'info')
            else:
                # clean the image
                clean_img, header = Clean.clean_img(file, file_name)

            rmdr = len(files) - idx - 1
            Utils.log(str(rmdr) + ' images remain to be cleaned.',  'info')

        # stop clock
        fn = time.time()
        dt = np.around((fn - st), decimals=2)
        Utils.log('Total image cleaning complete in ' + str(dt) + ' s.', 'info')

        return

    @staticmethod
    def clean_img(file, file_name=None):
        '''This function is the primary script to clean an image.

        :parameter file - The file path of the unprocessed image
        :parameter file_name - The file name of the unprocessed image

        :return img - 
        :return header - 
        '''

        Utils.log('Now cleaning ' + file + '.', 'info')

        # load image
        img, header = fits.getdata(file, header=True)

        # remove bias and dark
        if (Configuration.SUBTRACT_BIAS == 'Y') and (Configuration.SUBTRACT_DARK == 'Y'):
            st = time.time()
            bias, dark = Preprocessing.mk_combined_bias_and_dark(image_overwrite='N')
            img, header = Preprocessing.subtract_scaled_bias_dark(img, header)
            fn = time.time()
            dt1 = np.around((fn - st), decimals=2)
            Utils.log('Bias and dark subtracted in ' + str(dt1) + ' s.', 'info')
        else:
            Utils.log('Skipping bias and dark subtraction.', 'info')

        # flat divide
        if Configuration.DIVIDE_FLAT == 'Y':
            st = time.time()
            flat = Preprocessing.mk_flat(Configuration.FLAT_EXP, Configuration.DARK_EXP)
            img, header = Preprocessing.flat_divide(img, header)
            fn = time.time()
            dt2 = np.around((fn - st), decimals=2)
            Utils.log('Image flattened in ' + str(dt2) + ' s.', 'info')
        else:
            Utils.log('Skipping image flattening.', 'info')

        # clip image
        if Configuration.CLIP_IMAGE == 'Y':
            st = time.time()
            img, header = Preprocessing.clip_image(img, header)
            fn = time.time()
            dt3 = np.around((fn - st), decimals=2)
            Utils.log('Overscan removed in ' + str(dt3) + ' s.', 'info')
        else:
            Utils.log('Skipping overscan removal.', 'info')

        # sky subtract
        if Configuration.SUBTRACT_SKY == 'Y':
            st = time.time()
            img, header = Preprocessing.sky_subtract(img, header, Configuration.WRITE_SKY)
            fn = time.time()
            dt4 = np.around((fn - st), decimals=2)
            Utils.log('Sky subtracted in ' + str(dt4) + ' s.', 'info')
        else:
            Utils.log('Skipping sky subtraction.', 'info')

        # plate solve
        if Configuration.PLATE_SOLVE == 'Y':
            st = time.time()
            img, header = Preprocessing.correct_header(img, header, file_name)
            fn = time.time()
            dt5 = np.around((fn - st), decimals=2)
            Utils.log('Image plate solved in ' + str(dt5) + ' s.', 'info')

        Utils.log('Cleaning finished on image ' + file + '.', 'info')

        dt = np.around((dt1 + dt2 + dt3 + dt4 + dt5), decimals=2)
        Utils.log('Total time = ' + str(dt) + ' s.', 'info')

        return img, header