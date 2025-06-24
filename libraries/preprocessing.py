from config import Configuration
from libraries.utils import Utils
import os
from astropy.io import fits

class Preprocessing:

	@staticmethod
	def mk_bias(overwrite_bias, combine_type):
		''' This function will make the master bias frame using the provided image list and desired method.

		:parameter overwrite_bias - Y/N if you want to force an overwrite for the bias frame
		:parameter combine_type - Right now the combination type is mean, but it can be updated for whatever method you desire

		:return - The bias frame is returned and written to the calibration directory
		'''

		file_name = 'bias' + Configuration.FILE_EXTENSION

		if (os.path.isfile(Configuration.CALIBRATION_DIRECTORY + file_name) == 0) | (overwrite_bias == 'Y'):
			
			if combine_type == 'mean':
				
				chk_bias_list = Utils.get_file_list(Configuration.BIAS_DIRECTORY, Configuration.FILE_EXTENSION)

				images = pd.read_csv(Configuration.CALIBRATION_DIRECTORY + 'bias_list.csv', sep=',')
				image_list = images.apply(lambda x: Configuration.RAW_DIRECTORY + x.Date + x.File, axis=1).to_list()

				nfiles = len(image_list)

				if len(chk_bias_list) == 0:
					
					Utils.log('Generating a master bias frame from multiple files using a mean combination. There are ' + str(nfiles) + ' images to combine.', 'info')

					tmp_num = 0

					bias_header = fits.getheader(image_list[0])

					for kk in range(0, nfiles):

						bias_tmp = fits.getdata(image_list[kk]).astype('float')

						try:
							bias_image += bias_tmp
						except:
							bias_image = bias_tmp

						if (kk % 20) == 0:
							Utils.log(str(kk + 1) + ' files read in.' + str(nfiles - kk - 1) + ' files remain.', 'info')

						if ((kk % 100 == 0) & (kk > 0)) | (kk == nfiles - 1):
							Utils.log('Writing temporary file ' + str(tmp_num) + '. ' + str(nfiles - kk - 1) + ' files remain.', 'info')

							fits.writeto(Configuration.BIAS_DIRECTORY + str(tmp_num) + '_tmp_bias.fits', bias_image, bias_header, overwrite=True)

							tmp_num += 1

							del bias_image

						del bias_tmp

				else:
					Utils.log('Temporary bias files found. Skipping generating files. Delete them if you need to!', 'info')

				tmp_bias_list = Utils.get_file_list(Configuration.BIAS_DIRECTORY, Configuration.FILE_EXTENSION)

				for kk in range(0, len(tmp_bias_list)):

					bias_tmp = fits.getdata(Configuration.BIAS_DIRECTORY + tmp_bias_list[kk]).astype('float')

					try:
						bias_image_fin += bias_tmp
					except:
						bias_image_fin = bias_tmp

				del bias_tmp

				bias_image_mean = bias_image_fin / nfiles
				del bias_image_fin

				bias_header['BIAS_COMB'] = 'mean'
				bias_header['NUM_BIAS'] = nfiles

				fits.writeto(Configuration.CALIBRATION_DIRECTORY + file_name, bias_image_mean, bias_header, overwrite=True)

			else:
				Utils.log('Specific bias-combination method is unavailable, try again.', 'info')

		else:
			bias_image_mean = fits.getdata(Configuration.CALIBRATION_DIRECTORY + file_name, 0)

		return bias_image_mean