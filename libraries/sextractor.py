"""
SExtractor class

"""

import numpy as np
import subprocess

class Sextractor:

	@staticmethod
	def sextractor(object_frame):

		object_name = object_frame[:-4]

		Sextractor.sextractor_conv()
		Sextractor.sextractor_nnw()
		Sextractor.sextractor_param()
		Sextractor.sextractor_sex(object_name)

		subprocess.run(["source-extractor", object_frame])
		subprocess.run(["rm", "default.conv"])
		subprocess.run(["rm", "default.nnw"])

		seeing_pix_list = []
		seeing_sky_list = []
		growth_radius_list = []

		catalog = open(object_name + ".cat")

		for line in catalog:

			line = line.split()

			if line[0] == "#":
				pass

			else:
				if line[1] == "0.0000000" or line[2] == "+0.0000000":
					pass

				else:
					seeing_pix = float(line[8])
					seeing_pix_list.append(seeing_pix)

					seeing_sky = float(line[9]) * 3600
					seeing_sky_list.append(seeing_sky)

					growth_radius = float(line[7])
					growth_radius_list.append(growth_radius)

		mean_seeing_pix = np.mean(seeing_pix_list)
		mean_seeing_sky = np.mean(seeing_sky_list)
		mean_growth_radius = np.mean(growth_radius_list)

		return mean_seeing_pix, mean_seeing_sky, mean_growth_radius

	@staticmethod
	def sextractor_conv():

		conv_norm_1 = "CONV NORM\n"
		conv_norm_2 = "# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.\n"
		conv_norm_3 = "1 2 1\n"
		conv_norm_4 = "2 4 2\n"
		conv_norm_5 = "1 2 1"

		conv_norm_lines = [conv_norm_1, conv_norm_2, conv_norm_3, conv_norm_4, conv_norm_5]

		default_conv = open("default.conv", "w")
		for line in conv_norm_lines:
			default_conv.write(line)

		default_conv.close()

		return

	@staticmethod
	def sextractor_nnw():

		nnw_1 = "NNW\n"
		nnw_2 = "# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)\n"
		nnw_3 = "# inputs:  9 for profile parameters + 1 for seeing.\n"
		nnw_4 = "# outputs: ``Stellarity index'' (0.0 to 1.0)\n"
		nnw_5 = "# Seeing FWHM range: from 0.025 to 5.5'' (images must have 1.5 < FWHM < 5 pixels)\n"
		nnw_6 = "# Optimized for Moffat profiles with 2<= beta <= 4.\n\n"
		nnw_7 = " 3 10 10 1\n\n"
		nnw_8 = "-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01 -5.22453e-01  4.61342e-02  8.31957e-01  2.15505e+00  2.64769e-01\n"
		nnw_9 = " 3.03477e+00  2.69561e+00  3.16188e+00  3.34497e+00  3.51885e+00  3.65570e+00  3.74856e+00  3.84541e+00  4.22811e+00  3.27734e+00\n\n"
		nnw_10 = "-3.22480e-01 -2.12804e+00  6.50750e-01 -1.11242e+00 -1.40683e+00 -1.55944e+00 -1.84558e+00 -1.18946e-01  5.52395e-01 -4.36564e-01 -5.30052e+00\n"
		nnw_11 = " 4.62594e-01 -3.29127e+00  1.10950e+00 -6.01857e-01  1.29492e-01  1.42290e+00  2.90741e+00  2.44058e+00 -9.19118e-01  8.42851e-01 -4.69824e+00\n"
		nnw_12 = "-2.57424e+00  8.96469e-01  8.34775e-01  2.18845e+00  2.46526e+00  8.60878e-02 -6.88080e-01 -1.33623e-02  9.30403e-02  1.64942e+00 -1.01231e+00\n"
		nnw_13 = " 4.81041e+00  1.53747e+00 -1.12216e+00 -3.16008e+00 -1.67404e+00 -1.75767e+00 -1.29310e+00  5.59549e-01  8.08468e-01 -1.01592e-02 -7.54052e+00\n"
		nnw_14 = " 1.01933e+01 -2.09484e+01 -1.07426e+00  9.87912e-01  6.05210e-01 -6.04535e-02 -5.87826e-01 -7.94117e-01 -4.89190e-01 -8.12710e-02 -2.07067e+01\n"
		nnw_15 = "-5.31793e+00  7.94240e+00 -4.64165e+00 -4.37436e+00 -1.55417e+00  7.54368e-01  1.09608e+00  1.45967e+00  1.62946e+00 -1.01301e+00  1.13514e-01\n"
		nnw_16 = " 2.20336e-01  1.70056e+00 -5.20105e-01 -4.28330e-01  1.57258e-03 -3.36502e-01 -8.18568e-02 -7.16163e+00  8.23195e+00 -1.71561e-02 -1.13749e+01\n"
		nnw_17 = " 3.75075e+00  7.25399e+00 -1.75325e+00 -2.68814e+00 -3.71128e+00 -4.62933e+00 -2.13747e+00 -1.89186e-01  1.29122e+00 -7.49380e-01  6.71712e-01\n"
		nnw_18 = "-8.41923e-01  4.64997e+00  5.65808e-01 -3.08277e-01 -1.01687e+00  1.73127e-01 -8.92130e-01  1.89044e+00 -2.75543e-01 -7.72828e-01  5.36745e-01\n"
		nnw_19 = "-3.65598e+00  7.56997e+00 -3.76373e+00 -1.74542e+00 -1.37540e-01 -5.55400e-01 -1.59195e-01  1.27910e-01  1.91906e+00  1.42119e+00 -4.35502e+00\n\n"
		nnw_20 = "-1.70059e+00 -3.65695e+00  1.22367e+00 -5.74367e-01 -3.29571e+00  2.46316e+00  5.22353e+00  2.42038e+00  1.22919e+00 -9.22250e-01 -2.32028e+00\n\n\n"
		nnw_21 = " 0.00000e+00\n"
		nnw_22 = " 1.00000e+00"

		default_nnw_lines = [nnw_1, nnw_2, nnw_3, nnw_4, nnw_5, nnw_6, nnw_7, nnw_8, nnw_9, nnw_10, nnw_11, nnw_12, nnw_13, nnw_14, nnw_15, nnw_16, nnw_17, nnw_18, nnw_19, nnw_20, nnw_21, nnw_22]

		default_nnw = open("default.nnw", "w")
		for line in default_nnw_lines:
			default_nnw.write(line)
		default_nnw.close()

		return

	@staticmethod
	def sextractor_param():

		param_number = "NUMBER"
		param_alphapeak_j2000 = "ALPHAPEAK_J2000"
		param_deltapeak_j2000 = "DELTAPEAK_J2000"
		param_xpeak_image = "XPEAK_IMAGE"
		param_ypeak_image = "YPEAK_IMAGE"
		param_flux_growth = "FLUX_GROWTH"
		param_fluxerr_best = "FLUXERR_BEST"
		param_flux_growthstep = "FLUX_GROWTHSTEP"
		param_fwhm_image = "FWHM_IMAGE"
		param_fwhm_world = "FWHM_WORLD"
		
		default_parameters = [param_number, param_alphapeak_j2000, param_deltapeak_j2000, param_xpeak_image, param_ypeak_image, param_flux_growth, param_fluxerr_best, param_flux_growthstep, param_fwhm_image, param_fwhm_world]

		default_param = open("default.param", "w")
		for param in default_parameters:
			default_param.write(param + "\n")
		default_param.close()

		return

	@staticmethod
	def sextractor_sex(object_name, inverse_gain=1.39, pixel_scale=0.63):

		# Default configuration file for SExtractor 2.12.4
		# EB 2010-10-10

		catalog_name = ["CATALOG_NAME", object_name + ".cat"]
		catalog_type = ["CATALOG_TYPE", "ASCII_HEAD"]
		parameters_name = ["PARAMETERS_NAME", "default.param"]

		detect_type = ["DETECT_TYPE", "CCD"] # CCD (linear) or PHOTO (with gamma correction)
		detect_minarea = ["DETECT_MINAREA", "3"] # min. # of pixels above threshold
		detect_thresh = ["DETECT_THRESH", "5.0"] # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
		analysis_thresh = ["ANALYSIS_THRESH", "5.0"] # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
		detection_filter = ["FILTER", "Y"] # apply filter for detection (Y or N)?
		detection_filter_name = ["FILTER_NAME", "default.conv"] # name of the file containing the filter
		deblend_nthresh = ["DEBLEND_NTHRESH", "32"] # Number of deblending sub-thresholds
		deblend_mincont = ["DEBLEND_MINCONT", "0.005"] # Minimum contrast parameter for deblending
		clean = ["CLEAN", "Y"] # Clean spurious detections? (Y or N)?
		clean_param = ["CLEAN_PARAM", "1.0"] # Cleaning efficiency
		
		weight_type = ["WEIGHT_TYPE", "NONE"] # type of WEIGHTing: NONE, BACKGROUND, MAP_RMS, MAP_VAR or MAP_WEIGHT
		weight_image = ["WEIGHT_IMAGE", "weight.fits"] # weight-map filename
		
		flag_image = ["FLAG_IMAGE", "flag.fits"] # filename for an input FLAG-image
		flag_type = ["FLAG_TYPE", "OR"] # flag pixel combination: OR, AND, MIN, MAX or MOST
		
		phot_apertures = ["PHOT_APERTURES", "10"] # MAG_APER aperture diameter(s) in pixels
		phot_autoparams = ["PHOT_AUTOPARAMS", "2.5", "3.5"] # MAG_AUTO parameters: <Kron_fact>,<min_radius>
		phot_petroparams = ["PHOT_PETROPARAMS", "2.0", "3.5"] # MAG_PETRO parameters: <Petrosian_fact>,<min_radius>
		phot_autoapers = ["PHOT_AUTOAPERS", "0.0", "0.0"] # <estimation>,<measurement> minimum apertures for MAG_AUTO and MAG_PETRO
		satur_level = ["SATUR_LEVEL", "65535.0"] # level (in ADUs) at which arises saturation
		satur_key = ["SATUR_KEY", "SATURATE"] # keyword for saturation level (in ADUs)
		mag_zeropoint = ["MAG_ZEROPOINT", "0.0"] # magnitude zero-point
		mag_gamma = ["MAG_GAMMA", "4.0"] # gamma of emulsion (for photographic scans)
		gain = ["GAIN", str(inverse_gain)] # keyword for detector gain in e-/ADU
		gain_key = ["GAIN_KEY", "GAIN"] # keyword for detector gain in e-/ADU
		pixel_scale = ["PIXEL_SCALE", str(pixel_scale)] # size of pixel in arcsec (0=use FITS WCS info)
		
		seeing_fwhm = ["SEEING_FWHM", "1.2"] # stellar FWHM in arcsec
		starnnw_name = ["STARNNW_NAME", "default.nnw"] # Neural-Network_Weight table filename
		
		back_type = ["BACK_TYPE", "AUTO"] # AUTO or MANUAL
		back_value = ["BACK_VALUE", "0.0"] # Default background value in MANUAL mode
		back_size = ["BACK_SIZE", "64"] # Background mesh: <size> or <width>,<height>
		back_filtersize = ["BACK_FILTERSIZE", "3"] # Background filter: <size> or <width>,<height>
		
		checkimage_type = ["CHECKIMAGE_TYPE", "NONE"] # can be NONE, BACKGROUND, BACKGROUND_RMS, MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND, FILTERED, OBJECTS, -OBJECTS, SEGMENTATION, or APERTURES
		checkimage_name = ["CHECKIMAGE_NAME", "check.fits"] # Filename for the check-image
		
		memory_objstack = ["MEMORY_OBJSTACK", "3000"] # number of objects in stack
		memory_pixstack = ["MEMORY_PIXSTACK", "300000"] # number of pixels in stack
		memory_bufsize = ["MEMORY_BUFSIZE", "1024"] # number of lines in buffer
		
		assoc_name = ["ASSOC_NAME", "sky.list"] # name of the ASCII file to ASSOCiate
		assoc_data = ["ASSOC_DATA", "2", "3", "4"] # columns of the data to replicate (0=all)
		assoc_params = ["ASSOC_PARAMS", "2", "3", "4"] # columns of xpos,ypos[,mag]
		assoc_radius = ["ASSOC_RADIUS", "2.0"] # cross-matching radius (pixels)
		assoc_type = ["ASSOC_TYPE", "NEAREST"] # ASSOCiation method: FIRST, NEAREST, MEAN, MAG_MEAN, SUM, MAG_SUM, MIN or MAX
		assocselec_type = ["ASSOCSELEC_TYPE", "MATCHED"] # ASSOC selection type: ALL, MATCHED or -MATCHED
		
		verbose_type = ["VERBOSE_TYPE", "NORMAL"] # can be QUIET, NORMAL or FULL
		header_suffix = ["HEADER_SUFFIX", ".head"] # Filename extension for additional headers
		write_xml = ["WRITE_XML", "N"]# Write XML file (Y/N)?
		xml_name = ["XML_NAME", "sex.xml"] # Filename for XML output
		xsl_url = ["XSL_URL", "file:///usr/local/share/sextractor/sextractor.xsl"] # Filename for XSL style-sheet

		config_parameters = [catalog_name, catalog_type, parameters_name, detect_type, detect_minarea, detect_thresh,
								analysis_thresh, detection_filter, detection_filter_name, deblend_nthresh, deblend_mincont, 
								clean, clean_param, weight_type, weight_image, flag_image, flag_type, phot_apertures, 
								phot_autoparams, phot_petroparams, phot_autoapers, satur_level, satur_key, mag_zeropoint,
								mag_gamma, gain, gain_key, pixel_scale, seeing_fwhm, starnnw_name, back_type, back_value,
								back_size, back_filtersize, checkimage_type, checkimage_name, memory_objstack, memory_pixstack,
								memory_bufsize, assoc_name, assoc_data, assoc_params, assoc_radius, assoc_type, assocselec_type,
								verbose_type, header_suffix, write_xml, xml_name, xsl_url]

		default_sex = open("default.sex", "w")
		for param in config_parameters:
			for item in param:
				default_sex.write(item + " ")
			default_sex.write("\n")
		default_sex.close()

		return