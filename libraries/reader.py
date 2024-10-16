from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from astropy.wcs import WCS
from tabulate import tabulate

#import warnings
#warnings.filterwarnings("ignore")

class Reader:

	@staticmethod
	def config_camera(name="PL16803"):

		params = {}

		if name == "ST8300":
			dx = 3352
			dy = 2532
			gain = 2.48
			inverse_gain = 0.403
			read_noise = 28.5

		else:
			dx = 4096
			dy = 4096
			gain = 0.72
			inverse_gain = 1.39
			read_noise = 11.4	

		params["dx"] = dx
		params["dy"] = dy
		params["gain"] = gain
		params["inverse_gain"] = inverse_gain
		params["read_noise"] = read_noise

		return params

	@staticmethod
	def read_directory(path, verbose=True):

		os.chdir(path)
		frame_list = []

		for item in glob.glob("*.fit"):
			frame_list.append(item)
		frame_list = sorted(frame_list)

		data_table = []

		for item in frame_list:
			param_list = []

			frame = fits.open(item)
			frame_header = frame[0].header

			binning = frame_header["XBINNING"]
			binning = str(binning) + "x" + str(binning)

			exptime = frame_header["EXPTIME"]
			exptime = str(exptime)

			dateobs = frame_header["DATE-OBS"]
			obstime = Time(dateobs, format="isot")

			try:
				filtertype = frame_header["FILTER"]
			except KeyError:
				filtertype = "n/a"

			imagetype = frame_header["IMAGETYP"]

			param_list.append(item)
			param_list.append(imagetype)
			param_list.append(obstime)
			param_list.append(exptime)
			param_list.append(filtertype)
			param_list.append(binning)

			data_table.append(param_list)

			if verbose == True:
				tabulated = tabulate(data_table, headers=["Item", "Type", "Time", "Exposure", "Filter", "Binning"], tablefmt="github")
				print()
				print(tabulated)
			
		return frame_list, data_table

	@staticmethod
	def read_frame(fits_frame):
		
		object_name = fits_frame[:-4]
		sigma = 3.0

		params = {}

		frame = fits.open(fits_frame)
		frame_data = frame[0].data
		frame_header = frame[0].header

		wcs = WCS(frame_header)

		dateobs = frame_header["DATE-OBS"]
		jd = frame_header["JD"]
		obs_time = Time(dateobs)

		ra = frame_header["CRVAL1"]
		dec = frame_header["CRVAL2"]
		x = frame_header["CRPIX1"]
		y = frame_header["CRPIX2"]

		exptime = frame_header["EXPTIME"]

		mean, median, std = sigma_clipped_stats(frame_data, sigma=sigma)

		params["dateobs"] = dateobs
		params["jd"] = jd
		params["ra"] = ra
		params["dec"] = dec
		params["x"] = x
		params["y"] = y
		params["exptime"] = exptime
		params["mean"] = mean
		params["median"] = median
		params["std"] = std

		return params