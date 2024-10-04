from config import Configuration
from libraries.utils import Utils

from astropy.coordinates import AltAz, SkyCoord
from scipy.optimize import curve_fit
import math
import numpy as np

class Calculator:

	@staticmethod
	def angular_distance(ra1, dec1, ra2, dec2):
		''' This function will determine the angular distance between two points.

		:parameter ra1 - The right ascension of point 1
		:parameter dec1 - The declination of point 1
		:parameter ra2 - The right ascension of point 2
		:parameter dec2 - The declination of point 2

		:return ang_dist_deg - The angular distance of the two points in degrees

		'''

		# convert to radians
		ra1_rad = np.deg2rad(ra1)
		dec1_rad = np.deg2rad(dec1)
		ra2_rad = np.deg2rad(ra2)
		dec2_rad = np.deg2rad(dec2)

		# the angular distance in radians
		ang_dist_rad = np.arccos((np.sin(dec1_rad) * np.sin(dec2_rad)) + (np.cos(dec1_rad) * np.cos(dec2_rad) * np.cos(ra1_rad - ra2_rad)))

		# the angular distance in degrees
		ang_dist_deg = np.rad2deg(ang_dist_rad)

		return ang_dist_deg

	@staticmethod
	def calculate_altitude(location, time, ra, dec):

		aa = AltAz(location=location, obstime=time)
		coord = SkyCoord(str(ra), str(dec), unit="deg")
		coord_transform = coord.transform_to(aa)

		altitude = coord_transform.alt.degree

		return altitude

	@staticmethod
	def calculate_airmass(altitude, method="ky1998"):

		zenith_distance = 90 - altitude

		if method == "p2002":
			airmass = 1 / math.sin(altitude + (244/(47*altitude**1.1))) # Pickering (2002) method

		else:
			airmass = 1 / math.cos(math.radians(zenith_distance)) + 0.50572*((6.07995 + 90 - zenith_distance)**-1.6364) # Kasten and Young (1994) method

		return airmass

	@staticmethod
	def sigma_clip(x_list, y_list, sigma):

		x_array = np.asarray(x_list)
		y_array = np.asarray(y_list)

		mean_x = np.mean(x_array)
		mean_y = np.mean(y_array)

		std_x = np.std(x_array)
		std_y = np.std(y_array)

		new_x_list = []
		new_y_list = []

		for i in range(len(x_list)):

			if y_list[i] < (mean_y - (sigma * std_y)):

				pass

			else:

				new_x_list.append(x_list[i])
				new_y_list.append(y_list[i])

		return new_x_list, new_y_list

	@staticmethod
	def unweighted_fit(x_list, y_list):

		x_array = np.asarray(x_list)
		y_array = np.asarray(y_list)

		def f(x, m, b):
			y = (m*x) + b
			return y

		popt, pcov = curve_fit(f, x_array, y_array)

		yfit = f(x_array, *popt)

		slope = popt[0]
		delta_slope = math.sqrt(pcov[0][0])

		intercept = popt[1]
		delta_intercept = math.sqrt(pcov[1][1])

		return yfit, slope, intercept, delta_slope, delta_intercept