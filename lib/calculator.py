from config import Configuration
from lib.utility import Utility

from astropy.coordinates import AltAz, SkyCoord
from scipy.optimize import curve_fit
import numpy as np

class Calculator:

	@staticmethod
	def airmass(h_deg, method=Configuration.AIRMASS_METHOD):
		'''This function calculates the airmass of a given altitude (in degrees).

		:parameter h_deg - The altitude (in degrees)
		:parameter method - A string for the calculation method

		:return x - The airmass
		'''

		z_deg = 90. - h_deg
		z_rad = np.deg2rad(z_deg)

		if method == 'ky1989':
			x = 1 / (np.cos(z_rad) + (0.50572*((6.07995 + h_deg)**-1.6364)))
		else:
			x = 1 / np.cos(z_rad)

		return x

	@staticmethod
	def angular_distance(ra1, de1, ra2, de2):
		'''This function determines the angular distance between two sky coordiantes.

		:parameter ra1 [float, deg] - The right ascension of point 1
		:parameter de1 [float, deg] - The declination of point 1
		:parameter ra2 [float, deg] - The right ascension of point 2
		:parameter de2 [float, deg] - The declination of point 2

		:return ang_dist [float, deg] - The angular distance between the two points
		'''

		# convert to radians
		ra1_rad = np.deg2rad(ra1)
		de1_rad = np.deg2rad(de1)
		ra2_rad = np.deg2rad(ra2)
		de2_rad = np.deg2rad(de2)

		# the angular distance in radians
		ang_dist_rad = np.arccos((np.sin(de1_rad) * np.sin(de2_rad)) + (np.cos(de1_rad) * np.cos(de2_rad) * np.cos(ra1_rad - ra2_rad)))

		# the angular distance in degrees
		ang_dist = np.rad2deg(ang_dist_rad)

		return ang_dist

	@staticmethod
	def sigma_clip(xlist, ylist, sigma):

		xarr = np.asarray(xlist)
		yarr = np.asarray(ylist)

		xm = np.mean(xarr)
		ym = np.mean(yarr)

		xs = np.std(xarr)
		ys = np.std(yarr)

		new_xlist = []
		new_ylist = []

		for i in range(len(xlist)):
			if ylist[i] < (ym - (sigma * ys)):
				pass
			else:
				new_xlist.append(xlist[i])
				new_ylist.append(ylist[i])

		return new_xlist, new_ylist

	@staticmethod
	def sky_to_altitude(location, time, ra, dec):

		aa = AltAz(location=location, obstime=time)

		coord = SkyCoord(str(ra), str(dec), unit='deg')
		coord_transform = coord.transform_to(aa)

		altitude = coord_transform.alt.degree

		return altitude

	@staticmethod
	def unweighted_fit(xlist, ylist):

		xarr = np.asarray(xlist)
		yarr = np.asarray(ylist)

		def f(x, m, b):
			y = (m*x) + b
			return y

		popt, pcov = curve_fit(f, xarr, yarr)
		yfit = f(xarr, *popt)

		m = popt[0]
		dm = np.sqrt(pcov[0][0])

		b = popt[1]
		db = np.sqrt(pcov[1][1])

		return yfit, m, dm, b, db