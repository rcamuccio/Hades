'''
Λογιστής - The Calculator
'''

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, SkyCoord
from scipy.optimize import curve_fit

from conf import Configuration

class Logistes:

	@staticmethod
	def get_airmass(h_deg, method=Configuration.AIRMASS_METHOD):
		'''This function calculates the airmass of a given altitude.

		:parameter h_deg [float] - the altitude [deg]
		:parameter method [str] - the method of calculation

		:return x [float] - the airmass
		'''

		# calculate the zenith distance
		z_deg = 90. - h_deg

		# convert the zenith distance units to radians
		z_rad = np.deg2rad(z_deg)

		# calculate the airmass based on the method
		if method == 'ky1989':
			x = 1 / (np.cos(z_rad) + (0.50572*((6.07995 + h_deg)**-1.6364)))
		else:
			x = 1 / np.cos(z_rad)

		return x

	@staticmethod
	def get_altitude(location, time, ra, dec):
		'''This function calculates the altitude of a sky position for a given location and time.

		:parameter location [astropy.coordinates.EarthLocation] - the location (lat, lon, elev) of the observation [deg, deg, m]
		:parameter time [astropy.time.Time] - the time of the observation
		:parameter ra [float] - the right ascension of the sky position [deg]
		:parameter dec [float] - the declination of the sky position [deg]

		:return alt [float] - the altitude of the sky position [deg]
		'''

		# generate an AltAz object for the observer
		alt_az = AltAz(location=location, obstime=time)

		# generate a SkyCoord object for the sky position
		coord_sky = SkyCoord(str(ra), str(dec), unit='deg')

		# transform the sky position to a local position
		coord_loc = coord.transform_to(alt_az)

		# get the altitude of the local position
		alt = coord_transform.alt.degree

		return alt

	@staticmethod
	def get_angular_distance(ra_1_deg, dec_1_deg, ra_2_deg, dec_2_deg):
		'''This function calculates the angular distance between two sky coordinates.

		:parameter ra_1_deg [float] - the right ascension of point #1 [deg]
		:parameter dec_1_deg [float] - the declination of point #1 [deg]
		:parameter ra_2_deg [float] - the right ascension of point #2 [deg]
		:parameter dec_2_deg [float] - the declination of point #2 [deg]

		:return ang_dist_deg [float] - the angular distance between the two points [deg]
		'''

		# convert the inputs to radians
		ra_1_rad = np.deg2rad(ra_1_deg)
		dec_1_rad = np.deg2rad(dec_1_deg)
		ra_2_rad = np.deg2rad(ra_2_deg)
		dec_2_rad = np.deg2rad(dec_2_deg)

		# calculate the angular distance in radians
		ang_dist_rad = np.arccos((np.sin(dec_1_rad) * np.sin(dec_2_rad)) + (np.cos(dec_1_rad) * np.cos(dec_2_rad) * np.cos(ra_1_rad - ra_2_rad)))

		# convert the angular distance units to degrees
		ang_dist_deg = np.rad2deg(ang_dist_rad)

		return ang_dist_deg

	@staticmethod
	def get_galactic_center():

		gc = SkyCoord(l=0*u.deg, b=0*u.deg, frame='galactic').icrs

		gc_ra = gc.ra.degree
		gc_dec = gc.dec.degree

		gc_ra_shift = np.remainder(gc_ra + 180, 360) - 180

		gc_ra_rad = np.radians(gc_ra_shift)
		gc_dec_rad = np.radians(gc_dec)

		return gc_ra_rad, gc_dec_rad

	@staticmethod
	def get_galactic_plane(contour_level):

		# generate the contour
		l = np.linspace(0, 360, 2000) * u.deg
		if contour_level == 0.:
			b = np.zeros_like(l)
		else:
			b = np.full_like(l.value, contour_level) * u.deg
		gal_plane = SkyCoord(l=l, b=b, frame='galactic')

		# convert the plane to equatorial coordinates
		eq_plane = gal_plane.icrs
		eq_plane_ra = eq_plane.ra.degree
		eq_plane_dec = eq_plane.dec.degree

		# convert the right ascension values into range [-180, +180]
		ra_shift = np.remainder(eq_plane_ra + 180, 360) - 180
		ra_rad = np.radians(ra_shift)
		dec_rad = np.radians(eq_plane_dec)

		# fix the discontinuity
		jumps = np.where(np.abs(np.diff(ra_rad)) > np.pi/2)[0]
		ra_plot = ra_rad.copy()
		dec_plot = dec_rad.copy()
		for idx in reversed(jumps):
			ra_plot = np.insert(ra_plot, idx + 1, np.nan)
			dec_plot = np.insert(dec_plot, idx + 1, np.nan)

		return ra_plot, dec_plot

	@staticmethod
	def make_unweighted_fit(x_list, y_list):
		'''This function calculates the unweighted linear fit (y=mx+b) of a set of data lists (x, y).

		:parameter x_list [list] - the x input array
		:parameter y_list [list] - the y input array

		:return y_fit [numpy.ndarray] - the linear fit
		:return m [float] - the best-fit slope
		:return dm [float] - the slope error
		:return b [float] - the best-fit intercept
		:return db [float] - the intercept error
		'''

		# convert the lists to numpy arrays
		x_array = np.asarray(x_list)
		y_array = np.asarray(y_list)

		# generate a model linear function for the fitting
		def f(x, m, b):
			y = (m*x) + b
			return y

		# perform the linear fit
		popt, pcov = curve_fit(f, x_array, y_array)

		# generate a linear function from the fit
		y_fit = f(x_array, *popt)

		# get the slope and its error
		m = popt[0]
		dm = np.sqrt(pcov[0][0])

		# get the intercept and its error
		b = popt[1]
		db = np.sqrt(pcov[1][1])

		return y_fit, m, dm, b, db		