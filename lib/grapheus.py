'''
Γραφεύς - The Painter
'''

import astropy.units as u
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.visualization import ZScaleInterval

from conf import Configuration
from lib.logistes import Logistes
from lib.skopeus import Skopeus

import matplotlib as mpl
#mpl.rc('font', family=Configuration.FONT_NAME, size=Configuration.FONT_SIZE)

class Grapheus:

	@staticmethod
	def draw_airmass(time_list, airmass_list, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws an air mass diagram.

		:p time_list [list][float] - the time list [d]
		:p airmass_list [list][float] - the air mass list
		:p save_path [str] - the diagram save path
		:p save_figure [str] - a toggle to save or display the diagram

		:r - nothing
		'''

		# draw the diagram
		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.plot(time_list, airmass_list, s=1, color='black')
		plt.xlabel('Time [JD]')
		plt.ylabel('Air mass')

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def draw_color_magnitude(color_list, delta_list, linear_fit, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a color-magnitude diagram.

		:p color_list [list][float] - the color list [mag]
		:p delta_list [list][float] - the delta magnitude list [mag]
		:p linear_fit [list][float] - the linear fit [mag/mag]
		:p save_path [str] - the diagram save path
		:p save_figure [str] - a toggle to save or display the diagram

		:r - nothing
		'''

		# draw the diagram
		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.scatter(color_list, delta_list, s=1, color='gray')
		plt.plot(color_list, linear_fit, ls='--', lw=0.75, color='black')
		plt.xlabel('Color index [mag]')
		plt.ylabel('Delta magnitude [mag]')

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def draw_extinction(airmass_list, mag_list, linear_fit, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws an extinction diagram.

		:p airmass_list [list][float] - the air mass list
		:p mag_list [list][float] - the magnitude list [mag]
		:p linear_fit [list][float] - the linear fit [mag/airmass]
		:p save_path [str] - the diagram save path
		:p save_figure [str] - a toggle to save or display the diagram

		:r - nothing
		'''

		# draw the diagram
		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.scatter(airmass_list, mag_list, s=1, color='gray')
		plt.plot(airmass_list, linear_fit, ls='--', lw=0.75, color='black')
		plt.xlabel('Air mass')
		plt.ylabel('Instrumental magnitude [mag]')

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def draw_field_diagram(frame_data, apertures, boxes, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a field with annotations.

		:p frame_data [numpy.ndarray] - the image data array
		:p apertures [photutils.aperture.circle.CircularAperture] - the set of circular apertures
		:p boxes [dict][photutils.aperture.rectangle.RectangularAperture] - the set of mask boxes
		:p save_path [str] - the diagram save path
		:p save_figure [str] - a toggle to save or display the diagram

		:return - nothing
		'''

		# calculate the display limits
		interval = ZScaleInterval()
		vmin, vmax = interval.get_limits(frame_data)

		# draw the diagram
		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.imshow(frame_data, cmap=Configuration.COLORMAP, origin='lower', vmin=vmin, vmax=vmax)
		plt.xlabel('x [px]')
		plt.ylabel('y [px]')
		plt.colorbar()

		# draw the apertures
		if apertures != None:
			apertures.plot(color='lime', lw=0.5, alpha=0.5)

		# draw the mask boxes
		if boxes != None:
			eml_box = boxes['eml_box']
			emt_box = boxes['emt_box']
			emr_box = boxes['emr_box']
			emb_box = boxes['emb_box']

			eml_box.plot(color='cyan', ls='dashed')
			emt_box.plot(color='lime', ls='dashed')
			emr_box.plot(color='yellow', ls='dashed')
			emb_box.plot(color='orange', ls='dashed')

		# save or display the figure
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def draw_growth_radius(time_list, radius_list, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a growth radius diagram.

		:p time_list [list][float] - the time list [JD]
		:p radius_list [list][float] - the radius list [px]
		:p save_path [str] - the diagram save path
		:p save_figure [str] - a toggle to save or display the diagram

		:return - nothing
		'''

		# draw the diagram
		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.plot(time_list, radius_list, s=1, color='black')
		plt.xlabel('Time [JD]')
		plt.ylabel('Growth radius [px]')

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def draw_histogram_field(frame_data, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a histogram of data array pixel values.

		:p frame_data [numpy.ndarray] - the image data array
		:p save_path [str] - the diagram save path
		:p save_figure [str] - a toggle to save or display the diagram

		:return - nothing
		'''

		# set the histogram limits
		data_mn, data_md, data_sd = sigma_clipped_stats(frame_data, sigma=Configuration.SIGMA_BACKGROUND)
		xmin = data_md - Configuration.HISTOGRAM_LIMIT * data_sd
		xmax = data_md + Configuration.HISTOGRAM_LIMIT * data_sd

		# draw the diagram
		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.hist(frame_data.flatten(), bins=Configuration.HISTOGRAM_BINS, range=(xmin, xmax), histtype=Configuration.HISTOGRAM_TYPE)
		plt.yscale(Configuration.HISTOGRAM_SCALE)
		plt.xlabel('Pixel value [ADU]')
		plt.ylabel('Count')

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def draw_histogram_flux(flux_list, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a histogram of stellar fluxes.

		:p flux_list [list][float] - the flux list [ADU]
		:p save_path [str] - the diagram save path
		:p save_figure [str] - a toggle to save or display the diagram

		:return - nothing
		'''

		# draw the diagram
		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.hist(flux_list, bins=Configuration.HISTOGRAM_BINS, range=(0, 2.5e5), histtype=Configuration.HISTOGRAM_TYPE)
		plt.yscale(Configuration.HISTOGRAM_SCALE)
		plt.xlabel('Flux [ADU]')
		plt.ylabel('Count')

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def draw_histogram_gaia(r_mag_list, g_mag_list, b_mag_list, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a histogram of stellar Gaia magnitudes.

		:p r_mag_list [list][float] - the Gaia R magnitude list [mag]
		:p g_mag_list [list][float] - the Gaia G magnitude list [mag]
		:p b_mag_list [list][float] - the Gaia B magnitude list [mag]
		:p save_path [str] - the diagram save path
		:p save_figure [str] - a toggle to save or display the diagram

		:return - nothing
		'''

		# dras the diagram
		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.hist(r_mag_list, bins=Configuration.HISTOGRAM_BINS, histtype=Configuration.HISTOGRAM_TYPE, color='red', label='Gaia R')
		plt.hist(g_mag_list, bins=Configuration.HISTOGRAM_BINS, histtype=Configuration.HISTOGRAM_TYPE, color='green', label='Gaia G')
		plt.hist(b_mag_list, bins=Configuration.HISTOGRAM_BINS, histtype=Configuration.HISTOGRAM_TYPE, color='blue', label='Gaia B')
		plt.yscale(Configuration.HISTOGRAM_SCALE)
		plt.xlabel('Gaia magnitude [mag]')
		plt.ylabel('Count')
		plt.legend()

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def draw_histogram_magnitude(magnitude_list, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a histogram of stellar fluxes.

		:p magnitude_list [list][float] - the magnitude list [mag]
		:p save_path [str] - the diagram save path
		:p save_figure [str] - a toggle to save or display the diagram

		:return - nothing
		'''

		# draw the diagram
		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.hist(magnitude_list, bins=Configuration.HISTOGRAM_BINS, range=(-12, 4), histtype=Configuration.HISTOGRAM_TYPE)
		plt.yscale(Configuration.HISTOGRAM_SCALE)
		plt.xlabel('Instrumental magnitude [mag]')
		plt.ylabel('Count')

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def draw_lightcurve(time_list, mag_list, err_list, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a light curve diagram.

		:p time_list [list][float] - the time list [JD]
		:p mag_list [list][float] - the magnitude list [mag]
		:p err_list [list][float] - the magnitude error list [mag]
		:p save_path [str] - the diagram save path
		:p save_figure [str] - a toggle to save or display the diagram

		:return - nothing
		'''

		# draw the diagram
		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.errorbar(time_list, mag_list, yerr=err_list, fmt='o', linewidth=0.5, markersize=0.5, capsize=2, capthick=0.5)
		plt.xlabel('Time [JD]')
		plt.ylabel('Instrumental magnitude [mag]')

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def draw_photometric_precision(mag_list, rms_list, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a photometric precision diagram.

		:p mag_list [list][float] - the magnitude list [mag]
		:p rms_list [list][float] - the rms list [mag]
		:p save_path [str] - the diagram save path
		:p save_figure [str] - a toggle to save or display the diagram

		:return - nothing
		'''

		# draw the diagram
		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.scatter(mag_list, rms_list, s=1, color='gray')
		plt.yscale('log')
		plt.xlabel('Magnitude [mag]')
		plt.ylabel('Root mean square [mag]')

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def draw_seeing(time_list, seeing_list, type, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a seeing diagram.

		:p time_list [list][float] - the time list [JD]
		:p seeing_list [list][float] - the seeing list [px or arcsec]
		:p type [str] - a toggle for either 'physical' or 'sky' units
		:p save_path [str] - the diagram save path
		:p save_figure [str] - a toggle to save or display the diagram

		:return - nothing
		'''

		# draw the diagram
		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.plot(time_list, seeing_list, s=1, color='black')
		plt.xlabel('Time [JD]')
		if type == 'physical':
			plt.ylabel('Seeing [px]')
		elif type == 'sky':
			plt.ylabel('Seeing [arcsec]')

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def draw_default_skymap(ra_deg, dec_deg, save_path, save_figure):

		# set up the figure
		fig = plt.figure(figsize=Configuration.FIGURE_SIZE)
		ax = fig.add_subplot(111, projection=Configuration.PROJECTION)

		# draw the coordinates
		if (ra_deg != None) and (dec_deg != None):
			# ensure the input coordinates are array-like
			ra_deg = np.asarray(ra_deg)
			dec_deg = np.asarray(dec_deg)

			# convert the right ascension values into range [-180, +180]
			ra_shift = np.remainder(ra_deg + 180, 360) - 180

			# convert the coordinates to radians
			ra_rad = np.radians(ra_shift)
			dec_rad = np.radians(dec_deg)

			# draw the objects
			ax.scatter(ra_rad, dec_rad, c='black', s=3)

		# draw the galactic plane
		if Configuration.PLOT_GALACTIC_PLANE:
			ra_plot, dec_plot = Logistes.get_galactic_plane(contour_level=0.)
			ax.plot(ra_plot, dec_plot, ls='--', lw=0.75, c='black', label='galactic plane')

		# draw the galactic contours
		if Configuration.PLOT_GALACTIC_CONTOUR:
			ra_neg_plot, dec_neg_plot = Logistes.get_galactic_plane(contour_level=-15)
			ra_pos_plot, dec_pos_plot = Logistes.get_galactic_plane(contour_level=15.)
			ax.plot(ra_neg_plot, dec_neg_plot, ls='--', lw=0.5, c='black')
			ax.plot(ra_pos_plot, dec_pos_plot, ls='--', lw=0.5, c='black')

		# draw the galactic center
		if Configuration.PLOT_GALACTIC_CENTER:
			gc_ra_rad, gc_dec_rad = Logistes.get_galactic_center()
			ax.scatter(gc_ra_rad, gc_dec_rad, marker='x', s=100, color='black', label='galactic center')

		# draw the meridians
		ra_lines = np.arange(-150, 181, 30)
		for ra in ra_lines:
			dec = np.linspace(-90, 90, 500)
			ra_array = np.full_like(dec, ra)
			ax.plot(np.radians(ra_array), np.radians(dec), color='gray', alpha=0.75, linewidth=0.5)

		# draw the parallels
		dec_lines = np.arange(-75, 90, 15)
		for dec in dec_lines:
			ra = np.linspace(-180, 180, 1000)
			dec_array = np.full_like(ra, dec)
			ax.plot(np.radians(ra), np.radians(dec_array), color='gray', alpha=0.75, linewidth=0.5)

		# set the tick labels
		tick_labels = ['14h', '16h', '18h', '20h', '22h', '0h', '2h', '4h', '6h', '8h', '10h']
		ax.set_xticklabels(tick_labels)		

		########## BEGIN TESTING ##########
		'''
		rectangle = Grapheus.make_rectangle(270, -30, 30, 15)
		ax.add_patch(rectangle)

		circle = Grapheus.make_circle(52.759, -35.611, 3)
		ax.add_patch(circle)

		rectangle_ra, rectangle_dec = Grapheus.make_rectangular_footprint(15, -30, 30, 15)
		ax.plot(np.radians(rectangle_ra), np.radians(rectangle_dec), color='red', linewidth=1)
		ax.fill(np.radians(rectangle_ra), np.radians(rectangle_dec), color='gray')

		tuc47_ra, tuc47_dec = Grapheus.make_circular_footprint(6.207, -72.292, 0.727)
		ax.plot(np.radians(tuc47_ra), np.radians(tuc47_dec), color='red', linewidth=1.5)
		ax.fill(np.radians(tuc47_ra), np.radians(tuc47_dec), color='gray', alpha=0.5)

		circle_ra, circle_dec = Grapheus.make_circular_footprint(120, -30, 5)
		ax.plot(np.radians(circle_ra), np.radians(circle_dec), color='red', linewidth=1.5)
		ax.fill(np.radians(circle_ra), np.radians(circle_dec), color='gray', alpha=0.5)

		survey_ra = np.linspace(-np.pi, np.pi, 1000)
		limit_dec = Skopeus.get_declination_limit()
		ax.fill_between(survey_ra, -np.pi/2, np.radians(limit_dec), color='tab:blue', alpha=0.25, edgecolor='none')
		ax.plot(survey_ra, np.full_like(survey_ra, np.radians(limit_dec)), color='tab:red', linewidth=1.25, linestyle='-', label='declination limit', zorder=9)
		'''
		########## END TESTING ##########

		plt.legend()
		plt.tight_layout()

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def make_circular_footprint(ra, dec, radius):

		center = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
		radius = radius * u.deg

		# generate position angles around the circle
		position_angles = np.linspace(0, 360, 361) * u.deg

		# generate spherical circle coordinates
		circle_coordinates = center.directional_offset_by(position_angles, radius)

		# transform right ascension into skymap projection coordinates
		circle_ra = np.remainder(circle_coordinates.ra.deg + 180, 360) - 180
		circle_dec = circle_coordinates.dec.deg

		return circle_ra, circle_dec

	@staticmethod
	def make_rectangular_footprint(center_ra, center_dec, width, height):

		center = SkyCoord(ra=center_ra*u.deg, dec=center_dec*u.deg, frame='icrs')

		offset_frame = SkyOffsetFrame(origin=center)

		x = np.linspace(-width/2, width/2, 100)
		y = np.linspace(-height/2, height/2, 100)

		lon = np.concatenate([x, np.full(100, width/2), x[::-1], np.full(100, -width/2)])
		lat = np.concatenate([np.full(100, -height/2), y, np.full(100, height/2), y[::-1]])

		footprint = SkyCoord(lon=lon*u.deg, lat=lat*u.deg, frame=offset_frame).icrs

		rectangle_ra = np.remainder(footprint.ra.deg + 180, 360) - 180
		rectangle_dec = footprint.dec.deg

		return rectangle_ra, rectangle_dec

	@staticmethod
	def make_circle(center_ra, center_dec, radius):

		circle_ra = np.remainder(center_ra + 180, 360) - 180
		circle_ra = np.radians(circle_ra)

		circle_dec = np.radians(center_dec)

		center = (circle_ra, circle_dec)

		radius = np.radians(radius)

		circle = mpatches.Circle(center, radius, linewidth=0.5, edgecolor='r', facecolor='gray', alpha=0.5)

		return circle

	@staticmethod
	def make_rectangle(anchor_ra, anchor_dec, width, height):

		rectangle_ra = np.remainder(anchor_ra + 180, 360) - 180
		rectangle_ra = np.radians(rectangle_ra)

		rectangle_dec = np.radians(anchor_dec)

		anchor = (rectangle_ra, rectangle_dec)

		width = np.radians(width)
		height = np.radians(height)

		rectangle = mpatches.Rectangle(anchor, width, height, linewidth=0.5, edgecolor='r', facecolor='gray', alpha=0.5)

		return rectangle

	@staticmethod
	def draw_commissioning_skymap(save_path, save_figure):

		# set up the figure
		fig = plt.figure(figsize=Configuration.FIGURE_SIZE)
		ax = fig.add_subplot(111, projection=Configuration.PROJECTION)

		# draw the galactic plane
		if Configuration.PLOT_GALACTIC_PLANE:
			ra_plot, dec_plot = Logistes.get_galactic_plane(contour_level=0.)
			ax.plot(ra_plot, dec_plot, ls='--', lw=0.75, c='black')

		# draw the galactic contours
		if Configuration.PLOT_GALACTIC_CONTOUR:
			ra_neg_plot, dec_neg_plot = Logistes.get_galactic_plane(contour_level=-15)
			ra_pos_plot, dec_pos_plot = Logistes.get_galactic_plane(contour_level=15.)
			ax.plot(ra_neg_plot, dec_neg_plot, ls='--', lw=0.5, c='black')
			ax.plot(ra_pos_plot, dec_pos_plot, ls='--', lw=0.5, c='black')

		# draw the galactic center
		if Configuration.PLOT_GALACTIC_CENTER:
			gc_ra_rad, gc_dec_rad = Logistes.get_galactic_center()
			ax.scatter(gc_ra_rad, gc_dec_rad, marker='x', s=100, color='black', zorder=10)

		# draw the meridians
		ra_lines = np.arange(-150, 181, 30)
		for ra in ra_lines:
			dec = np.linspace(-90, 90, 500)
			ra_array = np.full_like(dec, ra)
			ax.plot(np.radians(ra_array), np.radians(dec), color='gray', alpha=0.75, linewidth=0.5)

		# draw the parallels
		dec_lines = np.arange(-75, 90, 15)
		for dec in dec_lines:
			ra = np.linspace(-180, 180, 1000)
			dec_array = np.full_like(ra, dec)
			ax.plot(np.radians(ra), np.radians(dec_array), color='gray', alpha=0.75, linewidth=0.5)

		# set the tick labels
		tick_labels = ['14h', '16h', '18h', '20h', '22h', '0h', '2h', '4h', '6h', '8h', '10h']
		ax.set_xticklabels(tick_labels)

		# draw the survey footprint
		survey_ra = np.linspace(-np.pi, np.pi, 1000)
		limit_dec = Skopeus.get_declination_limit()
		ax.fill_between(survey_ra, -np.pi/2, np.radians(limit_dec), color='tab:blue', alpha=0.25, edgecolor='none')
		ax.plot(survey_ra, np.full_like(survey_ra, np.radians(limit_dec)), color='tab:red', linewidth=1.25, linestyle='-', zorder=9)

		# draw the commissioning fields
		com_field_ra_list, com_field_dec_list, gw_field_ra_list, gw_field_dec_list, grb_field_ra_list, grb_field_dec_list, ot_field_ra_list, ot_field_dec_list, sn_field_ra_list, sn_field_dec_list, xrb_field_ra_list, xrb_field_dec_list = Skopeus.get_commissioning_fields()
		ax.scatter(com_field_ra_list, com_field_dec_list, marker='s', s=50, color='black', edgecolor='none', label='Comm')
		ax.scatter(grb_field_ra_list, grb_field_dec_list, marker='s', s=50, color='green', edgecolor='none', label='GRB')
		ax.scatter(gw_field_ra_list, gw_field_dec_list, marker='s', s=50, color='purple', edgecolor='none', label='GW')
		ax.scatter(ot_field_ra_list, ot_field_dec_list, marker='s', s=50, color='red', edgecolor='none', label='OT')
		ax.scatter(sn_field_ra_list, sn_field_dec_list, marker='s', s=50, color='orange', edgecolor='none', label='SN')
		ax.scatter(xrb_field_ra_list, xrb_field_dec_list, marker='s', s=50, color='blue', edgecolor='none', label='XRB')

		plt.legend()
		plt.tight_layout()

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def draw_survey_skymap(save_path, save_figure):

		# set up the figure
		fig = plt.figure(figsize=Configuration.FIGURE_SIZE)
		ax = fig.add_subplot(111, projection=Configuration.PROJECTION)

		# draw the galactic plane
		if Configuration.PLOT_GALACTIC_PLANE:
			ra_plot, dec_plot = Logistes.get_galactic_plane(contour_level=0.)
			ax.plot(ra_plot, dec_plot, ls='--', lw=0.75, c='black')

		# draw the galactic contours
		if Configuration.PLOT_GALACTIC_CONTOUR:
			ra_neg_plot, dec_neg_plot = Logistes.get_galactic_plane(contour_level=-15.)
			ra_pos_plot, dec_pos_plot = Logistes.get_galactic_plane(contour_level=15.)
			ax.plot(ra_neg_plot, dec_neg_plot, ls='--', lw=0.5, c='black')
			ax.plot(ra_pos_plot, dec_pos_plot, ls='--', lw=0.5, c='black')

		# draw the galactic center
		if Configuration.PLOT_GALACTIC_CENTER:
			gc_ra_rad, gc_dec_rad = Logistes.get_galactic_center()
			ax.scatter(gc_ra_rad, gc_dec_rad, marker='x', s=100, color='black', zorder=10)

		# draw the meridians
		ra_lines = np.arange(-150, 181, 30)
		for ra in ra_lines:
			dec = np.linspace(-90, 90, 500)
			ra_array = np.full_like(dec, ra)
			ax.plot(np.radians(ra_array), np.radians(dec), color='gray', alpha=0.75, linewidth=0.5)

		# draw the parallels
		dec_lines = np.arange(-75, 90, 15)
		for dec in dec_lines:
			ra = np.linspace(-180, 180, 1000)
			dec_array = np.full_like(ra, dec)
			ax.plot(np.radians(ra), np.radians(dec_array), color='gray', alpha=0.75, linewidth=0.5)

		# set the tick labels
		tick_labels = ['14h', '16h', '18h', '20h', '22h', '0h', '2h', '4h', '6h', '8h', '10h']
		ax.set_xticklabels(tick_labels)

		# draw the survey fields
		gal_ra_list, gal_dec_list, egal_ra_list, egal_dec_list = Skopeus.get_survey_fields()
		ax.scatter(gal_ra_list, gal_dec_list, marker='.', s=10, color='tab:pink', label='galactic fields')
		ax.scatter(egal_ra_list, egal_dec_list, marker='.', s=10, color='tab:blue', label='extragalactic fields')

		# draw the survey footprint
		survey_ra = np.linspace(-np.pi, np.pi, 1000)
		limit_dec = Skopeus.get_declination_limit()
		ax.plot(survey_ra, np.full_like(survey_ra, np.radians(limit_dec)), color='tab:red', linewidth=1.25, linestyle='-', zorder=9)

		plt.legend()
		plt.tight_layout()

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def draw_tns_skymap(file_path, save_path, save_figure):

		# set up the figure
		fig = plt.figure(figsize=Configuration.FIGURE_SIZE)
		ax = fig.add_subplot(111, projection=Configuration.PROJECTION)

		# draw the galactic plane
		if Configuration.PLOT_GALACTIC_PLANE:
			ra_plot, dec_plot = Logistes.get_galactic_plane(contour_level=0.)
			ax.plot(ra_plot, dec_plot, ls='--', lw=0.75, c='black')

		# draw the galactic contours
		if Configuration.PLOT_GALACTIC_CONTOUR:
			ra_neg_plot, dec_neg_plot = Logistes.get_galactic_plane(contour_level=-15.)
			ra_pos_plot, dec_pos_plot = Logistes.get_galactic_plane(contour_level=15.)
			ax.plot(ra_neg_plot, dec_neg_plot, ls='--', lw=0.5, c='black')
			ax.plot(ra_pos_plot, dec_pos_plot, ls='--', lw=0.5, c='black')

		# draw the galactic center
		if Configuration.PLOT_GALACTIC_CENTER:
			gc_ra_rad, gc_dec_rad = Logistes.get_galactic_center()
			ax.scatter(gc_ra_rad, gc_dec_rad, marker='x', s=100, color='black', zorder=10)

		# draw the meridians
		ra_lines = np.arange(-150, 181, 30)
		for ra in ra_lines:
			dec = np.linspace(-90, 90, 500)
			ra_array = np.full_like(dec, ra)
			ax.plot(np.radians(ra_array), np.radians(dec), color='gray', alpha=0.75, linewidth=0.5)

		# draw the parallels
		dec_lines = np.arange(-75, 90, 15)
		for dec in dec_lines:
			ra = np.linspace(-180, 180, 1000)
			dec_array = np.full_like(ra, dec)
			ax.plot(np.radians(ra), np.radians(dec_array), color='gray', alpha=0.75, linewidth=0.5)

		# set the tick labels
		tick_labels = ['14h', '16h', '18h', '20h', '22h', '0h', '2h', '4h', '6h', '8h', '10h']
		ax.set_xticklabels(tick_labels)

		tns_file = open(file_path, 'r')
		for line in tns_file:
			line = line.split(",")

			if line[0] != '"ID"':
				ra = line[2].strip('"')
				dec = line[3].strip('"')
				type = line[4].strip('"')

				coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')

				# convert the right ascension values into range [-180, +180]
				ra_shift = np.remainder(coord.ra.deg + 180, 360) - 180

				# convert the coordinates to radians
				ra_rad = np.radians(ra_shift)
				dec_rad = np.radians(coord.dec.deg)

				# draw the objects
				if (type == 'SN Ia') or (type == 'SN Ia-91bg-like') or (type == 'SN Ia-91T-like') or (type == 'SN Ia-CSM'):
					ax.scatter(ra_rad, dec_rad, c='tab:red', s=25, label='Thermonuclear')
				elif (type == 'SN II') or (type == 'SN Iax[02cx-like]') or (type == 'SN Ic') or (type == 'SN Ib/c') or (type == 'SN IIP') or (type == 'SN Ib') or (type == 'SN IIb') or (type == 'SN Ic-BL') or (type == 'SN IIn'):
					ax.scatter(ra_rad, dec_rad, c='tab:orange', s=25, label='CCSNe')
				elif type == 'SLSN-I':
					ax.scatter(ra_rad, dec_rad, c='gold', marker='*', s=50, edgecolor='red', linewidth=0.5, label='SLSNe')
				else:
					ax.scatter(ra_rad, dec_rad, c='tab:gray', s=3, label='OT')

		plt.legend()
		plt.tight_layout()

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()
