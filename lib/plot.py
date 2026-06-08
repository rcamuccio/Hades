from config import Configuration
from lib.calculator import Calculator

from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.visualization import ZScaleInterval
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
mpl.rc('font', family=Configuration.FONT_NAME, size=Configuration.FONT_SIZE)

class Plot:

	@staticmethod
	def airmass(xlist, ylist):

		plt.clf()
		font = {'fontname':'Monospace', 'size':12}

		plt.figure(figsize=(10, 10))
		plt.plot(xlist, ylist, color='black')

		plt.title('Air Mass', **font)
		plt.xlabel('Time [JD]', **font)
		plt.ylabel('Air Mass', **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig('plot-airmass.png', dpi=400)
		plt.close()

	@staticmethod
	def color_magnitude(color_list, delta_list, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a color-magnitude diagram.

		:parameter color_list - The list of catalog colors
		:parameter delta_list - The list of difference magnitudes
		:parameter save_path - The save path of the diagram
		:parameter save_figure - A toggle to save or display the diagram

		:return - Nothing is returned
		'''

		# calculate the linear fit
		linear_fit, slope, intercept, delta_slope, delta_intercept = Calculator.unweighted_fit(color_list, delta_list)

		# draw the diagram
		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.scatter(color_list, delta_list, s=1, color='gray')
		plt.plot(color_list, linear_fit, ls='--', lw=0.75, color='black')
		plt.title('Color-Magnitude Diagram')
		plt.xlabel('Color index [mag]')
		plt.ylabel('Delta magnitude [mag]')

		# save or display the diagram
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def extinction(xlist, ylist, yfit):

		plt.clf()
		font = {'fontname':'Monospace', 'size':12}

		plt.figure(figsize=(10, 10))
		plt.scatter(xlist, ylist, s=1, color='gray')
		plt.plot(xlist, yfit, color='blue')

		plt.title('Extinction', **font)
		plt.xlabel('Air mass', **font)
		plt.ylabel('Instrumental magnitude [mag]', **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig('plot-extinction.png', dpi=400)
		plt.close()

	@staticmethod
	def field(frame_data, apertures, boxes, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws an image with annotations.

		:parameter frame_data - The image data array to be drawn
		:parameter apertures - A set of apertures for optional display
		:parameter boxes - A set of mask boxes for optional display
		:parameter save_path - The save path of the diagram
		:parameter save_figure - A toggle to save or display the diagram

		:return - Nothing is returned
		'''

		# calculate the range limits
		interval = ZScaleInterval()
		vmin, vmax = interval.get_limits(frame_data)

		# draw the diagram
		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.imshow(frame_data, cmap=Configuration.CMAP, origin='lower', vmin=vmin, vmax=vmax)
		plt.xlabel('x Pixel')
		plt.ylabel('y Pixel')
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
	def growth_radius(xlist, ylist):

		plt.clf()
		font = {'fontname':'Monospace', 'size':12}

		plt.figure(figsize=(10, 10))
		plt.plot(xlist, ylist, color='black')

		plt.title('Growth Radius', **font)
		plt.xlabel('Time [JD]', **font)
		plt.ylabel('Radius [px]', **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig('plot-growth-radius.png', dpi=400)
		plt.close()

	@staticmethod
	def stellar_flux_histogram(flux_list, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a histogram of stellar fluxes.

		:parameter flux_list - The list of stellar fluxes
		:parameter save_path - The save path of the diagram
		:parameter save_figure - A toggle to save or display the diagram

		:return - Nothing is returned
		'''

		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.hist(flux_list, bins=Configuration.HISTOGRAM_BINS, range=(0, 2.5e5), histtype=Configuration.HISTOGRAM_TYPE)
		plt.yscale(Configuration.HISTOGRAM_SCALE)
		plt.xlabel('Flux [ADU]')
		plt.ylabel('Count')
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def stellar_gaia_histogram(r_mag_list, g_mag_list, b_mag_list, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a histogram of stellar Gaia magnitudes.

		:parameter r_mag_list - The list of Gaia Rp magnitudes
		:parameter g_mag_list - The list of Gaia Gp magnitudes
		:parameter b_mag_list - The list of Gaia Bp magnitudes
		:parameter save_path - The save path of the diagram
		:parameter save_figure - A toggle to save or display the diagram

		:return - Nothing is returned
		'''

		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.hist(r_mag_list, bins=Configuration.HISTOGRAM_BINS, histtype=Configuration.HISTOGRAM_TYPE, color='red', label='Gaia Rp')
		plt.hist(g_mag_list, bins=Configuration.HISTOGRAM_BINS, histtype=Configuration.HISTOGRAM_TYPE, color='green', label='Gaia G')
		plt.hist(b_mag_list, bins=Configuration.HISTOGRAM_BINS, histtype=Configuration.HISTOGRAM_TYPE, color='blue', label='Gaia Bp')
		plt.yscale(Configuration.HISTOGRAM_SCALE)
		plt.xlabel('Gaia magnitude [mag]')
		plt.ylabel('Count')
		plt.legend()
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def stellar_magnitude_histogram(mag_list, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a histogram of stellar magnitudes.

		:parameter mag_list - The list of stellar magnitudes
		:parameter save_path - The save path of the diagram
		:parameter save_figure - A toggle to save or display the diagram

		:return - Nothing is returned
		'''

		plt.figure(figsize=Configuration.FIGURE_SIZE)
		plt.hist(mag_list, bins=Configuration.HISTOGRAM_BINS, range=(-12, 4), histtype=Configuration.HISTOGRAM_TYPE)
		plt.yscale(Configuration.HISTOGRAM_SCALE)
		plt.xlabel('Instrumental magnitude [mag]')
		plt.ylabel('Count')
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def stellar_histogram(flux_list, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a histogram of stellar fluxes.

		:parameter flux_list - The list of stellar fluxes
		:parameter save_path - The save path of the diagram
		:parameter save_figure - A toggle to save or display the diagram

		:return - Nothing is returned
		'''

		plt.figure(figsize=Configuration.FIGURE_SIZE)

		if Configuration.HISTOGRAM_FLUX == 'flux':
			plt.hist(flux_list, bins=Configuration.HISTOGRAM_BINS, range=(0, 2.5e5), histtype=Configuration.HISTOGRAM_TYPE)
			plt.xlabel('Flux [ADU]')
		elif Configuration.HISTOGRAM_FLUX == 'mag':
			plt.hist(flux_list, bins=Configuration.HISTOGRAM_BINS, histtype=Configuration.HISTOGRAM_TYPE)
			plt.xlabel('Instrumental magnitude [mag]')
		elif Configuration.HISTOGRAM_FLUX == 'gaia':
			plt.hist(flux_list, bins=Configuration.HISTOGRAM_BINS, histtype=Configuration.HISTOGRAM_TYPE)
			plt.xlabel('Gaia magnitude [mag]')
		plt.yscale(Configuration.HISTOGRAM_SCALE)
		plt.xlabel('Flux [ADU]')
		plt.ylabel('Count')
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def frame_histogram(frame_data, save_path, save_figure=Configuration.SAVE_FIGURE):
		'''This function draws a histogram of pixel values.

		:parameter frame_data - The image data array
		:parameter save_path - The save path of the diagram
		:parameter save_figure - A toggle to save or display the diagram

		:return - Nothing is returned
		'''

		# calculate the histogram limits
		data_mn, data_md, data_sd = sigma_clipped_stats(frame_data, sigma=Configuration.SIG_BKG)
		xmin = data_md - Configuration.HISTOGRAM_LIMIT * data_sd
		xmax = data_md + Configuration.HISTOGRAM_LIMIT * data_sd
		
		# draw the plot
		plf.figure(figsize=Configuration.FIGURE_SIZE)
		plt.hist(frame_data.flatten(), bins=Configuration.HISTOGRAM_BINS, range=(xmin, xmax), histtype=Configuration.HISTOGRAM_TYPE)
		plt.yscale(Configuration.HISTOGRAM_SCALE)
		plt.xlabel('Pixel Value [ADU]', **font)
		plt.ylabel('Count', **font)
		if save_figure:
			plt.savefig(save_path, dpi=Configuration.DPI)
			plt.close()
		else:
			plt.show()

	@staticmethod
	def lightcurve(xlist, ylist, elist):

		plt.clf()
		font = {'fontname':'Monospace', 'size':12}

		plt.figure(figsize=(10, 10))
		plt.errorbar(xlist, ylist, yerr=elist, fmt='o', linewidth=0.5, markersize=0.5, capsize=2, capthick=0.5)

		#plt.ylim(13, 19)
		plt.title('Time Series of ' + name, **font)
		plt.xlabel('Time [JD]', **font)
		plt.ylabel('Magnitude [mag]', **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig('plot-' + name + '-timeseries.png', dpi=400)
		plt.close()

	@staticmethod
	def seeing(xlist, ylist, type):

		plt.clf()
		font = {'fontname':'Monospace', 'size':12}

		plt.figure(figsize=(10, 10))
		plt.plot(xlist, ylist, color='black')

		if type == 'physical':
			plt.title('Seeing (physical)', **font)
			plt.xlabel('Time [JD]', **font)
			plt.ylabel('Seeing [px]', **font)
		elif type == 'sky':
			plt.title('Seeing (sky)', **font)
			plt.xlabel('Time [JD]', **font)
			plt.ylabel('Seeing [arcsec]', **font)

		plt.xticks(**font)
		plt.yticks(**font)
		
		plt.savefig('plot-seeing.png', dpi=400)
		plt.close()

	@staticmethod
	def skymap(ra_deg=None, dec_deg=None):
		'''This function plots an all-sky map of coordinates.

		:parameter ra_deg - The right ascension of one or more objects in degrees [0, 360)
		:parameter dec_deg - The declination of one or more objects in degrees [-90, 90]

		:return - Nothing is returned
		'''

		# generate the plot
		plt.clf()
		font = {'fontname':Configuration.FONT_NAME, 'size':Configuration.FONT_SIZE}
		fig = plt.figure(figsize=Configuration.FIGURE_SIZE)
		ax = fig.add_subplot(111, projection=Configuration.PROJECTION)

		# draw the coordinates
		if (ra_deg != None) and (dec_deg != None):
			# ensure the input coordinates are array-like
			ra_deg = np.asarray(ra_deg)
			dec_deg = np.asarray(dec_deg)

			# convert RA into range [-180, +180]
			ra_shift = np.remainder(ra_deg + 180, 360) - 180

			# convert the coordinates to radians
			ra_rad = np.radians(ra_shift)
			dec_rad = np.radians(dec_deg)

			# draw the objects
			ax.scatter(ra_rad, dec_rad, c='black', s=3)

		# draw the galactic plane
		if Configuration.GALACTIC_PLANE:
			# generate the plane
			l = np.linspace(0, 360, 2000) * u.deg
			b = np.zeros_like(l)
			gal_plane = SkyCoord(l=l, b=b, frame='galactic')

			# convert the plane to equatorial coordinates
			eq_plane = gal_plane.icrs
			eq_plane_ra = eq_plane.ra.degree
			eq_plane_dec = eq_plane.dec.degree

			# convert RA into range [-180, +180]
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

			# draw the plane
			ax.plot(ra_plot, dec_plot, ls='--', lw=0.75, c='black', zorder=9)

		# draw the galactic center
		if Configuration.GALACTIC_CENTER:
			gc = SkyCoord(l=0*u.deg, b=0*u.deg, frame='galactic').icrs
			gc_ra = gc.ra.degree
			gc_dec = gc.dec.degree
			gc_ra_shift = np.remainder(gc_ra + 180, 360) - 180
			gc_ra_rad = np.radians(gc_ra_shift)
			gc_dec_rad = np.radians(gc_dec)
			ax.scatter(gc_ra_rad, gc_dec_rad, marker='+', s=200, color='black', zorder=10)

		if Configuration.SKYMAP_MODE == 'survey_fields':
			field_path = Configuration.PLOUTON_DIRECTORY + 'fields/toros_fields' + Configuration.TABLE_EXTENSION
			field_file = open(field_path, 'r')
			gal_ra_list = []
			gal_dec_list = []
			egal_ra_list = []
			egal_dec_list = []
			for ln in field_file:
				ln = ln.split()
				if ln[0] != 'field_id':
					ra = float(ln[1])
					ra = np.remainder(ra + 180, 360) - 180
					ra = np.radians(ra)
					dec = float(ln[2])
					dec = np.radians(dec)
					l = float(ln[3])
					b = float(ln[4])
					if abs(b) <= Configuration.GALACTIC_PLANE:
						gal_ra_list.append(ra)
						gal_dec_list.append(dec)
					else:
						egal_ra_list.append(ra)
						egal_dec_list.append(dec)
			ax.scatter(gal_ra_list, gal_dec_list, marker='.', s=10, color='red')
			ax.scatter(egal_ra_list, egal_dec_list, marker='.', s=10, color='blue')

		elif Configuration.SKYMAP_MODE == 'comm_fields':
			field_path = Configuration.PLOUTON_DIRECTORY + 'fields/comm_fields' + Configuration.TABLE_EXTENSION
			field_file = open(field_path, 'r')
			com_field_ra_list = []
			com_field_dec_list = []
			grb_field_ra_list = []
			grb_field_dec_list = []
			gw_field_ra_list = []
			gw_field_dec_list = []
			ot_field_ra_list = []
			ot_field_dec_list = []
			sn_field_ra_list = []
			sn_field_dec_list = []
			xrb_field_ra_list = []
			xrb_field_dec_list = []
			for ln in field_file:
				ln = ln.split()
				field_id = ln[0]
				field_ra = float(ln[1])
				field_ra = np.remainder(field_ra + 180, 360) - 180
				field_ra = float(np.radians(field_ra))
				field_dec = float(ln[2])
				field_dec = float(np.radians(field_dec))
				field_type = ln[3]
				event_type = ln[4]
				if field_type == 'C':
					com_field_ra_list.append(field_ra)
					com_field_dec_list.append(field_dec)
				else:
					if event_type == 'GW':
						gw_field_ra_list.append(field_ra)
						gw_field_dec_list.append(field_dec)
					elif (event_type == 'LGRB') or (event_type == 'SGRB'):
						grb_field_ra_list.append(field_ra)
						grb_field_dec_list.append(field_dec)
					elif event_type == 'OT':
						ot_field_ra_list.append(field_ra)
						ot_field_dec_list.append(field_dec)
					elif (event_type == 'SNIA') or (event_type == 'SNII'):
						sn_field_ra_list.append(field_ra)
						sn_field_dec_list.append(field_dec)
					elif event_type == 'XRB':
						xrb_field_ra_list.append(field_ra)
						xrb_field_dec_list.append(field_dec)
			ax.scatter(com_field_ra_list, com_field_dec_list, marker='s', s=100, color='black', label='Commissioning')
			ax.scatter(grb_field_ra_list, grb_field_dec_list, marker='s', s=100, color='green', label='GRB')
			ax.scatter(ot_field_ra_list, ot_field_dec_list, marker='s', s=100, color='red', label='OT')
			ax.scatter(sn_field_ra_list, sn_field_dec_list, marker='s', s=100, color='orange', label='SN')
			ax.scatter(xrb_field_ra_list, xrb_field_dec_list, marker='s', s=100, color='blue', label='XRB')

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

		plt.legend()
		plt.tight_layout()
		plt.show()