from config import Configuration
from libraries.calculator import Calculator
from libraries.plotter import Plotter
from libraries.reducer import Reducer

from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import hstack, Table
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from photutils.aperture import aperture_photometry, CircularAnnulus, CircularAperture
import math


class Photometer:

	@staticmethod
	def extract_pix_sources(object_frame):

		object_name = object_frame[:-4]
		frame = fits.open(object_frame)
		frame_data = frame[0].data
		frame_header = frame[0].header

		mask, boxes = Reducer.make_mask(object_frame)

		print('Extracting sources (xp, yp) for', object_frame)
		mean, median, std = sigma_clipped_stats(frame_data, mask=mask, sigma=Configuration.SIGMA_BKG)

		fwhm = 6.0
		threshold = Configuration.SIGMA_SRC * std

		daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)
		table = daofind(frame_data, mask=mask)

		return table

	@staticmethod
	def extract_sky_sources(object_frame):

		object_name = object_frame[:-4]
		frame = fits.open(object_frame)
		frame_data = frame[0].data
		frame_header = frame[0].header

		wcs = WCS(frame_header)

		mask, boxes = Reducer.make_mask(object_frame)

		print('Extracting sources (\u03B1, \u03B4) for', object_frame)
		mean, median, std = sigma_clipped_stats(frame_data, mask=mask, sigma=Configuration.SIGMA_BKG)

		fwhm = 6.0
		threshold = Configuration.SIGMA_SRC * std

		daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)
		table = daofind(frame_data, mask=mask)

		pix_coord_list = []
		ra_list = []
		dec_list = []

		for item in table:
			xp = item['xcentroid']
			yp = item['ycentroid']

			sky_positions = pixel_to_skycoord(xp, yp, wcs=wcs)

			ra = sky_positions.ra.deg
			dec = sky_positions.dec.deg

			ra_list.append(ra)
			dec_list.append(dec)

		table['ra'] = ra_list
		table['dec'] = dec_list

		return table

	@staticmethod
	def match_catalogs(source_table, query_table):

		print('Matching catalogs')
		coo_source = SkyCoord(source_table['ra']*u.deg, source_table['dec']*u.deg)

		if (Configuration.CATALOG == 'gaia-cone') or (Configuration.CATALOG == 'gaia-square'):
			coo_query = SkyCoord(query_table['ra'], query_table['dec'])

		elif (Configuration.CATALOG == 'ps1'):
			coo_query = SkyCoord(query_table['raMean']*u.deg, query_table['decMean']*u.deg)

		idx, d2d, d3d = coo_query.match_to_catalog_sky(coo_source)

		t = d2d < 1*u.arcsec

		query_table['sep_flag'] = t
		query_table['idx'] = idx

		ct = query_table.dtype
		names, dtypes = zip(*ct.descr)
		match_table = Table(names=names, dtype=dtypes)

		for row, sep_flag in enumerate(query_table['sep_flag']):
			if sep_flag:
				match_table.add_row(query_table[row].as_void())

		return match_table

	@staticmethod
	def photometry_field(object_frame, match_table, survey='gaia'):

		print('Performing photometery on', object_frame)
		object_name = object_frame[:-4]
		frame = fits.open(object_frame)
		frame_data = frame[0].data
		frame_header = frame[0].header
		
		exptime = frame_header['EXPTIME']
		wcs = WCS(frame_header)

		# --- Create frame mask
		mask, boxes = Reducer.make_mask(object_frame)

		# --- Calculate background statistics
		mean, median, std = sigma_clipped_stats(frame_data, mask=mask, sigma=Configuration.SIGMA_BKG)
		print(mean, median, std)

		# --- Grab sky coordinates from catalog
		sky_coord_list = []
		for item in match_table:

			ra = item['ra']
			dec = item['dec']

			coord_tuple = (ra, dec)
			sky_coord_list.append(coord_tuple)

		# --- Convert sky coordinates to pixel coordinates
		sky_positions = SkyCoord(sky_coord_list, unit='deg')
		pix_positions = skycoord_to_pixel(sky_positions, wcs=wcs)

		pix_coord_list = []
		for px in range(len(pix_positions[0])):
			xp = pix_positions[0][px]
			yp = pix_positions[1][px]
			coord_tuple = (xp, yp)
			pix_coord_list.append(coord_tuple)

		# --- Create apertures and annuli
		apertures = CircularAperture(pix_coord_list, r=Configuration.RAD_AP)
		annuli = CircularAnnulus(pix_coord_list, r_in=Configuration.RAD_AN_IN, r_out=Configuration.RAD_AN_OUT)
		apers = [apertures, annuli]

		# --- Conduct aperture photometry at all catalog positions
		phot_table = aperture_photometry(frame_data, apers, mask=mask, wcs=wcs)

		aperture_area = apertures.area
		annulus_area = annuli.area

		bkg_mean = phot_table['aperture_sum_1'] / annulus_area
		phot_table['annulus_mean'] = bkg_mean

		bkg_sum = bkg_mean * aperture_area
		phot_table['aperture_bkg_sum'] = bkg_sum

		final_sum = phot_table['aperture_sum_0'] - bkg_sum
		phot_table['res_aperture_sum'] = final_sum

		# --- Merge catalog and photometry tables into master table
		master_table = hstack([match_table, phot_table])

		# --- Calculate fluxes and magnitudes
		flux_list = []
		flux_error_list = []
		flux_flag_list = []

		inst_mag_list = []
		inst_mag_error_list = []

		bkg_val = median + std

		for item in master_table:

			# --- Filter flux values below mean (and zero)
			if item['res_aperture_sum'] <= bkg_val:
				flux_flag_list.append(True)
				flux = bkg_val

			else:
				flux_flag_list.append(False)
				flux = item['res_aperture_sum']

			print(flux)

			flux_list.append(flux)

			# --- Calculate flux error (Poisson + RMS)
			flux_error = math.sqrt(flux + (aperture_area * (1 + (math.pi * aperture_area) / (2 * annulus_area))*(std**2)))
			flux_error_list.append(flux_error)

			# --- Calculate instrumental magnitude (exposure time correction)
			inst_mag = (-2.5*math.log10(flux)) + (2.5*math.log10(exptime))
			inst_mag_list.append(inst_mag)

			# --- Calculate instrumental magnitude error
			inst_mag_error = (2.5 * flux_error) / (math.log(10) * flux)
			inst_mag_error_list.append(inst_mag_error)

		master_table['flux'] = flux_list
		master_table['flux_error'] = flux_error_list
		master_table['flux_flag'] = flux_flag_list
		master_table['inst_mag'] = inst_mag_list
		master_table['inst_mag_error'] = inst_mag_error_list

		# --- Calculate colors and delta magnitudes
		color_list = []
		color_flag_list = []
		delta_mag_list = []

		for item in master_table:

			if survey == "gaia":

				if item["phot_bp_n_obs"] == 0 or item["phot_rp_n_obs"] == 0:
					
					color_flag_list.append(True)
					color = 0
					cat_mag = 0

				else:
					
					color_flag_list.append(False)
					b = item["phot_bp_mean_mag"]
					#b = b.unmasked.value
					r = item["phot_rp_mean_mag"]
					#r = r.unmasked.value

					color = b - r
					cat_mag = item['phot_g_mean_mag']

			else:

				if item["nr"] == 0 or item["ni"] == 0:

					color_flag_list.append(True)
					color = 0
					cat_mag = 0

				else:

					color_flag_list.append(False)
					r = item["rMeanPSFMag"]
					i = item["iMeanPSFMag"]

					color = r - i
					cat_mag = item["rMeanPSFMag"]

			color_list.append(color)

			inst_mag = item["inst_mag"]
			delta_mag = cat_mag - inst_mag
			delta_mag_list.append(delta_mag)

		master_table["color"] = color_list
		master_table["color_flag"] = color_flag_list
		master_table["delta_mag"] = delta_mag_list

		# --- Calculate transform and zero point
		color_list = []
		delta_mag_list = []

		for item in master_table:

			if (item["flux_flag"] == True) or (item["color_flag"] == True):
				pass

			else:
				color = item["color"]
				color_list.append(color)

				delta_mag = item["delta_mag"]
				delta_mag_list.append(delta_mag)

		yfit, slope, intercept, delta_slope, delta_intercept = Calculator.unweighted_fit(color_list, delta_mag_list)

		print("Transform:", slope, "+/-", delta_slope)
		print("ZP:", intercept, "+/-", delta_intercept)

		Plotter.plot_colormag(object_name, color_list, delta_mag_list, yfit)

		return master_table

	@staticmethod
	def photometry_pix_point(object_frame, xp, yp):

		object_name = object_frame[:-4]
		frame = fits.open(object_frame)
		frame_data = frame[0].data
		frame_header = frame[0].header
		wcs = WCS(frame_header)

		sigma = 3.0

		mask, boxes = reducer.make_mask(object_frame)
		mean, median, std = sigma_clipped_stats(frame_data, mask=mask, sigma=sigma)

		exptime = frame_header["EXPTIME"]
		dateobs = frame_header["DATE-OBS"]
		obstime = Time(dateobs)
		jd = frame_header["JD"]

		latitude = 25.995789
		longitude = -97.568954
		height = 11.5
		location = EarthLocation(lat=latitude, lon=longitude, height=height)

		altitude = 59.1
		airmass = Calculator.calculate_airmass(altitude)

		pix_coord_tuple = (xp, yp)

		ap = 11
		an_in = 18
		an_out = 26

		aperture = CircularAperture(pix_coord_tuple, r=ap)
		annulus = CircularAnnulus(pix_coord_tuple, r_in=an_in, r_out=an_out)
		apers = [aperture, annulus]

		target_table = aperture_photometry(frame_data, apers, mask=mask)

		aperture_area = aperture.area
		annulus_area = annulus.area

		bkg_mean = target_table["aperture_sum_1"] / annulus_area
		target_table["annulus_mean"] = bkg_mean

		bkg_sum = bkg_mean * aperture_area
		target_table["aperture_bkg_sum"] = bkg_sum

		final_sum = target_table["aperture_sum_0"] - bkg_sum
		target_table["res_aperture_sum"] = final_sum

		if target_table["res_aperture_sum"] <= mean:
			flux_flag = True
			if mean <= 0.0:
				flux = std
			else:
				flux = mean
		else:
			flux_flag = False
			flux = final_sum

		flux_error = math.sqrt(flux + (aperture_area * (1 + (math.pi * aperture_area) / (2 * annulus_area))*(std**2)))
		inst_mag = (-2.5*math.log10(flux)) + (2.5*math.log10(exptime))
		inst_mag_error = (2.5 * flux_error) / (math.log(10) * flux)

		target_table["flux"] = flux
		target_table["flux_error"] = flux_error
		target_table["flux_flag"] = flux_flag
		target_table["inst_mag"] = inst_mag
		target_table["inst_mag_error"] = inst_mag_error

		for item in target_table:
			flux = item["flux"]
			flux_error = item["flux_error"]
			inst_mag = item["inst_mag"]
			inst_mag_error = item["inst_mag_error"]

		return flux, flux_error, inst_mag, inst_mag_error, airmass, jd

	@staticmethod
	def photometry_sky_point(object_frame, ra, dec):

		object_name = object_frame[:-4]
		frame = fits.open(object_frame)
		frame_data = frame[0].data
		frame_header = frame[0].header
		wcs = WCS(frame_header)

		sigma = 3.0

		mask, boxes = Reducer.make_mask(object_frame)
		mean, median, std = sigma_clipped_stats(frame_data, mask=mask, sigma=sigma)
		
		exptime = frame_header["EXPTIME"]
		dateobs = frame_header["DATE-OBS"]
		obstime = Time(dateobs)
		jd = frame_header["JD"]

		latitude = 25.995789
		longitude = -97.568954
		height = 11.5
		location = EarthLocation(lat=latitude, lon=longitude, height=height)
		altitude = Calculator.calculate_altitude(location, obstime, ra, dec)
		airmass = Calculator.calculate_airmass(altitude)

		sky_position = SkyCoord(ra, dec, unit="deg")
		pix_position = skycoord_to_pixel(sky_position, wcs=wcs)
		xp = pix_position[0]
		yp = pix_position[1]
		pix_coord_tuple = (xp, yp)

		ap = 11
		an_in = 18
		an_out = 26

		aperture = CircularAperture(pix_coord_tuple, r=ap)
		annulus = CircularAnnulus(pix_coord_tuple, r_in=an_in, r_out=an_out)
		apers = [aperture, annulus]

		target_table = aperture_photometry(frame_data, apers, mask=mask, wcs=wcs)

		aperture_area = aperture.area
		annulus_area = annulus.area

		bkg_mean = target_table["aperture_sum_1"] / annulus_area
		target_table["annulus_mean"] = bkg_mean

		bkg_sum = bkg_mean * aperture_area
		target_table["aperture_bkg_sum"] = bkg_sum

		final_sum = target_table["aperture_sum_0"] - bkg_sum
		target_table["res_aperture_sum"] = final_sum

		if target_table["res_aperture_sum"] <= mean:
			flux_flag = True
			if mean <= 0.0:
				flux = std
			else:
				flux = mean
		else:
			flux_flag = False
			flux = final_sum

		flux_error = math.sqrt(flux + (aperture_area * (1 + (math.pi * aperture_area) / (2 * annulus_area))*(std**2)))
		inst_mag = (-2.5*math.log10(flux)) + (2.5*math.log10(exptime))
		inst_mag_error = (2.5 * flux_error) / (math.log(10) * flux)

		target_table["flux"] = flux
		target_table["flux_error"] = flux_error
		target_table["flux_flag"] = flux_flag
		target_table["inst_mag"] = inst_mag
		target_table["inst_mag_error"] = inst_mag_error

		for item in target_table:
			flux = item["flux"]
			flux_error = item["flux_error"]
			inst_mag = item["inst_mag"]
			inst_mag_error = item["inst_mag_error"]

		return flux, flux_error, inst_mag, inst_mag_error, airmass, jd

	@staticmethod
	def photometry_timeseries(object_list, ra, dec):

		flux_list = []
		flux_error_list = []
		inst_mag_list = []
		inst_mag_error_list = []
		airmass_list = []
		jd_list = []

		for item in object_list:
			flux, flux_error, inst_mag, inst_mag_error, airmass, jd = photometry_point(item, ra, dec)

			flux_list.append(flux)
			flux_error_list.append(flux_error)
			inst_mag_list.append(inst_mag)
			inst_mag_error_list.append(inst_mag_error)
			airmass_list.append(airmass)
			jd_list.append(jd)

		Plotter.plot_lightcurve(item, jd_list, inst_mag_list, inst_mag_error_list)

		yfit, slope, intercept, delta_slope, delta_intercept = unweighted_fit(airmass_list, inst_mag_list)
		Plotter.plot_extinction(item, airmass_list, inst_mag_list, yfit)

		return slope, intercept, delta_slope, delta_intercept