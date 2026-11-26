'''
Ζητητής - The Inquirer
'''

import os
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier

from conf import Configuration

from astropy import log
log.setLevel('WARNING')

class Zetetes:

	@staticmethod
	def match_catalog(source_table, query_table, match_table_path, query_type):
		'''This function matches an extracted source catalog with a queried source catalog.

		:p source_table [astropy.table] - the extracted source catalog
		:p query_table [astropy.table] - the queried source catalog
		:p match_table_path [str] - the matched catalog path
		:p query_type [str] - a toggle for which query is used [aavso, gaia]

		:r match_table [astropy.table] - the matched source catalog
		'''

		line_preamble = '\t\033[1m[Ζητητής | match_catalog]\033[0m :'

		# match the tables
		if not os.path.exists(match_table_path):
			# get the queried source positions
			if query_type == 'aavso':
				print(line_preamble, 'Matching AAVSO-VSX catalog')
				query_pos = SkyCoord(query_table['RAJ2000'], query_table['DEJ2000'], unit='deg')
			elif query_type == 'gaia':
				print(line_preamble, 'Matching Gaia DR3 catalog')
				query_pos = SkyCoord(query_table['ra'], query_table['dec'], unit='deg')

			# get the extracted source positions
			source_pos = SkyCoord(source_table['ra'], source_table[
				'dec'], unit='deg')

			# find the nearest catalog matches for the extracted sources
			idx, d2d, d3d = query_pos.match_to_catalog_sky(source_pos)

			# define the match tolerance
			t = d2d < Configuration.MATCH_TOLERANCE * u.arcsec

			# add the match tolerance flags and matched indices
			query_table['sep_flag'] = t
			query_table['idx'] = idx

			# generate the matched table
			dt = query_table.dtype
			names, dtypes = zip(*dt.descr)
			match_table = Table(names=names, dtype=dtypes)
			for row, sep_flag in enumerate(query_table['sep_flag']):
				if sep_flag:
					match_table.add_row(query_table[row].as_void())

			# save the matched table
			match_table.write(match_table_path, format=Configuration.TABLE_FORMAT)

		# read the existing matched table
		else:
			if query_type =='aavso':
				print(line_preamble, 'Reading matched AAVSO-VSX catalog')
				match_table = Table.read(match_table_path, format=Configuration.TABLE_FORMAT)
			elif query_type == 'gaia':
				print(line_preamble, 'Reading matched Gaia DR3 catalog')
				match_table = Table.read(match_table_path, format=Configuration.TABLE_FORMAT)

		return match_table

	@staticmethod
	def parse_tns_file(file_path):

		tns_file = open(file_path, 'r')

		for line in tns_file:

			line = line.split(',')

			if line[0] != '"ID"':
				ra = line[2].strip('"')
				dec = line[3].strip('"')
				type = line[4].strip('"')

				coords = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')

				ra_deg = coords.ra.deg
				dec_deg = coords.dec.deg

				ra_shift = np.remainder(ra_deg + 180, 360) - 180

				ra_rad = np.radians(ra_shift)
				dec_rad = np.radians(dec_deg)

				slsni_ra_list = []
				slsni_dec_list = []

				snia_ra_list = []
				snia_dec_list = []

				snia91bg_ra_list = []
				snia91bg_dec_list = []

				snia91t_ra_list = []
				snia91t_dec_list = []

				sniacsm_ra_list = []
				sniacsm_dec_list = []

				sniax02cx_ra_list = []
				sniax02cx_dec_list = []
				
				snib_ra_list = []
				snib_dec_list = []
				
				snibc_ra_list = []
				snibc_dec_list =-[]
				
				snic_list = []
				
				snicbl_list = []

				
				snii_list = []
				
				sniib_list = []
				
				sniin_list = []
				
				sniip_list = []




	@staticmethod
	def query_aavso_vsx_cone(query_ra, query_dec, query_radius, table_path):
		'''This function queries a sky region and returns a catalog of AAVSO VSX sources.

		:p query_ra [float] - the query right ascension [deg]
		:p query_dec [float] - the query declination [deg]
		:p query_radius [float] - the query radius [deg]
		:p table_path [str] - the table save path

		:r query_table [astropy.table] - the AAVSO VSX source table
		'''

		if not os.path.exists(table_path):
			print('Submitting a ' + str(query_radius) + '-deg AAVSO VSX cone search at position (' + str(query_ra) + ', ' + str(query_dec) + ')')

			# configure the query
			vizier = Vizier()
			vizier.ROW_LIMIT = Configuration.ROW_LIMIT

			# set up the query geometry
			search_coord = SkyCoord(str(query_ra), str(query_dec), unit=(u.deg, u.deg), frame='icrs')
			search_radius = u.Quantity(query_radius, u.deg)

			# run the query
			search = vizier.query_region(search_coord, radius=search_radius, catalog='B/vsx')
			query_table = search[0]

			# save the table to file
			query_table.write(table_path, format=Configuration.TABLE_FORMAT)

		else:
			print('Reading existing AAVSO VSX cone search table')
			query_table = Table.read(table_path, format=Configuration.TABLE_FORMAT)

		return query_table

	@staticmethod
	def query_gaia_dr3_cone(query_ra, query_dec, query_radius, table_path):
		'''This function queries a sky region and returns a catalog of Gaia DR3 sources.

		:p query_ra [float] - the query right ascension [deg]
		:p query_dec [float] - the query declination [deg]
		:p query_radius [float] - the query radius [deg]
		:p table_path [str] - the table save path

		:r query_table [astropy.table] - the Gaia DR3 source catalog
		'''

		if not os.path.exists(table_path):
			print('Submitting a ' + str(query_radius) + '-deg Gaia cone search at position (' + str(query_ra) + ', ' + str(query_dec) + ')')

			# configure the query
			Gaia.MAIN_GAIA_TABLE = 'gaiadr3.gaia_source'
			Gaia.ROW_LIMIT = Configuration.ROW_LIMIT

			# set up the query geometry
			search_coord = SkyCoord(str(query_ra), str(query_dec), unit=(u.deg, u.deg), frame='icrs')
			search_radius = u.Quantity(query_radius, u.deg)

			# run the query
			search = Gaia.cone_search_async(search_coord, radius=search_radius)
			query_table = search.get_results()

			# save the table to file
			query_table.write(table_path, format=Configuration.TABLE_FORMAT)

		else:
			print('Reading existing Gaia cone search table')
			query_table = Table.read(table_path, format=Configuration.TABLE_FORMAT)

		return query_table

	@staticmethod
	def glade_plus_cone(query_ra, query_dec, query_radius, table_path):
		'''This function queries a sky region and returns a catalog of GLADE+ sources.

		:p query_ra [float] - the query right ascension [deg]
		:p query_dec [float] - the query declination [deg]
		:p query_radius [float] - the query radius [deg]
		:p table_path [str] - the table save path

		:r query_table [astropy.table] - the GLADE+ source catalog
		'''

		if not os.path.exists(table_path):
			print('Submitting a ' + str(query_radius) + '-deg GLADE+ cone search at position (' + str(query_ra) + ',' + str(query_dec) + ')')

			# configure the query
			vizier = Vizier()
			vizier.ROW_LIMIT = Configuration.ROW_LIMIT

			# set up the query geometry
			search_coord = SkyCoord(str(query_ra), str(query_dec), unit=(u.deg, u.deg), frame='icrs')
			search_radius = u.Quantity(query_radius, u.deg)

			# run the query
			search = vizier.query_region(search_coord, radius=search_radius, catalog='VII/291')
			query_table = search[0]

			# save the table to file
			query_table.write(table_path, format=Configuration.TABLE_FORMAT)

		else:
			print('Reading existing GLADE+ cone search table')
			query_table = Table.read(table_path, format=Configuration.TABLE_FORMAT)

		return query_table

	@staticmethod
	def transform_sky_to_pixel(frame_header, source_table, query_type):
		'''This function creates a list of pixel coordinate tuples (xp, yp) from a set of sky coordinate tuples (ra, dec) and a FITS frame header.

		:p frame_header [astropy.io.fits.header] - the image header
		:p source_table [astropy.table] - the source catalog
		:p query_type [str] - a toggle for which query is used [aavso, gaia]

		:r pix_coord_list [list][tuple] - the list of pixel coordinates
		'''

		# process the WCS information
		frame_wcs = WCS(frame_header)

		# read the sky coordinates
		sky_coord_list = []

		for src in source_table:
			if query_type == 'aavso':
				ra = src['RAJ2000']
				dec = src['DEJ2000']
			elif query_type == 'gaia':
				ra = src['ra']
				dec = src['dec']

			coord_tuple = (ra, dec)
			sky_coord_list.append(coord_tuple)

		sky_pos = SkyCoord(sky_coord_list, unit='deg')

		# convert the sky coordinates to pixel coordinates
		pix_pos = skycoord_to_pixel(sky_pos, wcs=frame_wcs)
		pix_coord_list = []

		for pxl in range(len(pix_pos[0])):
			xp = pix_pos[0][pxl]
			yp = pix_pos[1][pxl]

			coord_tuple = (xp, yp)
			pix_coord_list.append(coord_tuple)

		return pix_coord_list