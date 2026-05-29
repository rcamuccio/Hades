from config import Configuration
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import os

from astropy import log
log.setLevel('WARNING')

class Query:

	@staticmethod
	def aavso_vsx_cone(query_ra, query_dec, query_radius, table_path):
		'''This function queries a sky region and returns a catalog of AAVSO VSX sources.

		:parameter query_ra - The query right ascension [deg]
		:parameter query_dec - The query declination [deg]
		:parameter query_radius - The query radius [deg]
		:parameter table_path - The save path for the table

		:return query_table - The table of AAVSO VSX sources
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
	def gaia_cone(query_ra, query_dec, query_radius, table_path):
		'''This function queries a sky region and returns a catalog of Gaia DR3 sources.

		:parameter query_ra - The query right ascension [deg]
		:parameter query_dec - The query declination [deg]
		:parameter query_radius - The query radius [deg]
		:parameter table_path - The save path for the table

		:return query_table - The table of Gaia DR3 sources
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

		:parameter query_ra - The query right ascension [deg]
		:parameter query_dec - The query declination [deg]
		:parameter query_radius - The query radius [deg]
		:parameter table_path - The save path for the table

		:return query_table - The table of GLADE+ sources
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