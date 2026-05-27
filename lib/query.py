from config import Configuration
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.gaia import Gaia
import os

class Query:

	@staticmethod
	def gaia_cone(ra, dec, query_radius, table_path):

		if not os.path.exists(table_path):
			print('Submitting a ' + str(query_radius) + '-deg Gaia cone search at position (' + str(ra) + ', ' + str(dec) + ')')

			Gaia.MAIN_GAIA_TABLE = 'gaiadr3.gaia_source'
			Gaia.ROW_LIMIT = Configuration.ROW_LIMIT

			search_coo = SkyCoord(str(ra), str(dec), unit=(u.deg, u.deg), frame='icrs')
			search_rad = u.Quantity(query_radius, u.deg)
			search = Gaia.cone_search_async(search_coo, radius=search_rad)

			query_table = search.get_results()
			query_table.write(table_path, format=Configuration.TABLE_FORMAT)

		else:
			print('Reading existing Gaia cone search table')
			query_table = Table.read(table_path, format=Configuration.TABLE_FORMAT)

		return query_table