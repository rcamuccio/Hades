from config import Configuration
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.gaia import Gaia
import os

class Query:

	@staticmethod
	def gaia_cone(ra, de, qry_rad, tbl_path):

		if not os.path.exists(tbl_path):
			print('Submitting Gaia cone search', qry_rad, ra, de)

			Gaia.MAIN_GAIA_TABLE = 'gaiadr3.gaia_source'
			Gaia.ROW_LIMIT = -1

			search_coo = SkyCoord(str(ra), str(de), unit=(u.deg, u.deg), frame='icrs')
			search_rad = u.Quantity(qry_rad, u.deg)
			search = Gaia.cone_search_async(search_coo, radius=search_rad)

			tbl = search.get_results()
			tbl.write(tbl_path, format=Configuration.TABLE_FORMAT)

		else:
			print('Reading existing Gaia cone search table')
			tbl = Table.read(tbl_path, format=Configuration.TABLE_FORMAT)

		return tbl