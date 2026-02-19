from config import Configuration
from lib.utilities import Utils
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

class Query:

	@staticmethod
	def query_gaia_cone(ra, de, radius):

		Gaia.MAIN_GAIA_TABLE = 'gaiadr3.gaia_source'
		Gaia.ROW_LIMIT = -1

		Utils.log('Submitting Gaia cone search of ' + str(radius) + ' deg at position (RA = ' + str(ra) + ' deg, DE = ' + str(de) + ' deg).', 'info')

		search_coord = SkyCoord(str(ra), str(de), unit=(u.deg, u.deg), frame='icrs')
		search_radius = u.Quantity(radius, u.deg)

		search = Gaia.cone_search_async(search_coord, radius=search_radius)

		table = search.get_results()

		Utils.log('Gaia table returned with ' + str(len(table)) + ' sources.', 'info')
		Utils.log('Writing Gaia table to file.', 'info')
		table.write('table-query.cat', format='ascii.fixed_width', overwrite=True)

		return table