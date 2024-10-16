from config import Configuration
from libraries.photometer import *
from libraries.querier import *

from astropy.table import Table
import os
import time

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

os.system('clear')
print('Running limiting magnitude')
start_time = time.time()

photometer = Photometer()
querier = Querier()

stack = 'stack.fit'
stack_path = os.path.join(Configuration.WORKING_DIR, stack)

# --- Source extract field
source_path = os.path.join(Configuration.WORKING_DIR, 'source.cat')

if os.path.isfile(source_path):
	print('Reading extant source table')
	source_table = Table.read(source_path, format='ascii.tab')

else:
	source_table = photometer.extract_sky_sources(stack_path)
	source_table.write(source_path, format='ascii.tab', overwrite=True)

# --- Query field
if Configuration.CATALOG == 'gaia-cone':
	query_path = os.path.join(Configuration.WORKING_DIR, 'gaiadr3-cone.cat')

	if os.path.isfile(query_path):
		print('Reading extant query table (gaia-cone)')
		query_table = Table.read(query_path, format='ascii.tab')

	else:
		query_table = querier.query_gaia_cone(Configuration.FIELD_RA, Configuration.FIELD_DEC, Configuration.RAD_QUERY)
		query_table.write(query_path, format='ascii.tab', overwrite=True)

elif Configuration.CATALOG == 'gaia-square':
	pass

elif Configuration.CATALOG == 'ps1':
	query_path = os.path.join(Configuration.WORKING_DIR, 'ps1-cone.cat')

	if os.path.isfile(query_path):
		print('Reading extant query table (ps1-cone)')
		query_table = Table.read(query_path, format='ascii.tab')

	else:
		query_table = querier.query_ps1_cone(Configuration.FIELD_RA, Configuration.FIELD_DEC, Configuration.RAD_QUERY)
		query_table.write(query_path, format='ascii.tab', overwrite=True)

else:
	pass

# --- Match extracted and queried fields
match_path = os.path.join(Configuration.WORKING_DIR, 'new-match.cat')

if os.path.isfile(match_path):
	print('Reading extant matched catalog')
	match_table = Table.read(match_path, format='ascii.tab')

else:
	match_table = photometer.match_catalogs(source_table, query_table)
	match_table.write(match_path, format='ascii.tab', overwrite=True)

# --- Relative photometry on field
catalog_path = os.path.join(Configuration.WORKING_DIR, 'catalog.cat')

if os.path.isfile(catalog_path):
	print('Reading extant master catalog')
	catalog_table = Table.read(catalog_path, format='ascii.tab')

else:
	catalog_table = photometer.photometry_field(stack_path, match_table, survey='gaia')
	catalog_table.write(catalog_path, format='ascii.tab', overwrite=True)

end_time = time.time()
total_time = end_time - start_time
print()
print('Script ended in', '%.1f' % total_time, 'seconds')