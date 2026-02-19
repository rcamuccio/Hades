from config import Configuration
from lib.clean import Clean
from lib.differencing import Diff
from lib.photometry import Photometry
from lib.master import Master
from lib.utilities import Utils
import os

if __name__ == '__main__':

	os.system('clear')
	Utils.log('Running main pipeline.', 'info')

	Utils.config_output_dir(create=True)
	Utils.config_data_dir(create=True)

	if Configuration.CLEAN_SKIP == 'N':
		Clean.clean_images()
	else:
		Utils.log('Skipping image cleaning.', 'info')

	if Configuration.MASTER_SKIP == 'N':
		master, star_list = Master.pull_master()
	else:
		Utils.log('Skipping master frame generation.','info')

	if Configuration.DIFFERENCE_SKIP == 'N':
		Diff.difference_images()
	else:
		Utils.log('Skipping image differencing.', 'info')

	if Configuration.PHOTOMETRY_SKIP == 'N':
		Photometry.generate_flux_files()
	else:
		Utils.log('Skipping photometry.', 'info')

	if Configuration.LIGHTCURVE_SKIP == 'N':
		Photometry.mk_raw_lightcurves()
	else:
		Utils.log('Skipping light curve generation.', 'info')

	Utils.log('Pipeline terminated.', 'info')