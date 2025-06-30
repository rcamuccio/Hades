from config import Configuration
from libraries.clean import Clean
from libraries.difference import Difference
from libraries.lightcurves import Lightcurves
from libraries.master import Master
from libraries.utils import Utils

if Configuration.CLEAN_SKIP == 'N':
	Clean.clean_images(clip_image='Y', subtract_bias='Y', subtract_dark='Y', divide_flat='Y', subtract_sky='Y', plate_solve='Y')
else:
	Utils.log('Skipping image cleaning.', 'info')

if Configuration.MASTER_SKIP == 'N':
	master, star_list = Master.pull_master()
else:
	Utils.log('Skipping master frame generation.', 'info')

if Configuration.DIFFERENCE_SKIP == 'N':
	Difference.difference_images(star_list)
else:
	Utils.log('Skipping image differencing.', 'info')

if Configuration.PHOTOMETRY_SKIP == 'N':
	Lightcurves.generate_flux_files(star_list)
else:
	Utils.log('Skipping photometry.', 'info')

if Configuration.LIGHTCURVE_SKIP == 'N':
	Lightcurves.mk_raw_lightcurves(star_list)
else:
	Utils.log('Skipping making raw light curves.', 'info')

Utils.log('All done! See ya later, alligator.', 'info')