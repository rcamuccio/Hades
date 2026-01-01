from config import Configuration
from libraries.clean import Clean
from libraries.difference import BigDiff
from libraries.lightcurves import Lightcurves
from libraries.master import Master
from libraries.utils import Utils

if __name__ == '__main__':

    Utils.create_directories(Configuration.DIRECTORIES)

    if Configuration.CLEAN_SKIP == 'N':
        Clean.clean_images()
    else:
        Utils.log('Skipping image cleaning.', 'info')

    if Configuration.MASTER_SKIP == 'N':
        master, star_list = Master.pull_master()
    else:
        Utils.log('Skipping master frame generation.', 'info')

    if Configuration.DIFFERENCE_SKIP == 'N':
        BigDiff.difference_images()
    else:
        Utils.log('Skipping image differencing.', 'info')

    if Configuration.PHOTOMETRY_SKIP == 'N':
        Lightcurves.generate_flux_files()
    else:
        Utils.log('Skipping photometry.', 'info')

    if Configuration.LIGHTCURVE_SKIP == 'N':
        Lightcurves.mk_raw_lightcurves()
    else:
        Utils.log('Skipping light curve generation.', 'info')

    Utils.log('All done! See ya later, alligator.', 'info')