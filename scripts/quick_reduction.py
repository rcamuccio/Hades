from config import Configuration
from libraries.reducer import Reducer
from libraries.utils import Utils

import os
import time

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

os.system('clear')
print('Running quick reduction')
start_time = time.time()

dark_dir = os.path.join(Configuration.WORKING_DIR, 'dark')
dark_flat_dir = os.path.join(dark_dir, 'flat')
dark_light_dir = os.path.join(dark_dir, 'light')
flat_dir = os.path.join(Configuration.WORKING_DIR, 'flat')
light_dir = os.path.join(Configuration.WORKING_DIR, Configuration.OBJECT)

raw_dir, cal_dir, wcs_dir, align_dir = Utils.create_directories(light_dir)

Reducer.make_dark(dark_flat_dir)
Reducer.make_dark(dark_light_dir)
Reducer.make_flat(flat_dir, dark_flat_dir)
Reducer.reduce_objects(light_dir, flat_dir, dark_light_dir, bkg_method=Configuration.BKG_METHOD)
Reducer.solve_plates(light_dir)
Reducer.align_frames(light_dir)
Reducer.make_stack(light_dir)

end_time = time.time()
total_time = end_time - start_time
print()
print('Script ended in', '%.1f' % total_time, 'seconds')