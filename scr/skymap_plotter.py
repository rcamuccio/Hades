import numpy as np
import os
import time

from conf import Configuration
from lib.boethos import Boethos
from lib.grapheus import Grapheus

# write the preamble
os.system('clear')
script_preamble = Boethos.get_script_preamble('skymap_plotter')
print(script_preamble)
line_preamble = Boethos.get_line_preamble('script', 'skymap_plotter')

# start the clock
st = time.time()

# set up plouton
Boethos.setup_plouton()
plouton_directory = Configuration.PLOUTON_DIRECTORY
skymap_directory = plouton_directory + 'skymap/'

# generate a blank survey map
print(line_preamble, 'Drawing default skymap')
default_skymap_path = skymap_directory + 'default' + Configuration.IMAGE_EXTENSION
Grapheus.draw_default_skymap(ra_deg=None, dec_deg=None, save_path=default_skymap_path, save_figure=Configuration.SAVE_FIGURE)

# generate a survey skymap
print(line_preamble, 'Drawing survey skymap')
survey_skymap_path = skymap_directory + 'survey' + Configuration.IMAGE_EXTENSION
Grapheus.draw_survey_skymap(save_path=survey_skymap_path, save_figure=Configuration.SAVE_FIGURE)

# generate a commissioning skymap
print(line_preamble, 'Drawing commissioning skymap')
commissioning_skymap_path = skymap_directory + 'commissioning' + Configuration.IMAGE_EXTENSION
Grapheus.draw_commissioning_skymap(save_path=commissioning_skymap_path, save_figure=Configuration.SAVE_FIGURE)

# generate a TNS skymap
print(line_preamble, 'Drawing TNS skymap')
tns_file_path = skymap_directory + 'tns_search.csv'
tns_skymap_path = skymap_directory + 'tns' + Configuration.IMAGE_EXTENSION
Grapheus.draw_tns_skymap(file_path=tns_file_path, save_path=tns_skymap_path, save_figure=Configuration.SAVE_FIGURE)

# write the postamble
script_postamble = Boethos.get_script_postamble('skymap_plotter')
fn = time.time()
dt = np.around(fn - st, decimals=2)
print(script_postamble, '[' + str(dt) + ' s]')