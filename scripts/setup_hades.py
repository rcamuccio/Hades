from config import Configuration
from libraries.utils import Utils
import os

# set up tartarus directory
if not os.path.exists(Configuration.TARTARUS_DIR):
	os.mkdir(Configuration.TARTARUS_DIR)

# set up logs directory
logs_dir = os.path.join(Configuration.TARTARUS_DIR, 'logs/')
if not os.path.exists(logs_dir):
	os.mkdir(logs_dir)

# set up alerts directory
alerts_dir = os.path.join(Configuration.TARTARUS_DIR, 'alerts/')
if not os.path.exists(alerts_dir):
	os.mkdir(alerts_dir)

alerts_gw_dir = os.path.join(alerts_dir, 'gw/')
if not os.path.exists(alerts_gw_dir):
	os.mkdir(alerts_gw_dir)

alerts_gw_json_dir = os.path.join(alerts_gw_dir, 'json/')
if not os.path.exists(alerts_gw_json_dir):
	os.mkdir(alerts_gw_json_dir)

alerts_gw_moc_dir = os.path.join(alerts_gw_dir, 'moc/')
if not os.path.exists(alerts_gw_moc_dir):
	os.mkdir(alerts_gw_moc_dir)

alerts_gw_param_dir = os.path.join(alerts_gw_dir, 'param/')
if not os.path.exists(alerts_gw_param_dir):
	os.mkdir(alerts_gw_param_dir)

alerts_gw_skymap_dir = os.path.join(alerts_gw_dir, 'skymap/')
if not os.path.exists(alerts_gw_skymap_dir):
	os.mkdir(alerts_gw_skymap_dir)

# set up analysis directory
analysis_dir = os.path.join(Configuration.TARTARUS_DIR, 'analysis/')
if not os.path.exists(analysis_dir):
	os.mkdir(analysis_dir)

analysis_survey_dir = os.path.join(analysis_dir, 'survey/')
if not os.path.exists(analysis_survey_dir):
	os.mkdir(analysis_survey_dir)

analysis_gw_dir = os.path.join(analysis_dir, 'gw/')
if not os.path.exists(analysis_gw_dir):
	os.mkdir(analysis_gw_dir)

analysis_gw_fields_dir = os.path.join(analysis_gw_dir, 'fields/')
if not os.path.exists(analysis_gw_fields_dir):
	os.mkdir(analysis_gw_fields_dir)

analysis_gw_galaxies_dir = os.path.join(analysis_gw_dir, 'galaxies/')
if not os.path.exists(analysis_gw_galaxies_dir):
	os.mkdir(analysis_gw_galaxies_dir)

# set up catalogs directory
catalogs_dir = os.path.join(Configuration.TARTARUS_DIR, 'catalogs/')
if not os.path.exists(catalogs_dir):
	os.mkdir(catalogs_dir)

Utils.log('HADES setup complete.', 'info')