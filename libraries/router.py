from config import Configuration
from libraries.alerts import Alerts
from libraries.priority import Priority
from libraries.utils import Utils

from astropy import units as u
from astropy.table import QTable, Table
from base64 import b64decode
from io import BytesIO
import astropy_healpix as ah
import json
import numpy as np
import os

class Router:

	@staticmethod
	def decode_classic_notice(value, verbose=False):

		contents = value.decode('utf-8').split('\n')

		keys = []
		items = []

		for content in contents:
			if len(content.split(': ')) > 1:
				key = content.split(': ')[0].strip()
				keys.append(key)

				item = content.split(': ')[1].strip()
				items.append(item)

				if verbose:
					print(key, item)

		params = dict(zip(keys, items))

		if verbose:
			print('----------')

		return params

	@staticmethod
	def route_alert(value, topic):

		Utils.log('Alert received from ' + topic + ' channel!', 'info')

		# process LVK kafka
		if topic == 'igwn.gwalert':
			Router.igwn_gwalert(value)

		# process IceCube kafka
		elif topic == 'gcn.notices.icecube.lvk_nu_track_search':
			Router.lvk_nu_track_search(value)

		# process Swift kafka
		elif topic == 'gcn.notices.swift.bat.guano':
			Router.swift_bat_guano(value)

		# process Fermi classic
		elif topic == 'gcn.classic.text.FERMI_GBM_ALERT':
			Router.process_fermi_classic_alert(value, topic)

		elif topic == 'gcn.classic.text.FERMI_GBM_FIN_POS':
			Router.process_fermi_classic_alert(value, topic)

		elif topic == 'gcn.classic.text.FERMI_GBM_FLT_POS':
			Router.process_fermi_classic_alert(value, topic)

		elif topic == 'gcn.classic.text.FERMI_GBM_GND_POS':
			Router.process_fermi_classic_alert(value, topic)

		elif topic == 'gcn.classic.text.FERMI_GBM_POS_TEST':
			Router.process_fermi_classic_alert(value, topic)

		elif topic == 'gcn.classic.text.FERMI_GBM_SUBTHRESH':
			Router.process_fermi_classic_alert(value, topic)

		elif topic == 'gcn.classic.text.FERMI_POINTDIR':
			Router.process_fermi_classic_alert(value, topic)

		# process Swift classic
		elif topic == 'gcn.classic.text.SWIFT_ACTUAL_POINTDIR':
			Router.process_swift_classic_alert(value, topic)

		elif topic == 'gcn.classic.text.SWIFT_BAT_GRB_POS_TEST':
			Router.process_swift_classic_alert(value, topic)

		elif topic == 'gcn.classic.text.SWIFT_POINTDIR':
			Router.process_swift_classic_alert(value, topic)

		# process unknown alert
		else:
			print(topic, type(topic))
			params = Router.decode_classic_notice(value, verbose=True)
			print(params)

	@staticmethod
	def process_fermi_classic_alert(value, topic):

		params = Router.decode_classic_notice(value)

		if topic == 'gcn.classic.text.FERMI_GBM_ALERT':

			''' G1) GBM_Alert Notice (r-t) starts the sequence of GBM messages.
			It occurs when the GBM instrument first triggers.
			It does NOT contain an RA,Dec location of a burst;
			only a date, time, and trigger criteria and trigger detection significance,
			and the algorithm used to make the detection.
			It is issued only once per trigger, (but may not be present for all triggers
			if there is TDRSS bit/frame sync-up delay -- that's the reason the Alerts exist,
			to start the TDRSS sync-up process so it's ready by the time
			the GBM_Pos Notice comes along).

			'''
			
			title = params['TITLE']
			notice_date = params['NOTICE_DATE']
			notice_type = params['NOTICE_TYPE']
			record_num = params['RECORD_NUM']
			trigger_num = params['TRIGGER_NUM']
			grb_date = params['GRB_DATE']
			grb_time = params['GRB_TIME']
			trigger_signif = params['TRIGGER_SIGNIF']
			trigger_dur = params['TRIGGER_DUR']
			e_range = params['E_RANGE']
			algorithm = params['ALGORITHM']
			detectors = params['DETECTORS']
			comments = params['COMMENTS']

			print('Title:', title)
			print('Notice Date:', notice_date)
			print('Notice Type:', notice_type)
			print('Record Num:', record_num)
			print('Trigger Num:', trigger_num)
			print('GRB Date:', grb_date)
			print('GRB Time:', grb_time)
			print('Trigger Signif:', trigger_signif)
			print('Trigger Dur:', trigger_dur)
			print('E Range:', e_range)
			print('Algorithm:', algorithm)
			print('Detectors:', detectors)
			print('Comments:', comments)
			print()

		elif topic == 'gcn.classic.text.FERMI_GBM_FIN_POS':

			''' G4) GBM_Final_Position Notice (delayed) And now humans are involved in the analysis.
			The humans can make better selections on the data used.
			There can be 0, 1, or possibly more instances of this Notice Type per trigger;
			They will be issued for the 10-20% brightest of the GRBs.
			(There may be others, less bright, that also get a Final notice,
			but they will have a larger location uncertainty -- see table below.)

			'''
			
			title = params['TITLE']
			notice_date = params['NOTICE_DATE']
			notice_type = params['NOTICE_TYPE']
			record_num = params['RECORD_NUM']
			trigger_num = params['TRIGGER_NUM']
			grb_ra = params['GRB_RA']
			grb_dec = params['GRB_DEC']
			grb_error = params['GRB_ERROR']
			grb_date = params['GRB_DATE']
			grb_time = params['GRB_TIME']
			grb_phi = params['GRB_PHI']
			grb_theta = params['GRB_THETA']
			e_range = params['E_RANGE']
			loc_algorithm = params['LOC_ALGORITHM']
			lc_url = params['LC_URL']
			loc_url = params['LOC_URL']
			sun_postn = params['SUN_POSTN']
			sun_dist = params['SUN_DIST']
			moon_postn = params['MOON_POSTN']
			moon_dist = params['MOON_DIST']
			moon_illum = params['MOON_ILLUM']
			gal_coords = params['GAL_COORDS']
			ecl_coords = params['ECL_COORDS']
			comments = params['COMMENTS']

			print('Title:', title)
			print('Notice Date:', notice_date)
			print('Notice Type:', notice_type)
			print('Record Num:', record_num)
			print('Trigger Num:', trigger_num)
			print('GRB RA:', grb_ra)
			print('GRB Dec:', grb_dec)
			print('GRB Error:', grb_error)
			print('GRB Date:', grb_date)
			print('GRB Time:', grb_time)
			print('GRB Phi:', grb_phi)
			print('GRB Theta:', grb_theta)
			print('E Range:', e_range)
			print('Loc Algorithm:', loc_algorithm)
			print('LC URL:', lc_url)
			print('Loc URL:', loc_url)
			print('Sun Postn:', sun_postn)
			print('Sun Dist:', sun_dist)
			print('Moon Postn:', moon_postn)
			print('Moon Dist:', moon_dist)
			print('Moon Illum:', moon_illum)
			print('Gal Coords:', gal_coords)
			print('Ecl Coords:', ecl_coords)
			print('Comments:', comments)
			print()

		elif topic == 'gcn.classic.text.FERMI_GBM_FLT_POS':

			''' G2) GBM_Flight_Position Notice (r-t) contains the RA,Dec location for the burst detected by GBM.
			The positions are calculated by the on-board Flight software.
			They are issued 1 to 5 times per burst.
			The Flight Position Notice comes second in the sequence of Notices on the triggers/bursts.
			Since the position is based on the least-possible amount of data in the processing
			(just a small initial portion of the burst's lightcurve),
			it has the lowest significance in the localization process. Even so,
			the uncertainty in the position is ~20 deg (1-sigma, radius) for the at-threshold bursts
			and ~10 deg for the bright bursts. (Both the 'at-threshold' and 'bright'
			refer to only the amount of photons in the time interval of the trigger-sampling interval
			of the burst lightcurve, not to the total burst duration.)
			It should be noted that even though this Notice type is called the
			GBM_FLT_POSITION, it contains detections of both GRBs and hard x-ray transients.
			There is both flight- and ground-software in place to correctly identify
			bursts from transients, but this identification is not 100% perfect.

			'''
			
			title = params['TITLE']
			notice_date = params['NOTICE_DATE']
			notice_type = params['NOTICE_TYPE']
			record_num = params['RECORD_NUM']
			trigger_num = params['TRIGGER_NUM']
			grb_ra = params['GRB_RA']
			grb_dec = params['GRB_DEC']
			grb_error = params['GRB_ERROR']
			grb_inten = params['GRB_INTEN']
			data_signif = params['DATA_SIGNIF']
			integ_time = params['INTEG_TIME']
			grb_date = params['GRB_DATE']
			grb_time = params['GRB_TIME']
			grb_phi = params['GRB_PHI']
			grb_theta = params['GRB_THETA']
			hard_ratio = params['HARD_RATIO']
			loc_algorithm = params['LOC_ALGORITHM']
			most_likely = params['MOST_LIKELY']
			secmost_likely = params['2nd_MOST_LIKELY']
			detectors = params['DETECTORS']
			lc_url = params['LC_URL']
			sun_postn = params['SUN_POSTN']
			sun_dist = params['SUN_DIST']
			moon_postn = params['MOON_POSTN']
			moon_dist = params['MOON_DIST']
			moon_illum = params['MOON_ILLUM']
			gal_coords = params['GAL_COORDS']
			ecl_coords = params['ECL_COORDS']
			comments = params['COMMENTS']

			print('Title:', title)
			print('Notice Date:', notice_date)
			print('Notice Type:', notice_type)
			print('Record Num:', record_num)
			print('Trigger Num:', trigger_num)
			print('GRB RA:', grb_ra)
			print('GRB Dec:', grb_dec)
			print('GRB Error:', grb_error)
			print('GRB Inten:', grb_inten)
			print('Data Signif:', data_signif)
			print('Integ Time:', integ_time)
			print('GRB Date:', grb_date)
			print('GRB Time:', grb_time)
			print('GRB Phi:', grb_phi)
			print('GRB Theta:', grb_theta)
			print('Hard Ratio:', hard_ratio)
			print('Loc Algorithm:', loc_algorithm)
			print('Most Likely:', most_likely)
			print('2nd Most Likely:', secmost_likely)
			print('Detectors:', detectors)
			print('LC URL:', lc_url)
			print('Sun Postn:', sun_postn)
			print('Sun Dist:', sun_dist)
			print('Moon Postn:', moon_postn)
			print('Moon Dist:', moon_dist)
			print('Moon Illum:', moon_illum)
			print('Gal Coords:', gal_coords)
			print('Ecl Coords:', ecl_coords)
			print('Comments:', comments)
			print()

		elif topic == 'gcn.classic.text.FERMI_GBM_GND_POS':

			''' G3) GBM_Ground_Position Notice (r-t) contains the RA,Dec location for the burst detected by GBM.
			The positions are calculated by automated ground software as soon as they are received on the ground.
			More sophisticated algorithms can be applied to the data to improve the location accuracy.
			Also, more data from the on-going progresion of the burst lightcurve is used in this calculation.
			There can be 0, 1, or more instances of this Notice Type per trigger.

			'''
			
			title = params['TITLE']
			notice_date = params['NOTICE_DATE']
			notice_type = params['NOTICE_TYPE']
			record_num = params['RECORD_NUM']
			trigger_num = params['TRIGGER_NUM']
			grb_ra = params['GRB_RA']
			grb_dec = params['GRB_DEC']
			grb_error = params['GRB_ERROR']
			data_signif = params['DATA_SIGNIF']
			data_interval = params['DATA_INTERVAL']
			grb_date = params['GRB_DATE']
			grb_time = params['GRB_TIME']
			grb_phi = params['GRB_PHI']
			grb_theta = params['GRB_THETA']
			e_range = params['E_RANGE']
			loc_algorithm = params['LOC_ALGORITHM']
			lc_url = params['LC_URL']
			sun_postn = params['SUN_POSTN']
			sun_dist = params['SUN_DIST']
			moon_postn = params['MOON_POSTN']
			moon_dist = params['MOON_DIST']
			moon_illum = params['MOON_ILLUM']
			gal_coords = params['GAL_COORDS']
			ecl_coords = params['ECL_COORDS']
			comments = params['COMMENTS']

			print('Title:', title)
			print('Notice Date:', notice_date)
			print('Notice Type:', notice_type)
			print('Record Num:', record_num)
			print('Trigger Num:', trigger_num)
			print('GRB RA:', grb_ra)
			print('GRB Dec:', grb_dec)
			print('GRB Error:', grb_error)
			print('Data Signif:', data_signif)
			print('Data Interval:', data_interval)
			print('GRB Date:', grb_date)
			print('GRB Time:', grb_time)
			print('GRB Phi:', grb_phi)
			print('GRB Theta:', grb_theta)
			print('E Range:', e_range)
			print('Loc Algorithm:', loc_algorithm)
			print('LC URL:', lc_url)
			print('Sun Postn:', sun_postn)
			print('Sun Dist:', sun_dist)
			print('Moon Postn:', moon_postn)
			print('Moon Dist:', moon_dist)
			print('Moon Illum:', moon_illum)
			print('Gal Coords:', gal_coords)
			print('Ecl coords:', ecl_coords)
			print('Comments:', comments)
			print()

		elif topic == 'gcn.classic.text.FERMI_GBM_POS_TEST':

			title = params['TITLE']
			notice_date = params['NOTICE_DATE']
			notice_type = params['NOTICE_TYPE']
			record_num = params['RECORD_NUM']
			trigger_num = params['TRIGGER_NUM']
			grb_ra = params['GRB_RA']
			grb_dec = params['GRB_DEC']
			grb_error = params['GRB_ERROR']
			grb_inten = params['GRB_INTEN']
			data_signif = params['DATA_SIGNIF']
			integ_time = params['INTEG_TIME']
			grb_date = params['GRB_DATE']
			grb_time = params['GRB_TIME']
			grb_phi = params['GRB_PHI']
			grb_theta = params['GRB_THETA']
			data_time_scale = params['DATA_TIME_SCALE']
			hard_ratio = params['HARD_RATIO']
			loc_algorithm = params['LOC_ALGORITHM']
			most_likely = params['MOST_LIKELY']
			secmost_likely = params['2nd_MOST_LIKELY']
			detectors = params['DETECTORS']
			sun_postn = params['SUN_POSTN']
			sun_dist = params['SUN_DIST']
			moon_postn = params['MOON_POSTN']
			moon_dist = params['MOON_DIST']
			moon_illum = params['MOON_ILLUM']
			gal_coords = params['GAL_COORDS']
			ecl_coords = params['ECL_COORDS']
			comments = params['COMMENTS']

			print('Title:', title)
			print('Notice Date:', notice_date)
			print('Notice Type:', notice_type)
			print('Record Num:', record_num)
			print('Trigger Num:', trigger_num)
			print('GRB RA:', grb_ra)
			print('GRB Dec:', grb_dec)
			print('GRB Error:', grb_error)
			print('GRB Inten:', grb_inten)
			print('Data Signif:', data_signif)
			print('Integ Time:', integ_time)
			print('GRB Date:', grb_date)
			print('GRB Time:', grb_time)
			print('GRB Phi:', grb_phi)
			print('GRB Theta:', grb_theta)
			print('Data Timescale:', data_time_scale)
			print('Hard Ratio:', hard_ratio)
			print('Loc Algorithm:', loc_algorithm)
			print('Most Likely:', most_likely)
			print('2nd Most Likely:', secmost_likely)
			print('Detectors:', detectors)
			print('Sun Postn:', sun_postn)
			print('Sun Dist:', sun_dist)
			print('Moon Postn:', moon_postn)
			print('Moon Dist:', moon_dist)
			print('Moon Illum:', moon_illum)
			print('Gal Coords:', gal_coords)
			print('Ecl Coords:', ecl_coords)
			print('Comments:', comments)
			print()

		elif topic == 'gcn.classic.text.FERMI_GBM_SUBTHRESH':

			''' G4) GBM_Subthreshold Notice (delayed) This is a separate stream of transients than described above for the Alert to the Final.
			These are produced from a ground pipeline procesing of the data
			looking for transients that were below the on-board s/w trigger threshold level.
			These "subthreshold" events are typtically used for coincidence searchs with other data streams.

			'''

			title = params['TITLE']
			notice_date = params['NOTICE_DATE']
			notice_type = params['NOTICE_TYPE']
			trans_num = params['TRANS_NUM']
			full_id_num = params['FULL_ID_NUM']
			trans_ra = params['TRANS_RA']
			trans_dec = params['TRANS_DEC']
			trans_error = params['TRANS_ERROR']
			trans_duration = params['TRANS_DURATION']
			trans_date = params['TRANS_DATE']
			trans_time = params['TRANS_TIME']
			earth_angle = params['EARTH_ANGLE']
			spectral_class = params['SPECTRAL_CLASS']
			type_class = params['TYPE_CLASS']
			reliability = params['RELIABILITY']
			healpix_url = params['HEALPIX_URL']
			map_url = params['MAP_URL']
			lc_url = params['LC_URL']
			sun_postn = params['SUN_POSTN']
			sun_dist = params['SUN_DIST']
			moon_postn = params['MOON_POSTN']
			moon_dist = params['MOON_DIST']
			moon_illum = params['MOON_ILLUM']
			gal_coords = params['GAL_COORDS']
			ecl_coords = params['ECL_COORDS']
			comments = params['COMMENTS']

			print('Title:', title)
			print('Notice Date:', notice_date)
			print('Notice Type:', notice_type)
			print('Trans Num:', trans_num)
			print('Full ID Num:', full_id_num)
			print('Trans RA:', trans_ra)
			print('Trans Dec:', trans_dec)
			print('Trans Error:', trans_error)
			print('Trans Duration:', trans_duration)
			print('Trans Date:', trans_date)
			print('Trans Time:', trans_time)
			print('Earth Angle:', earth_angle)
			print('Spectral Class:', spectral_class)
			print('Type Class:', type_class)
			print('Reliability:', reliability)
			print('HEALPIX URL:', healpix_url)
			print('Map URL:', map_url)
			print('LC URL:', lc_url)
			print('Sun Postn:', sun_postn)
			print('Sun Dist:', sun_dist)
			print('Moon Postn:', moon_postn)
			print('Moon Dist:', moon_dist)
			print('Moon Illum:', moon_illum)
			print('Gal Coords:', gal_coords)
			print('Ecl coords:', ecl_coords)
			print('Comments:', comments)
			print()

		elif topic == 'gcn.classic.text.FERMI_POINTDIR':

			title = params['TITLE']
			notice_date = params['NOTICE_DATE']
			notice_type = params['NOTICE_TYPE']
			curr_point_ra = params['CURR_POINT_RA']
			curr_point_dec = params['CURR_POINT_DEC']
			curr_date = params['CURR_DATE']
			curr_time = params['CURR_TIME']
			delta_time = params['DELTA_TIME']
			sun_postn = params['SUN_POSTN']
			sun_dist = params['SUN_DIST']
			moon_postn = params['MOON_POSTN']
			moon_dist = params['MOON_DIST']
			moon_illum = params['MOON_ILLUM']
			gal_coords = params['GAL_COORDS']
			ecl_coords = params['ECL_COORDS']
			future_ra_dec = params['FUTURE_RA_DEC']
			comments = params['COMMENTS']

			print('Title:', title)
			print('Notice Date:', notice_date)
			print('Notice Type:', notice_type)
			print('Current Point RA:', curr_point_ra)
			print('Current Point Dec:', curr_point_dec)
			print('Current Date:', curr_date)
			print('Current Time:', curr_time)
			print('Delta Time:', delta_time)
			print('Sun Postn:', sun_postn)
			print('Sun Dist:', sun_dist)
			print('Moon Postn:', moon_postn)
			print('Moon Dist:', moon_dist)
			print('Moon Illum:', moon_illum)
			print('Gal coords:', gal_coords)
			print('Ecl coords:', ecl_coords)
			print('Future RA/Dec:', future_ra_dec)
			print('Comments:', comments)
			print()

	@staticmethod
	def process_swift_classic_alert(value, topic):

		params = Router.decode_classic_notice(value)

		if topic == 'gcn.classic.text.SWIFT_ACTUAL_POINTDIR':

			title = params['TITLE']
			notice_date = params['NOTICE_DATE']
			notice_type = params['NOTICE_TYPE']
			curr_point_ra = params['CURR_POINT_RA']
			curr_point_dec = params['CURR_POINT_DEC']
			curr_point_roll = params['CURR_POINT_ROLL']
			slew_time = params['SLEW_TIME']
			slew_date = params['SLEW_DATE']
			tgt_num = params['TGT_NUM']
			sun_postn = params['SUN_POSTN']
			sun_dist = params['SUN_DIST']
			moon_postn = params['MOON_POSTN']
			moon_dist = params['MOON_DIST']
			moon_illum = params['MOON_ILLUM']
			gal_coords = params['GAL_COORDS']
			ecl_coords = params['ECL_COORDS']
			comments = params['COMMENTS']

			print('Title:', title)
			print('Notice Date:', notice_date)
			print('Notice Type:', notice_type)
			print('Current Point RA:', curr_point_ra)
			print('Current Point Dec:', curr_point_dec)
			print('Current Point Roll:', curr_point_roll)
			print('Slew Time:', slew_time)
			print('Slew Date:', slew_date)
			print('TGT Num:', tgt_num)
			print('Sun Postn:', sun_postn)
			print('Sun Dist:', sun_dist)
			print('Moon Postn:', moon_postn)
			print('Moon Dist:', moon_dist)
			print('Moon Illum:', moon_illum)
			print('Gal coords:', gal_coords)
			print('Ecl coords:', ecl_coords)
			print('Comments:', comments)
			print()

		elif topic == 'gcn.classic.text.SWIFT_BAT_GRB_POS_TEST':

			title = params['TITLE']
			notice_date = params['NOTICE_DATE']
			notice_type = params['NOTICE_TYPE']
			trigger_num = params['TRIGGER_NUM']
			grb_ra = params['GRB_RA']
			grb_dec = params['GRB_DEC']
			grb_error = params['GRB_ERROR']
			grb_inten = params['GRB_INTEN']
			trigger_dur = params['TRIGGER_DUR']
			trigger_index = params['TRIGGER_INDEX']
			bkg_inten = params['BKG_INTEN']
			bkg_time = params['BKG_TIME']
			bkg_dur = params['BKG_DUR']
			grb_date = params['GRB_DATE']
			grb_time = params['GRB_TIME']
			grb_phi = params['GRB_PHI']
			grb_theta = params['GRB_THETA']
			soln_status = params['SOLN_STATUS']
			rate_signif = params['RATE_SIGNIF']
			image_signif = params['IMAGE_SIGNIF']
			merit_params = params['MERIT_PARAMS']
			sun_postn = params['SUN_POSTN']
			sun_dist = params['SUN_DIST']
			moon_postn = params['MOON_POSTN']
			moon_dist = params['MOON_DIST']
			moon_illum = params['MOON_ILLUM']
			gal_coords = params['GAL_COORDS']
			ecl_coords = params['ECL_COORDS']
			comments = params['COMMENTS']

			print('Title:', title)
			print('Notice Date:', notice_date)
			print('Notice Type:', notice_type)
			print('Trigger Num:', trigger_num)
			print('GRB RA:', grb_ra)
			print('GRB Dec:', grb_dec)
			print('GRB Error:', grb_error)
			print('GRB Inten:', grb_inten)
			print('Trigger Dur:', trigger_dur)
			print('Trigger Index:', trigger_index)
			print('Bkg Inten:', bkg_inten)
			print('Bkg Time:', bkg_time)
			print('Bkg Dur:', bkg_dur)
			print('GRB Date:', grb_date)
			print('GRB Time:', grb_time)
			print('GBR Phi:', grb_phi)
			print('GRB Theta', grb_theta)
			print('Soln Status:', soln_status)
			print('Rate Signif:', rate_signif)
			print('Image Signif:', image_signif)
			print('Merit Params:', merit_params)
			print('Sun Postn:', sun_postn)
			print('Sun Dist:', sun_dist)
			print('Moon Postn:', moon_postn)
			print('Moon Dist:', moon_dist)
			print('Moon Illum:', moon_illum)
			print('Gal coords:', gal_coords)
			print('Ecl coords:', ecl_coords)
			print('Comments:', comments)
			print()

		elif topic == 'gcn.classic.text.SWIFT_POINTDIR':
			
			title = params['TITLE']
			notice_date = params['NOTICE_DATE']
			notice_type = params['NOTICE_TYPE']
			next_point_ra = params['NEXT_POINT_RA']
			next_point_dec = params['NEXT_POINT_DEC']
			next_point_roll = params['NEXT_POINT_ROLL']
			slew_time = params['SLEW_TIME']
			slew_date = params['SLEW_DATE']
			obs_time = params['OBS_TIME']
			tgt_name = params['TGT_NAME']
			tgt_num = params['TGT_NUM']
			merit = params['MERIT']
			inst_modes = params['INST_MODES']
			sun_postn = params['SUN_POSTN']
			sun_dist = params['SUN_DIST']
			moon_postn = params['MOON_POSTN']
			moon_dist = params['MOON_DIST']
			moon_illum = params['MOON_ILLUM']
			gal_coords = params['GAL_COORDS']
			ecl_coords = params['ECL_COORDS']
			comments = params['COMMENTS']

			print('Title:', title)
			print('Notice Date:', notice_date)
			print('Notice Type:', notice_type)
			print('Next Point RA:', next_point_ra)
			print('Next Point Dec:', next_point_dec)
			print('Next Point Roll:', next_point_roll)
			print('Slew Time:', slew_time)
			print('Slew Date:', slew_date)
			print('Obs Time:', obs_time)
			print('TGT Name:', tgt_name)
			print('TGT Num:', tgt_num)
			print('Merit:', merit)
			print('Inst Modes:', inst_modes)
			print('Sun Postn:', sun_postn)
			print('Sun Dist:', sun_dist)
			print('Moon Postn:', moon_postn)
			print('Moon Dist:', moon_dist)
			print('Moon Illum:', moon_illum)
			print('Gal coords:', gal_coords)
			print('Ecl coords:', ecl_coords)
			print('Comments:', comments)
			print()

	@staticmethod
	def lvk_nu_track_search(value):

		record = json.loads(value)
			
		alert_type = record['type']
		reference = record['reference']
		ref_id = record['ref_ID']
		alert_datetime = record['alert_datetime']
		trigger_time = record['trigger_time']
		observation_start = record['observation_start']
		observation_stop = record['observation_stop']
		observation_livetime = record['observation_livetime']
		pval_generic = record['pval_generic']
		pval_bayesian = record['pval_bayesian']
		n_events_coincident = record['n_events_coincident']
		coincident_events = record['coincident_events']
		most_probable_direction = record['most_probable_direction']
		flux_sensitivity = record['neutrino_flux_sensitivity_range']['flux_sensitivity']
		sensitive_energy_range = record['neutrino_flux_sensitivity_range']['sensitive_energy_range']

		print('Alert Type:', alert_type, type(alert_type))
		print('Reference:', reference, type(reference))
		print('Ref ID:', ref_id, type(ref_id))
		print('Alert Datetime:', alert_datetime, type(alert_datetime))
		print('Trigger Time:', trigger_time, type(trigger_time))
		print('Observation Start:', observation_start, type(observation_start))
		print('Observation Stop:', observation_stop, type(observation_stop))
		print('Observation Livetime:', observation_livetime, type(observation_livetime))
		print('pval Generic:', pval_generic, type(pval_generic))
		print('pval Bayesian:', pval_bayesian, type(pval_bayesian))
		print('n Events Coincident:', n_events_coincident, type(n_events_coincident))
		print('Coincident Events:', coincident_events, type(coincident_events), len(coincident_events))
		print('Most Probable Direction:', most_probable_direction, type(most_probable_direction))
		print('Flux Sensitivity:', flux_sensitivity, type(flux_sensitivity))
		print('Sensitive Energy Range:', sensitive_energy_range, type(sensitive_energy_range))
		print()

	@staticmethod
	def swift_bat_guano(value):

		record = json.loads(value)

		mission = record['mission']
		instrument = record['instrument']
		messenger = record['messenger']
		record_number = record['record_number']
		alert_datetime = record['alert_datetime']
		alert_tense = record['alert_tense']
		alert_type = record['alert_type']
		trigger_time = record['trigger_time']
		follow_up_event = record['follow_up_event']
		follow_up_type = record['follow_up_type']
		data_archive_page = record['data_archive_page']
		alert_id = record['id']

		Utils.log('Event BAT-GUANO/' + follow_up_event + ' generated at time ' + str(alert_datetime) + '.', 'info')
		Utils.log('Data Archive Page: ' + str(data_archive_page) + '.', 'info')

		# 1 - guano.example.json
		if record_number == 1:
			Utils.log('This is a BAT-GUANO (standard) alert for event ' + follow_up_event + '.', 'info')

			rate_snr = record['rate_snr']
			rate_duration = record['rate_duration']
			rate_energy_range = record['rate_energy_range']
			classification = record['classification']
			far = record['far']
		
		# 2 - guano.loc_map.example.json
		if record_number == 2:
			Utils.log('This is a BAT-GUANO (loc_map) alert for event ' + follow_up_event + '.', 'info')

			healpix_file = record['healpix_file']
			systematic_included = record['systematic_included']
			rate_snr = record['rate_snr']
			rate_duration = record['rate_duration']
			rate_energy_range = record['rate_energy_range']
			classification = record['classification']
			far = record['far']

		# 3 - guano.loc_arc_min.example.json
		elif record_number == 3:
			Utils.log('This is a BAT-GUANO (loc_arc_min) alert for event ' + follow_up_event + '.', 'info')

			ra = record['ra']
			dec = record['dec']
			ra_dec_error = record['ra_dec_error']
			containment_probability = record['containment_probability']
			systematic_included = record['systematic_included']
			image_snr = record['image_snr']
			image_duration = record['image_duration']
			image_energy_range = record['image_energy_range']
			rate_snr = record['rate_snr']
			rate_duration = record['rate_duration']
			rate_energy_range = record['rate_energy_range']
			classification = record['classification']
			far = record['far']

		# 4 - guano.retraction
		elif record_number == 4:
			Utils.log('BAT-GUANO alert for event ' + follow_up_event + ' was retracted!', 'info')

		else:
			print('Could not read BAT-GUANO notice.')
			print(topic)
			print(value)

	@staticmethod
	def igwn_gwalert(value):

		record = json.loads(value)

		superevent_id = record['superevent_id']
		time_created = record['time_created']
		url = record['urls']['gracedb']

		Utils.log('Event ' + str(superevent_id) + ' generated at time ' + str(time_created) + '.', 'info')
		Utils.log('GraceDB URL: ' + str(url) + '.', 'info')

		rm_filter = superevent_id[0]

		if (rm_filter == 'M') or (rm_filter == 'T'):
			Utils.log('This is a MOCK event!', 'info')

		elif rm_filter == 'S':
			Utils.log('This is a REAL event!', 'info')

		else:
			Utils.log('This is an UNKNOWN event (rb_filter == ' + str(rb_filter) + ').', 'info')

		field_path = Configuration.MAIN_DIR + 'analysis/gw/fields/' + superevent_id + '-fields.txt'
		galaxy_path = Configuration.MAIN_DIR + 'analysis/gw/galaxies/' + superevent_id + '-galaxies.txt'
		json_path = Configuration.MAIN_DIR + 'alerts/gw/json/' + superevent_id + '.json'
		moc_path = Configuration.MAIN_DIR + 'alerts/gw/moc/' + superevent_id + '-90percent.moc.fit'
		skymap_path = Configuration.MAIN_DIR + 'alerts/gw/skymap/' + superevent_id + '-skymap.fit'

		# save json
		if not os.path.isfile(json_path):
			with open(json_path, 'w') as outfile:
				json.dump(record, outfile)
				Utils.log('Saving JSON for ' + superevent_id + '.', 'info')
		else:
			Utils.log('JSON already saved for ' + superevent_id + '.', 'info')

		alert_type = record['alert_type']

		# filter alert type
		if alert_type == 'RETRACTION':
			Utils.log('Alert for ' + superevent_id + ' was retracted!', 'info')

			if os.path.isfile(field_path):
				Utils.log('Deleting fields for retracted alert ' + superevent_id + '.', 'info')
				os.system('rm ' + field_path)
			else:
				Utils.log('No fields for ' + superevent_id + ' to delete.', 'info')

			if os.path.isfile(galaxy_path):
				Utils.log('Deleting galaxies for retracted alert ' + superevent_id + '.', 'info')
				os.system('rm ' + galaxy_path)
			else:
				Utils.log('No galaxies for ' + superevent_id + ' to delete.', 'info')

			if os.path.isfile(json_path):
				Utils.log('Deleting JSON for retracted alert ' + superevent_id + '.', 'info')
				os.system('rm ' + json_path)
			else:
				Utils.log('No JSON for ' + superevent_id + ' to delete.', 'info')

			if os.path.isfile(skymap_path):
				Utils.log('Deleting skymap for ' + superevent_id + '.', 'info')
				os.system('rm ' + skymap_path)
			else:
				Utils.log('No skymap for ' + superevent_id + ' to delete.', 'info')

			if os.path.isfile(moc_path):
				Utils.log('Deleting 90 percent credible region for ' + superevent_id + '.', 'info')
				os.system('rm ' + moc_path)
			else:
				Utils.log('No 90 percent credible region for ' + superevent_id + ' to delete.', 'info')

			#Utils.log('Alerting team of retraction.', 'info')
			#Alerts.alert_team(alert_type=alert_type, event_name=superevent_id)

		elif (alert_type == 'PRELIMINARY') or (alert_type == 'INITIAL') or (alert_type == 'UPDATE'):
			Utils.log('This is a ' + alert_type + ' alert for ' + superevent_id + '!', 'info')

			#Utils.log('Alerting team of event.', 'info')
			#Alerts.alert_team(alert_type=alert_type, event_name=superevent_id)

			Utils.log('Reading event parameters and skymap.', 'info')

			time = record['event']['time'] # string
			far = record['event']['far'] # float
			significant = record['event']['significant'] # bool
			instruments = record['event']['instruments'] # list

			pipeline = record['event']['pipeline'] # string
			search = record['event']['search'] # string

			print('Time:', time)
			print('FAR:', far)
			print('Significant:', significant)
			print('Instruments:', instruments)
			print('Pipeline:', pipeline)
			print('Search:', search)
			print('--------------------')
			
			group = record['event']['group'] # string

			if group == 'Burst':

				Utils.log('This is a burst event!', 'info')

				duration = record['event']['duration'] # 
				central_frequency = record['event']['central_frequency'] # 

				print('Duration:', duration, type(duration))
				print('Central frequency:', central_frequency, type(central_frequency))
				print('--------------------')

			elif group == 'CBC':

				Utils.log('This is a merger event!', 'info')

				has_ns = record['event']['properties']['HasNS'] # float
				has_remnant = record['event']['properties']['HasRemnant'] # float
				has_mass_gap = record['event']['properties']['HasMassGap'] # float

				bns_class = record['event']['classification']['BNS'] # float
				nsbh_class = record['event']['classification']['NSBH'] # float
				bbh_class = record['event']['classification']['BBH'] # float
				terra_class = record['event']['classification']['Terrestrial'] # float

				print('HasNS:', has_ns)
				print('HasRemnant:', has_remnant)
				print('HasMassGap:', has_mass_gap)
				print('--------------------')
				print('BNS Classification:', bns_class)
				print('NSBH Classification:', nsbh_class)
				print('BBH Classification:', bbh_class)
				print('Terra Classification:', terra_class)
				print('--------------------')

			skymap_str = record.get('event', {}).pop('skymap')

			if skymap_str:

				skymap_bytes = b64decode(skymap_str)

				buffer = BytesIO(skymap_bytes)

				skymap = QTable.read(buffer)

				hdr_object = skymap.meta['OBJECT'] # unique identifier for this event
				hdr_date = skymap.meta['DATE-OBS'] # UTC of observation
				hdr_mjd = skymap.meta['MJD-OBS'] # MJD of observation

				try:
					hdr_distmean = skymap.meta['DISTMEAN'] # posterior mean distance [Mpc]
				except KeyError:
					print('Could not read DISTMEAN')

				try:
					hdr_diststd = skymap.meta['DISTSTD'] # posterior std distance [Mpc]
				except KeyError:
					print('Could not read DISTSTD')
					
				# most probable sky location
				i = np.argmax(skymap['PROBDENSITY'])
				uniq = skymap[i]['UNIQ']

				level, ipix = ah.uniq_to_level_ipix(uniq)
				nside = ah.level_to_nside(level)

				ra, dec = ah.healpix_to_lonlat(ipix, nside, order='nested')
				ra_max = ra.deg
				dec_max = dec.deg

				# sort the pixels of the skymap by descending probability density
				skymap.sort('PROBDENSITY', reverse=True)

				# find the area of each pixel
				level, ipix	= ah.uniq_to_level_ipix(skymap['UNIQ'])
				pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))

				# calculate the probability within each pixel: pixel area times probability density
				prob = pixel_area * skymap['PROBDENSITY']

				# calculate the cumulative sum of the probability
				cumprob = np.cumsum(prob)

				# find the pixel for which the probability sums to 0.9
				i = cumprob.searchsorted(0.9)

				# calculate the area of the 90 percent credible region: sum of the areas of the pixels up to that one
				area_90 = pixel_area[:i].sum()
				area_90.to_value(u.deg**2)

				# save skymap to FITS
				if not os.path.isfile(skymap_path):
					Utils.log('Saving skymap for event ' + superevent_id + '.', 'info')
					skymap.write(skymap_path, overwrite=True)

				else:
					Utils.log('Skymap already saved for ' + superevent_id + '.', 'info')

				# save 90 percent credible region to FITS
				if not os.path.isfile(moc_path):
					Utils.log('Saving 90 percent credible region for event ' + superevent_id + '.', 'info')

					# keep only the pixels that are within the 90 percent credible region
					skymap90 = skymap[:i]

					# sort the pixels by their UNIQ pixel index
					skymap90.sort('UNIQ')

					# delete all columns except for the UNIQ column
					skymap90 = skymap90['UNIQ',]
					skymap90.write(moc_path, overwrite=True)

				else:
					Utils.log('MOC 90 percent credible region already saved for ' + superevent_id + '.', 'info')

			# check to see if the field list from the skymap exists
			if not os.path.isfile(field_path):
				Utils.log('Generating field list from skymap for ' + superevent_id + '.', 'info')

				# pull in the survey fields
				survey_fields = Priority.field_generator(Configuration.FIELD_SIZE)

				# generate a ranked list of fields within the skymap
				fields_prio = Priority.sort_field_skymap(survey_fields, skymap, superevent_id)
				fields_prio.to_csv(field_path, sep=' ', index=False)

			else:
				Utils.log('Field list already generated from skymap for ' + superevent_id + '.', 'info')

			# check to see if the galaxy list from the skymap exists
			if not os.path.isfile(galaxy_path):
				Utils.log('Generating galaxy list from skymap for ' + superevent_id + '.', 'info')

				# pull in the survey fields
				survey_fields = Priority.field_generator(Configuration.FIELD_SIZE)

				# pull in the galaxy catalog
				catalog = Priority.generate_targets(detection_time=None)

				# generated a ranked list of galaxies within the skymap
				galaxies_prio = Priority.sort_galaxy_skymap(catalog, skymap, superevent_id)
				galaxies_prio.to_csv(galaxy_path, sep=' ', index=False)

			else:
				Utils.log('Galaxy list already generated from skymap for ' + superevent_id + '.', 'info')

			print()		

	@staticmethod
	def process_kafka_alert(value, topic):

		if topic == 'igwn.gwalert':
			
			Utils.log('Reading alert parameters for ' + topic + ' message...', 'info')

			record = json.loads(value)
			
			superevent_id = record['superevent_id']
			time_created = record['time_created']
			url = record['urls']['gracedb']

			Utils.log('Event ' + superevent_id + ' generated at time ' + str(time_created) + '.', 'info')
			Utils.log('GraceDB URL: ' + str(url) + '.', 'info')

			rb_filter = superevent_id[0]

			if (rb_filter == 'M') or (rb_filter == 'T'):
				Utils.log('Note: this is a MOCK event!', 'info')

			elif rb_filter == 'S':
				Utils.log('Note: this is a REAL event!', 'info')

			else:
				Utils.log('Note: this is an UNKNOWN event (rb_filter == ' + str(rb_filter) + ').', 'info')

			field_path = Configuration.ANALYSIS_DIRECTORY + superevent_id + '-skymap-fields.txt'
			galaxy_path = Configuration.ANALYSIS_DIRECTORY + superevent_id + '-skymap-galaxies.txt'
			json_path = Configuration.ALERTS_DIRECTORY + superevent_id + '.json'
			skymap_path = Configuration.ALERTS_DIRECTORY + superevent_id + '-skymap.fits'
			moc_path = Configuration.ALERTS_DIRECTORY + superevent_id + '-90percent.moc.fits'

			if not os.path.isfile(json_path):
				with open(json_path, 'w') as outfile:
					json.dump(record, outfile)
					Utils.log('Saving JSON for ' + superevent_id + '.', 'info')
			else:
				Utils.log('JSON already saved for ' + superevent_id + '.', 'info')
			
			alert_type = record['alert_type']

			# filter retraction notices
			if alert_type == 'RETRACTION':
				Utils.log('Alert for ' + superevent_id + ' was retracted!', 'info')

				if os.path.isfile(field_path):
					Utils.log('Deleting fields for retracted alert ' + superevent_id + '.', 'info')
					os.system('rm ' + field_path)
				else:
					Utils.log('No fields for ' + superevent_id + ' to delete.', 'info')

				if os.path.isfile(galaxy_path):
					Utils.log('Deleting galaxies for retracted alert ' + superevent_id + '.', 'info')
					os.system('rm ' + galaxy_path)
				else:
					Utils.log('No galaxies for ' + superevent_id + ' to delete.', 'info')

				if os.path.isfile(json_path):
					Utils.log('Deleting JSON for retracted alert ' + superevent_id + '.', 'info')
					os.system('rm ' + json_path)
				else:
					Utils.log('No JSON for ' + superevent_id + ' to delete.', 'info')

				if os.path.isfile(skymap_path):
					Utils.log('Deleting skymap for ' + superevent_id + '.', 'info')
					os.system('rm ' + skymap_path)
				else:
					Utils.log('No skymap for ' + superevent_id + ' to delete.', 'info')

				if os.path.isfile(moc_path):
					Utils.log('Deleting 90 percent credible region for ' + superevent_id + '.', 'info')
					os.system('rm ' + moc_path)
				else:
					Utils.log('No 90 percent credible region for ' + superevent_id + ' to delete.', 'info')

				#Utils.log('Alerting team of retraction.', 'info')
				#Alerts.alert_team(alert_type=alert_type, event_name=superevent_id)

			elif (alert_type == 'PRELIMINARY') or (alert_type == 'INITIAL') or (alert_type == 'UPDATE'):
				Utils.log('This is a ' + alert_type + ' alert for ' + superevent_id + '!', 'info')

				#Utils.log('Alerting team of event.', 'info')
				#Alerts.alert_team(alert_type=alert_type, event_name=superevent_id)

				Utils.log('Reading event parameters and skymap.', 'info')

				time = record['event']['time']										# string
				far = record['event']['far']										# float
				significant = record['event']['significant']						# bool
				instruments = record['event']['instruments']						# list

				pipeline = record['event']['pipeline']								# string
				search = record['event']['search']									# string

				print('Time:', time)
				print('FAR:', far)
				print('Significant:', significant)
				print('Instruments:', instruments)
				print('Pipeline:', pipeline)
				print('Search:', search)
				print('--------------------')
				
				group = record['event']['group']									# string

				if group == 'Burst':

					Utils.log('This is a burst event!', 'info')

					duration = record['event']['duration']							# 
					central_frequency = record['event']['central_frequency']		# 

					print('Duration:', duration, type(duration))
					print('Central frequency:', central_frequency, type(central_frequency))
					print('--------------------')

				elif group == 'CBC':

					Utils.log('This is a merger event!', 'info')

					has_ns = record['event']['properties']['HasNS']						# float
					has_remnant = record['event']['properties']['HasRemnant']			# float
					has_mass_gap = record['event']['properties']['HasMassGap']			# float

					bns_class = record['event']['classification']['BNS']				# float
					nsbh_class = record['event']['classification']['NSBH']				# float
					bbh_class = record['event']['classification']['BBH']				# float
					terra_class = record['event']['classification']['Terrestrial']		# float

					print('HasNS:', has_ns)
					print('HasRemnant:', has_remnant)
					print('HasMassGap:', has_mass_gap)
					print('--------------------')
					print('BNS Classification:', bns_class)
					print('NSBH Classification:', nsbh_class)
					print('BBH Classification:', bbh_class)
					print('Terra Classification:', terra_class)
					print('--------------------')

				skymap_str = record.get('event', {}).pop('skymap')

				if skymap_str:

					skymap_bytes = b64decode(skymap_str)

					buffer = BytesIO(skymap_bytes)

					skymap = QTable.read(buffer)

					hdr_object = skymap.meta['OBJECT']			# unique identifier for this event
					hdr_date = skymap.meta['DATE-OBS']			# UTC of observation
					hdr_mjd = skymap.meta['MJD-OBS']			# MJD of observation

					try:
						hdr_distmean = skymap.meta['DISTMEAN']		# posterior mean distance [Mpc]
					except KeyError:
						print('Could not read DISTMEAN')

					try:
						hdr_diststd = skymap.meta['DISTSTD']		# posterior std distance [Mpc]
					except KeyError:
						print('Could not read DISTSTD')
						
					# most probable sky location
					i = np.argmax(skymap['PROBDENSITY'])
					uniq = skymap[i]['UNIQ']

					level, ipix = ah.uniq_to_level_ipix(uniq)
					nside = ah.level_to_nside(level)

					ra, dec = ah.healpix_to_lonlat(ipix, nside, order='nested')
					ra_max = ra.deg
					dec_max = dec.deg

					# sort the pixels of the skymap by descending probability density
					skymap.sort('PROBDENSITY', reverse=True)

					# find the area of each pixel
					level, ipix	= ah.uniq_to_level_ipix(skymap['UNIQ'])
					pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))

					# calculate the probability within each pixel: pixel area times probability density
					prob = pixel_area * skymap['PROBDENSITY']

					# calculate the cumulative sum of the probability
					cumprob = np.cumsum(prob)

					# find the pixel for which the probability sums to 0.9
					i = cumprob.searchsorted(0.9)

					# calculate the area of the 90 percent credible region: sum of the areas of the pixels up to that one
					area_90 = pixel_area[:i].sum()
					area_90.to_value(u.deg**2)

					# save skymap to FITS
					if not os.path.isfile(skymap_path):
						Utils.log('Saving skymap for event ' + superevent_id + '.', 'info')
						skymap.write(skymap_path, overwrite=True)

					else:
						Utils.log('Skymap already saved for ' + superevent_id + '.', 'info')

					# save 90 percent credible region to FITS
					if not os.path.isfile(moc_path):
						Utils.log('Saving 90 percent credible region for event ' + superevent_id + '.', 'info')

						# keep only the pixels that are within the 90 percent credible region
						skymap90 = skymap[:i]

						# sort the pixels by their UNIQ pixel index
						skymap90.sort('UNIQ')

						# delete all columns except for the UNIQ column
						skymap90 = skymap90['UNIQ',]
						skymap90.write(moc_path, overwrite=True)

					else:
						Utils.log('MOC 90 percent credible region already saved for ' + superevent_id + '.', 'info')

				# check to see if the field list from the skymap exists
				if not os.path.isfile(field_path):
					Utils.log('Generating field list from skymap for ' + superevent_id + '.', 'info')

					# pull in the survey fields
					survey_fields = Priority.field_generator(Configuration.FIELD_SIZE)

					# generate a ranked list of fields within the skymap
					fields_prio = Priority.sort_field_skymap(survey_fields, skymap, superevent_id)
					fields_prio.to_csv(field_path, sep=' ', index=False)

				else:
					Utils.log('Field list already generated from skymap for ' + superevent_id + '.', 'info')

				# check to see if the galaxy list from the skymap exists
				if not os.path.isfile(galaxy_path):
					Utils.log('Generating galaxy list from skymap for ' + superevent_id + '.', 'info')

					# pull in the survey fields
					survey_fields = Priority.field_generator(Configuration.FIELD_SIZE)

					# pull in the galaxy catalog
					catalog = Priority.generate_targets(detection_time=None)

					# generated a ranked list of galaxies within the skymap
					galaxies_prio = Priority.sort_galaxy_skymap(catalog, skymap, superevent_id)
					galaxies_prio.to_csv(galaxy_path, sep=' ', index=False)

				else:
					Utils.log('Galaxy list already generated from skymap for ' + superevent_id + '.', 'info')

				print()

		elif topic == 'gcn.notices.icecube.lvk_nu_track_search':
			
			record = json.loads(value)
			
			alert_type = record['type']
			reference = record['reference']
			ref_id = record['ref_ID']
			alert_datetime = record['alert_datetime']
			trigger_time = record['trigger_time']
			observation_start = record['observation_start']
			observation_stop = record['observation_stop']
			observation_livetime = record['observation_livetime']
			pval_generic = record['pval_generic']
			pval_bayesian = record['pval_bayesian']
			n_events_coincident = record['n_events_coincident']
			coincident_events = record['coincident_events']
			most_probable_direction = record['most_probable_direction']
			flux_sensitivity = record['neutrino_flux_sensitivity_range']['flux_sensitivity']
			sensitive_energy_range = record['neutrino_flux_sensitivity_range']['sensitive_energy_range']

			print('Alert Type:', alert_type, type(alert_type))
			print('Reference:', reference, type(reference))
			print('Ref ID:', ref_id, type(ref_id))
			print('Alert Datetime:', alert_datetime, type(alert_datetime))
			print('Trigger Time:', trigger_time, type(trigger_time))
			print('Observation Start:', observation_start, type(observation_start))
			print('Observation Stop:', observation_stop, type(observation_stop))
			print('Observation Livetime:', observation_livetime, type(observation_livetime))
			print('pval Generic:', pval_generic, type(pval_generic))
			print('pval Bayesian:', pval_bayesian, type(pval_bayesian))
			print('n Events Coincident:', n_events_coincident, type(n_events_coincident))
			print('Coincident Events:', coincident_events, type(coincident_events), len(coincident_events))
			print('Most Probable Direction:', most_probable_direction, type(most_probable_direction))
			print('Flux Sensitivity:', flux_sensitivity, type(flux_sensitivity))
			print('Sensitive Energy Range:', sensitive_energy_range, type(sensitive_energy_range))
			print()

		elif topic == 'gcn.notices.swift.bat.guano':

			Utils.log('Reading alert parameters for ' + topic + ' message.', 'info')

			record = json.loads(value)

			mission = record['mission']
			instrument = record['instrument']
			messenger = record['messenger']
			record_number = record['record_number']
			alert_datetime = record['alert_datetime']
			alert_tense = record['alert_tense']
			alert_type = record['alert_type']
			trigger_time = record['trigger_time']
			follow_up_event = record['follow_up_event']
			follow_up_type = record['follow_up_type']
			data_archive_page = record['data_archive_page']
			alert_id = record['id']

			Utils.log('Event BAT-GUANO/' + follow_up_event + ' generated at time ' + str(alert_datetime) + '.', 'info')
			Utils.log('Data Archive Page: ' + str(data_archive_page) + '.', 'info')

			# 1 - guano.example.json
			if record_number == 1:
				Utils.log('This is a BAT-GUANO (standard) alert for event ' + follow_up_event + '.', 'info')

				rate_snr = record['rate_snr']
				rate_duration = record['rate_duration']
				rate_energy_range = record['rate_energy_range']
				classification = record['classification']
				far = record['far']
			
			# 2 - guano.loc_map.example.json
			if record_number == 2:
				Utils.log('This is a BAT-GUANO (loc_map) alert for event ' + follow_up_event + '.', 'info')

				healpix_file = record['healpix_file']
				systematic_included = record['systematic_included']
				rate_snr = record['rate_snr']
				rate_duration = record['rate_duration']
				rate_energy_range = record['rate_energy_range']
				classification = record['classification']
				far = record['far']

			# 3 - guano.loc_arc_min.example.json
			elif record_number == 3:
				Utils.log('This is a BAT-GUANO (loc_arc_min) alert for event ' + follow_up_event + '.', 'info')

				ra = record['ra']
				dec = record['dec']
				ra_dec_error = record['ra_dec_error']
				containment_probability = record['containment_probability']
				systematic_included = record['systematic_included']
				image_snr = record['image_snr']
				image_duration = record['image_duration']
				image_energy_range = record['image_energy_range']
				rate_snr = record['rate_snr']
				rate_duration = record['rate_duration']
				rate_energy_range = record['rate_energy_range']
				classification = record['classification']
				far = record['far']

			# 4 - guano.retraction
			elif record_number == 4:
				Utils.log('BAT-GUANO alert for event ' + follow_up_event + ' was retracted!', 'info')

			else:
				print('Could not read BAT-GUANO notice.')
				print(topic)
				print(value)