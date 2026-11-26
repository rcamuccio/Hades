'''
Ἄγγελος - The Messenger
'''

import json
import requests
import smtplib
from email.message import EmailMessage

from conf import Configuration

class Angelos:

	@staticmethod
	def decode_classic_notice(value):
		'''This function decodes the value of a classic text GCN notice.

		:parameter value [gcn_kafka.Consumer.value] - the notice value

		:return params [dict] - the dictionary of message value parameters
		'''

		# split the value into a list
		contents = value.decode('utf-8').split('\n')

		keys = []
		items = []

		for content in contents:
			if len(content.split(': ')) > 1:
				# extract the key
				key = content.split(': ')[0].strip()
				keys.append(key)

				# extract the item
				item = content.split(': ')[1].strip()
				items.append(item)

		# package the key/item pairs into a dictionary
		params = dict(zip(keys, items))

		return params

	@staticmethod
	def filter_notice(value, topic):
		'''This function redirects an incoming notice to the relevant processing function.

		:parameter value [gcn_kakfa.Consumer.value] - the notice value
		:parameter topic [gcn_kafka.Consumer.topic] - the notice topic

		:return - nothing
		'''

		line_header = '\033[1m[Ἄγγελος | filter_notice]\033[0m : '
		print(line_header + 'Alert received (' + str(topic) + ')')

		if topic == 'gcn.circulars':
			Angelos.process_gcn_circulars(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_GBM_ALERT':
			Angelos.process_fermi_gbm_alert(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_GBM_FIN_POS':
			Angelos.process_fermi_gbm_fin_pos(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_GBM_FLT_POS':
			Angelos.process_fermi_gbm_flt_pos(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_GBM_GND_POS':
			Angelos.process_fermi_gbm_gnd_pos(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_GBM_POS_TEST':
			Angelos.process_fermi_gbm_pos_test(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_GBM_SUBTHRESH':
			Angelos.process_fermi_gbm_subthresh(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_LAT_MONITOR':
			Angelos.process_fermi_lat_monitor(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_LAT_POS_DIAG':
			Angelos.process_fermi_lat_pos_diag(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_LAT_POS_INI':
			Angelos.process_fermi_lat_pos_ini(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_LAT_OFFLINE':
			Angelos.process_fermi_lat_pos_offline(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_LAT_POS_TEST':
			Angelos.process_fermi_lat_pos_test(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_LAT_POS_UPD':
			Angelos.process_fermi_lat_pos_upd(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_POINTDIR':
			Angelos.process_fermi_pointdir(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_ACTUAL_POINTDIR':
			Angelos.process_swift_actual_pointdir(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_BAT_GRB_LC':
			Angelos.process_swift_bat_grb_lc(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_BAT_GRB_POS_ACK':
			Angelos.process_swift_bat_grb_pos_ack(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_BAT_GRB_POS_TEST':
			Angelos.process_swift_bat_grb_pos_test(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_BAT_QL_POS':
			Angelos.process_swift_bat_grb_ql_pos(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_BAT_SCALEDMAP':
			Angelos.process_swift_bat_scaledmap(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_BAT_TRANS':
			Angelos.process_swift_bat_trans(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_FOM_OBS':
			Angelos.process_swift_fom_obs(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_POINTDIR':
			Angelos.process_swift_pointdir(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_SC_SLEW':
			Angelos.process_swift_sc_slew(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_TOO_FOM':
			Angelos.process_swift_too_fom(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_TOO_SC_SLEW':
			Angelos.process_swift_too_sc_slew(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_DBURST':
			Angelos.process_swift_uvot_dburst(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_DBURST_PROC':
			Angelos.process_swift_uvot_dburst_proc(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_EMERGENCY':
			Angelos.process_swift_uvot_emergency(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_FCHART':
			Angelos.process_swift_uvot_fchart(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_FCHART_PROC':
			Angelos.process_swift_uvot_fchart_proc(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_POS':
			Angelos.process_swift_uvot_pos(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_POS_NACK':
			Angelos.process_swift_uvot_pos_nack(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_CENTROID':
			Angelos.process_swift_xrt_centroid(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_EMERGENCY':
			Angelos.process_swift_xrt_emergency(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_IMAGE':
			Angelos.process_swift_xrt_image(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_IMAGE_PROC':
			Angelos.process_swift_xrt_image_proc(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_LC':
			Angelos.process_swift_xrt_lc(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_POSITION':
			Angelos.process_swift_xrt_position(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_SPECTRUM':
			Angelos.process_swift_xrt_spectrum(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_SPECTRUM_PROC':
			Angelos.process_swift_xrt_spectrum_proc(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_SPER':
			Angelos.process_swift_xrt_sper(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_SPER_PROC':
			Angelos.process_swift_xrt_sper_proc(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_THRESHPIX':
			Angelos.process_swift_xrt_threshpix(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_THRESHPIX_PROC':
			Angelos.process_swift_xrt_threshpix_proc(value, topic, False)

		elif topic == 'gcn.notices.chime.frb':
			Angelos.process_chime_frb(value, topic, False)

		elif topic == 'gcn.notices.dsa110.frb':
			Angelos.process_dsa110_frb(value, topic, False)

		elif topic == 'gcn.notices.einstein_probe.wxt.alert':
			Angelos.process_einstein_probe_wxt_alert(value, topic, False)

		elif topic == 'gcn.notices.icecube.gold_bronze_track_alerts':
			Angelos.process_icecube_gold_bronze_track_alerts(value, topic, False)

		elif topic == 'gcn.notices.icecube.lvk_nu_track_search':
			Angelos.process_icecube_lvk_nu_track_search(value, topic, False)

		elif topic == 'gcn.notices.maxi.gsc.known':
			Angelos.process_maxi_gsc_known(value, topic, False)

		elif topic == 'gcn.notices.maxi.gsc.unknown':
			Angelos.process_maxi_gsc_unknown(value, topic, False)

		elif topic == 'gcn.notices.superk.sn_alert':
			Angelos.process_superk_sn_alert(value, topic, False)

		elif topic == 'gcn.notices.svom.voevent.eclairs':
			Angelos.process_svom_eclairs(value, topic, False)

		elif topic == 'gcn.notices.swift.bat.guano':
			Angelos.process_swift_bat_guano(value, topic, False)

		elif topic == 'igwn.gwalert':
			Angelos.process_igwn_gwalert(value, topic, False)

	@staticmethod
	def get_notice_topics():
		'''This function reads a file listing available GCN Notice topics and returns them as a list.

		:parameter - nothing

		:return topic_list [list][str] - the list of available notice topics
		:return num_topics [int] - the number of available notice topics
		'''

		topic_path = Configuration.PLOUTON_DIRECTORY + 'alert/gcn_topics' + Configuration.TABLE_EXTENSION

		topic_file = open(topic_path, 'r')
		topic_list = []

		for ln in topic_file:
			if ln[0] == '#':
				pass
			else:
				ln = ln.rstrip('\n')
				topic_list.append(ln)

		num_topics = len(topic_list)

		return topic_list, num_topics

	@staticmethod
	def process_chime_frb(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_dsa110_frb(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_einstein_probe_wxt_alert(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_fermi_gbm_alert(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_fermi_gbm_fin_pos(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_fermi_gbm_flt_pos(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_fermi_gbm_gnd_pos(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_fermi_gbm_pos_test(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_fermi_gbm_subthresh(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_fermi_lat_monitor(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_fermi_lat_pos_diag(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_fermi_lat_pos_ini(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_fermi_lat_pos_offline(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_fermi_lat_pos_test(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_fermi_lat_pos_upd(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_fermi_pointdir(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_gcn_circulars(value, topic, alert=Configuration.SEND_ALERT):
		'''This function processes GCN Circulars via GCN Kafka.

		:parameter value [gcn_kakfa.Consumer.value] - the notice value
		:parameter topic [gcn_kafka.Consumer.topic] - the notice topic
		:parameter alert [bool] - a toggle for broadcasting the alert

		:return - nothing
		'''

		record = json.loads(value)
		body = ''

		rec_schema = record['$schema']
		str_schema = '\tSchema: ' + str(rec_schema) + '\n'
		body += str_schema

		rec_event_id = record['eventId']
		str_event_id = '\tEvent ID: ' + str(rec_event_id) + '\n'
		body += str_event_id

		rec_submitter = record['submitter']
		str_submitter = '\tSubmitter: ' + str(rec_submitter) + '\n'
		body += str_submitter

		rec_submitted_how = record['submittedHow']
		str_submitted_how = '\tSubmitted How: ' + str(rec_submitted_how) + '\n'
		body += str_submitted_how

		rec_subject = record['subject']
		str_subject = '\tSubject: ' + str(rec_subject) + '\n'
		body += str_subject

		rec_circular_id = record['circularId']
		str_circular_id = '\tCircular ID: ' + str(rec_circular_id) + '\n'
		body += str_circular_id

		rec_format = record['format']
		str_format = '\tFormat: ' + str(rec_format) + '\n'
		body += str_format

		rec_body = record['body']
		str_body = '\tBody: ' + str(rec_body) + '\n'
		body += str_body

		rec_created_on = record['createdOn']
		str_created_on = '\tCreated On: ' + str(rec_created_on) + '\n'
		body += str_created_on

		print(body)

		if alert:
			Angelos.send_alert(topic, body)

	@staticmethod
	def process_icecube_gold_bronze_track_alerts(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_icecube_lvk_nu_track_search(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_igwn_gwalert(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_maxi_gsc_known(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_maxi_gsc_unknown(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_superk_sn_alert(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_svom_eclairs(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_actual_pointdir(value, topic, alert=Configuration.SEND_ALERT):
		'''This function processes Swift Actual Pointing Direction Notices from the Neil Gehrels Swift Observatory via GCN Classic over Kafka.

		:parameter value [gcn.kafka.Consumer.value] - the notice value
		:parameter topic [gcn.kafka.Consumer.topic] - the notice topic
		:parameter alert [bool] - a toggle for broadcasting the alert

		:return - nothing
		'''

		record = Angelos.decode_classic_notice(value)
		body = ''

		rec_title = record['TITLE']
		str_title = '\tTitle: ' + str(rec_title) + '\n'
		body += str_title

		rec_notice_date = record['NOTICE_DATE']
		str_notice_date = '\tNotice Date: ' + str(rec_notice_date) + '\n'
		body += str_notice_date

		rec_notice_type = record['NOTICE_TYPE']
		str_notice_type = '\tNotice Type: ' + str(rec_notice_type) + '\n'
		body += str_notice_type

		rec_curr_point_ra = record['CURR_POINT_RA']
		str_curr_point_ra = '\tCurrent Point (RA): ' + str(rec_curr_point_ra) + '\n'
		body += str_curr_point_ra

		rec_curr_point_dec = record['CURR_POINT_DEC']
		str_curr_point_dec = '\tCurrent Point (Dec): ' + str(rec_curr_point_dec) + '\n'
		body += str_curr_point_dec

		rec_curr_point_roll = record['CURR_POINT_ROLL']
		str_curr_point_roll = '\tCurrent Point (Roll): ' + str(rec_curr_point_roll) + '\n'
		body += str_curr_point_roll

		rec_slew_time = record['SLEW_TIME']
		str_slew_time = '\tSlew Time: ' + str(rec_slew_time) + '\n'
		body += str_slew_time

		rec_slew_date = record['SLEW_DATE']
		str_slew_date = '\tSlew Date: ' + str(rec_slew_date) + '\n'
		body += str_slew_date

		rec_tgt_num = record['TGT_NUM']
		str_tgt_num = '\tTarget Number: ' + str(rec_tgt_num) + '\n'
		body += str_tgt_num

		rec_sun_postn = record['SUN_POSTN']
		str_sun_postn = '\tSun Position: ' + str(rec_sun_postn) + '\n'
		body += str_sun_postn
		
		rec_sun_dist = record['SUN_DIST']
		str_sun_dist = '\tSun Distance: ' + str(rec_sun_dist) + '\n'
		body += str_sun_dist
		
		rec_moon_postn = record['MOON_POSTN']
		str_moon_postn = '\tMoon Position: ' + str(rec_moon_postn) + '\n'
		body += str_moon_postn
		
		rec_moon_dist = record['MOON_DIST']
		str_moon_dist = '\tMoon Distance: ' + str(rec_moon_dist) + '\n'
		body += str_moon_dist
		
		rec_moon_illum = record['MOON_ILLUM']
		str_moon_illum = '\tMoon Illumination: ' + str(rec_moon_illum) + '\n'
		body += str_moon_illum

		rec_gal_coords = record['GAL_COORDS']
		str_gal_coords = '\tGalactic Coordinates: ' + str(rec_gal_coords) + '\n'
		body += str_gal_coords

		rec_ecl_coords = record['ECL_COORDS']
		str_ecl_coords = '\tEcliptic Coordinates: ' + str(rec_ecl_coords) + '\n'
		body += str_ecl_coords

		rec_comments = record['COMMENTS']
		str_comments = '\tComments: ' + str(rec_comments) + '\n'
		body += str_comments

		print(body)

		if alert:
			Angelos.send_alert(topic, body)

	@staticmethod
	def process_swift_bat_grb_lc(value, topic, alert=Configuration.SEND_ALERT):
		'''This function processes Swift BAT Lightcurve Notices from the Neil Gehrels Swift Observatory via GCN Classic over Kafka.

		:parameter value [gcn.kafka.Consumer.value] - the notice value
		:parameter topic [gcn.kafka.Consumer.topic] - the notice topic
		:parameter alert [bool] - a toggle for broadcasting the alert

		:return - nothing
		'''

		record = Angelos.decode_classic_notice(value)
		body = ''

		rec_title = record['TITLE']
		str_title = '\tTitle: ' + str(rec_title) + '\n'
		body += str_title

		rec_notice_date = record['NOTICE_DATE']
		str_notice_date = '\tNotice Date: ' + str(rec_notice_date) + '\n'
		body += str_notice_date

		rec_notice_type = record['NOTICE_TYPE']
		str_notice_type = '\tNotice Type: ' + str(rec_notice_type) + '\n'
		body += str_notice_type

		rec_trigger_num = record['TRIGGER_NUM']
		str_trigger_num = '\tTrigger Number: ' + str(rec_trigger_num) + '\n'
		body += str_trigger_num

		rec_grb_ra = record['GRB_RA']
		str_grb_ra = '\tGRB RA: ' + str(rec_grb_ra) + '\n'
		body += str_grb_ra

		rec_grb_dec = record['GRB_DEC']
		str_grb_dec = '\tGRB Dec: ' + str(rec_grb_dec) + '\n'
		body += str_grb_dec

		rec_grb_date = record['GRB_DATE']
		str_grb_date = '\tGRB Date: ' + str(rec_grb_date) + '\n'
		body += str_grb_date

		rec_grb_time = record['GRB_TIME']
		str_grb_time = '\tGRB Time: ' + str(rec_grb_time) + '\n'
		body += str_grb_time

		rec_trigger_index = record['TRIGGER_INDEX']
		str_trigger_index = '\tTrigger Index: ' + str(rec_trigger_index) + '\n'
		body += str_trigger_index

		rec_grb_phi = record['GRB_PHI']
		str_grb_phi = '\tGRB Phi: ' + str(rec_grb_phi) + '\n'
		body += str_grb_phi

		rec_grb_theta = record['GRB_THETA']
		str_grb_theta = '\tGRB Theta: ' + str(rec_grb_theta) + '\n'
		body += str_grb_theta

		rec_delta_time = record['DELTA_TIME']
		str_delta_time = '\tDelta Time: ' + str(rec_delta_time) + '\n'
		body += str_delta

		rec_lc_url = record['LC_URL']
		str_lc_url = '\tLightcurve URL: ' + str(rec_lc_url) + '\n'
		body += str_lc_url

		rec_sun_postn = record['SUN_POSTN']
		str_sun_postn = '\tSun Position: ' + str(rec_sun_postn) + '\n'
		body += str_sun_postn
		
		rec_sun_dist = record['SUN_DIST']
		str_sun_dist = '\tSun Distance: ' + str(rec_sun_dist) + '\n'
		body += str_sun_dist
		
		rec_moon_postn = record['MOON_POSTN']
		str_moon_postn = '\tMoon Position: ' + str(rec_moon_postn) + '\n'
		body += str_moon_postn
		
		rec_moon_dist = record['MOON_DIST']
		str_moon_dist = '\tMoon Distance: ' + str(rec_moon_dist) + '\n'
		body += str_moon_dist
		
		rec_moon_illum = record['MOON_ILLUM']
		str_moon_illum = '\tMoon Illumination: ' + str(rec_moon_illum) + '\n'
		body += str_moon_illum
		
		rec_gal_coords = record['GAL_COORDS']
		str_gal_coords = '\tGalactic Coordinates: ' + str(rec_gal_coords) + '\n'
		body += str_gal_coords
		
		rec_ecl_coords = record['ECL_COORDS']
		str_ecl_coords = '\tEcliptic Coordinates: ' + str(rec_ecl_coords) + '\n'
		body += str_ecl_coords
		
		rec_comments = record['COMMENTS']
		str_comments = '\tComments: ' + str(rec_comments) + '\n'
		body += str_comments

		print(body)

		if alert:
			Angelos.send_alert(topic, body)

	@staticmethod
	def process_swift_bat_grb_pos_ack(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_bat_grb_pos_test(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_bat_grb_ql_pos(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_bat_guano(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_bat_scaledmap(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_bat_trans(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_fom_obs(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_pointdir(value, topic, alert=Configuration.SEND_ALERT):
		'''This function processes Swift Pointing Direction Notices from the Neil Gehrels Swift Observatory via GCN Classic over Kafka.

		:parameter value [gcn.kafka.Consumer.value] - the notice value
		:parameter topic [gcn.kafka.Consumer.topic] - the notice topic
		:parameter alert [bool] - a toggle for broadcasting the alert

		:return - nothing
		'''

		record = Angelos.decode_classic_notice(value)
		body = ''

		rec_title = record['TITLE']
		str_title = '\tTitle: ' + str(rec_title) + '\n'
		body += str_title

		rec_notice_date = record['NOTICE_DATE']
		str_notice_date = '\tNotice Date: ' + str(rec_notice_date) + '\n'
		body += str_notice_date

		rec_notice_type = record['NOTICE_TYPE']
		str_notice_type = '\tNotice Type: ' + str(rec_notice_type) + '\n'
		body += str_notice_type
		
		rec_next_point_ra = record['NEXT_POINT_RA']
		str_next_point_ra = '\tNext Point (RA): ' + str(rec_next_point_ra) + '\n'
		body += str_next_point_ra
		
		rec_next_point_dec = record['NEXT_POINT_DEC']
		str_next_point_dec = '\tNext Point (Dec): ' + str(rec_next_point_dec) + '\n'
		body += str_next_point_dec
		
		rec_next_point_roll = record['NEXT_POINT_ROLL']
		str_next_point_roll = '\tNext Point (Roll): ' + str(rec_next_point_roll) + '\n'
		body += str_next_point_roll
		
		rec_slew_time = record['SLEW_TIME']
		str_slew_time = '\tSlew Time: ' + str(rec_slew_time) + '\n'
		body += str_slew_time
		
		rec_slew_date = record['SLEW_DATE']
		str_slew_date = '\tSlew Date: ' + str(rec_slew_date) + '\n'
		body += str_slew_date
		
		rec_obs_time = record['OBS_TIME']
		str_obs_time = '\tObservation Time: ' + str(rec_obs_time) + '\n'
		body += str_obs_time
		
		rec_tgt_name = record['TGT_NAME']
		str_tgt_name = '\tTarget Name: ' + str(rec_tgt_name) + '\n'
		body += str_tgt_name
		
		rec_tgt_num = record['TGT_NUM']
		str_tgt_num = '\tTarget Number: ' + str(rec_tgt_num) + '\n'
		body += str_tgt_num
		
		rec_merit = record['MERIT']
		str_merit = '\tMerit: ' + str(rec_merit) + '\n'
		body += str_merit
		
		rec_inst_modes = record['INST_MODES']
		str_inst_modes = '\tInstrument Modes: ' + str(rec_inst_modes) + '\n'
		body += str_inst_modes
		
		rec_sun_postn = record['SUN_POSTN']
		str_sun_postn = '\tSun Position: ' + str(rec_sun_postn) + '\n'
		body += str_sun_postn
		
		rec_sun_dist = record['SUN_DIST']
		str_sun_dist = '\tSun Distance: ' + str(rec_sun_dist) + '\n'
		body += str_sun_dist
		
		rec_moon_postn = record['MOON_POSTN']
		str_moon_postn = '\tMoon Position: ' + str(rec_moon_postn) + '\n'
		body += str_moon_postn
		
		rec_moon_dist = record['MOON_DIST']
		str_moon_dist = '\tMoon Distance: ' + str(rec_moon_dist) + '\n'
		body += str_moon_dist
		
		rec_moon_illum = record['MOON_ILLUM']
		str_moon_illum = '\tMoon Illumination: ' + str(rec_moon_illum) + '\n'
		body += str_moon_illum
		
		rec_gal_coords = record['GAL_COORDS']
		str_gal_coords = '\tGalactic Coordinates: ' + str(rec_gal_coords) + '\n'
		body += str_gal_coords
		
		rec_ecl_coords = record['ECL_COORDS']
		str_ecl_coords = '\tEcliptic Coordinates: ' + str(rec_ecl_coords) + '\n'
		body += str_ecl_coords
		
		rec_comments = record['COMMENTS']
		str_comments = '\tComments: ' + str(rec_comments) + '\n'
		body += str_comments

		print(body)

		if alert:
			Angelos.send_alert(topic, body)

	@staticmethod
	def process_swift_sc_slew(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_too_fom(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_too_sc_slew(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_uvot_dburst(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_uvot_dburst_proc(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_uvot_emergency(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_uvot_fchart(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_uvot_fchart_proc(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_uvot_pos(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_uvot_pos_nack(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_xrt_centroid(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_xrt_emergency(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_xrt_image(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_xrt_image_proc(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_xrt_lc(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_xrt_position(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_xrt_spectrum(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_xrt_spectrum_proc(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_xrt_sper(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_xrt_sper_proc(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_xrt_threshpix(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def process_swift_xrt_threshpix_proc(value, topic, alert=Configuration.SEND_ALERT):
		return

	@staticmethod
	def send_message(topic, body):
		'''This function creates and sends an alert message.

		:parameter topic [gcn.kafka.Consumer.topic] - the message topic
		:parameter body [str] - the message body

		:return - nothing
		'''

		# log into the server
		server = smtplib.SMTP_SSL(Configuration.EMAIL_SMTP, Configuration.EMAIL_PORT)
		server.ehlo()
		server.login(Configuration.EMAIL_ID, Configuration.EMAIL_PASSWORD)

		# construct the message
		message = EmailMessage()
		message['Subject'] = 'TOROS Alert: ' + str(topic)
		message['From'] = Configuration.EMAIL_ID
		message['To'] = ', '.join(Configuration.MAILING_LIST)
		message.set_content(body)

		# send the message
		server.send_message(message)
		server.quit()