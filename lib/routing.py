from config import Configuration
from lib.utilities import Utils
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import json
import requests
import smtplib

class Router:

	@staticmethod
	def filter_alert(value, topic):

		# GCN Kafka Alerts
		if topic == 'gcn.circulars':
			Utils.log('GCN Kafka Alert (' + str(topic) + ').', 'info')
			Router.gcn_circulars(value, topic)

		elif topic == 'gcn.heartbeat':
			Router.gcn_heartbeat(value, topic)

		elif topic == 'gcn.notices.einstein_probe.wxt.alert':
			Utils.log('GCN Kafka Alert (' + str(topic) + ').', 'info')
			Router.gcn_notices_einstein_probe_wxt_alert(value, topic)

		elif topic == 'gcn.notices.icecube.lvk_nu_track_search':
			Utils.log('GCN Kafka Alert (' + str(topic) + ').', 'info')
			Router.gcn_notices_icecube_lvk_nu_track_search(value, topic)

		elif topic =='gcn.notices.superk.sn_alert':
			Utils.log('GCN Kafka Alert (' + str(topic) + ').', 'info')
			Router.gcn_notices_superk_sn_alert(value, topic)

		elif topic == 'gcn.notices.swift.bat.guano':
			Utils.log('GCN Kafka Alert (' + str(topic) + ').', 'info')
			Router.gcn_notices_swift_bat_guano(value, topic)

		elif topic == 'igwn.gwalert':
			Utils.log('GCN Kafka Alert (' + str(topic) + ').', 'info')
			Router.igwn_gwalert(value, topic)

		# Unknown Alerts
		else:
			Utils.log('Unknown Alert (' + str(topic) + '). Printing value to terminal.', 'info')
			print(value)

	@staticmethod
	def gcn_circulars(value, topic):

		record = json.loads(value)

		try:
			schema = record['$schema']
		except KeyError:
			schema = None

		try:
			event_id = record['eventId']
		except KeyError:
			event_id = None

		try:
			submitter = record['submitter']
		except KeyError:
			submitter = None

		try:
			submitted_how = record['submittedHow']
		except KeyError:
			submitted_how = None

		try:
			subject = record['subject']
		except KeyError:
			subject = None

		try:
			circular_id = record['circularId']
		except KeyError:
			circular_id = None

		try:
			fformat = record['format']
		except KeyError:
			fformat = None

		try:
			body = record['body']
		except KeyError:
			body = None

		try:
			created_on = record['createdOn']
		except KeyError:
			created_on = None

		server = smtplib.SMTP_SSL(Configuration.SMTP, Configuration.PORT)
		server.ehlo()
		server.login(Configuration.EMAIL, Configuration.PAS)
		msg = MIMEMultipart()
		msg['From'] = Configuration.EMAIL
		msg['To'] = ', '.join(Configuration.MAILING_LIST)
		msg['Subject'] = 'Alert Received: GCN Circular #' + str(circular_id) + '\n'
		body = str(subject) + '\n\n' + str(body)
		msg.attach(MIMEText(body, 'plain'))
		sms = msg.as_string()
		server.sendmail(Configuration.EMAIL, Configuration.MAILING_LIST, sms)
		server.quit()

	@staticmethod
	def gcn_heartbeat(value, topic):

		record = json.loads(value)

	@staticmethod
	def gcn_notices_einstein_probe_wxt_alert(value, topic):

		record = json.loads(value)
		
		try:
			schema = record['$schema']
		except KeyError:
			schema = None
		try:
			instrument = record['instrument']
		except KeyError:
			instrument = None
		try:
			trigger_time = record['trigger_time']
		except KeyError:
			trigger_time = None
		try:
			iid = record['id']
		except KeyError:
			iid = None
		try:
			ra = record['ra']
		except KeyError:
			ra = None
		try:
			dec = record['dec']
		except KeyError:
			dec = None
		try:
			ra_dec_error = record['ra_dec_error']
		except KeyError:
			ra_dec_error = None
		try:
			image_energy_range = record['image_energy_range']
			image_energy_range_min = image_energy_range[0]
			image_energy_range_max = image_energy_range[1]
		except KeyError:
			image_energy_range = None
			image_energy_range_min = None
			image_energy_range_max = None
		try:
			net_count_rate = record['net_count_rate']
		except KeyError:
			net_count_rate = None
		try:
			image_snr = record['image_snr']
		except KeyError:
			image_snr = None
		try:
			additional_info = record['additional_info']
		except KeyError:
			additional_info = None

	@staticmethod
	def gcn_notices_icecube_lvk_nu_track_search(value, topic):

		record = json.loads(value)

		try:
			schema = record['$schema']
		except KeyError:
			schema = None
		try:
			ttype = record['type']
		except KeyError:
			ttype = None
		try:
			reference = record['reference']
			gcn_notices_lvk_alert = reference['gcn.notices.LVK.alert']
		except KeyError:
			reference = None
			gcn_notices_lvk_alert = None
		try:
			ref_id = record['ref_ID']
		except KeyError:
			ref_id = None
		try:
			alert_datetime = record['alert_datetime']
		except KeyError:
			alert_datetime = None
		try:
			trigger_time = record['trigger_time']
		except KeyError:
			trigger_time = None
		try:
			observation_start = record['observation_start']
		except KeyError:
			observation_start = None
		try:
			observation_stop = record['observation_stop']
		except KeyError:
			observation_stop = None
		try:
			observation_livetime = record['observation_livetime']
		except KeyError:
			observation_livetime = None
		try:
			pval_generic = record['pval_generic']
		except KeyError:
			pval_generic = None
		try:
			pval_bayesian = record['pval_bayesian']
		except KeyError:
			pval_bayesian = None
		try:
			n_events_coincident = record['n_events_coincident']
			if n_events_coincident > 0:
				
		except KeyError:
			n_events_coincident = None
		if (n_events_coincident != None) & (n_events_coincident > 0):
			try:
				coincident_events = record['coincident_events']
			except KeyError:
				coincident_events = None
			for i in range(n_events_coincident):
				coincident_event = coincident_events[i]
				try:
					event_dt = coincident_event['event_dt']
				except KeyError:
					event_dt = None
				try:
					localization = coincident_event['localization']
				except KeyError:
					localization = None
				if localization != None:
					try:
						ra = localization['ra']
					except KeyError:
						ra = None
					try:
						dec = localization['dec']
					except KeyError:
						dec = None
					try:
						ra_dec_error = localization['ra_dec_error']
					except KeyError:
						ra_dec_error = None
					try:
						containment_probability = localization['containment_probability']
					except KeyError:
						containment_probability = None
					try:
						systematic_included = localization['systematic_included']
					except KeyError:
						systematic_included = None
				try:
					iid = coincident_event['id']
				except KeyError:
					iid = None
				try:
					event_pval_generic = coincident_event['event_pval_generic']
				except KeyError:
					event_pval_generic = None
				try:
					event_pval_bayesian = coincident_event['event_pval_bayesian']
				except KeyError:
					event_pval_bayesian = None
		try:
			most_probable_direction = record['most_probable_direction']
			ra = most_probable_direction['ra']
			dec = most_probable_direction['dec']
		except KeyError:
			most_probable_direction = None
			ra = None
			dec = None
		try:
			neutrino_flux_sensitivity_range = record['neutrino_flux_sensitivity_range']
			try:	
				flux_sensitivity = neutrino_flux_sensitivity_range['flux_sensitivity']
				flux_sensitivity_min = flux_sensitivity[0]
				flux_sensitivity_max = flux_sensitivity[1]
			except KeyError:
				flux_sensitivity = None
				flux_sensitivity_min = None
				flux_sensitivity_max = None
			try:	
				sensitive_energy_range = neutrino_flux_sensitivity_range['sensitive_energy_range']
				sensitive_energy_range_min = sensitive_energy_range[0]
				sensitive_energy_range_max = sensitive_energy_range[1]
			except KeyError:
				sensitive_energy_range = None
				sensitive_energy_range_min = None
				sensitive_energy_range_max = None
		except KeyError:
			neutrino_flux_sensitivity_range = None
			flux_sensitivity = None
			flux_sensitivity_min = None
			flux_sensitivity_max = None
			sensitive_energy_range = None
			sensitive_energy_range_min = None
			sensitive_energy_range_max = None

	@staticmethod
	def gcn_notices_superk_sn_alert(value, topic):

		record = json.loads(value)
		
		try:
			schema = record['$schema']
		except KeyError:
			schema = None

		try:
			mission = record['mission']
		except KeyError:
			mission = None

		try:
			messenger = record['messenger']
		except KeyError:
			messenger = None

		try:
			iid = record['id']
		except KeyError:
			idd = None

		try:
			record_number = record['record_number']
		except KeyError:
			record_number = None

		try:
			trigger_number = record['trigger_number']
		except KeyError:
			trigger_number = None

		try:
			alert_datetime = record['alert_datetime']
		except KeyError:
			alert_datetime = None

		try:
			alert_tense = record['current']
		except KeyError:
			alert_tense = None

		try:
			alert_type = record['alert_type']
		except KeyError:
			alert_type = None

		try:
			trigger_time = record['trigger_time']
		except KeyError:
			trigger_time = None

		try:
			processed_sample = record['processed_sample']
		except KeyError:
			processed_sample = None

		try:
			pipeline = record['pipeline']
		except KeyError:
			pipeline = None

		try:
			n_events = record['n_events']
		except KeyError:
			n_events = None

		try:
			n_ibd_events = record['n_ibd_events']
		except KeyError:
			n_ibd_events = None

		try:
			detection_interval = record['detection_interval']
		except KeyError:
			detection_interval = None

		try:
			rate_energy_range = record['rate_energy_range']
			rate_energy_range_min = rate_energy_range[0]
			rate_energy_range_max = rate_energy_range[1]
		except KeyError:
			rate_energy_range = None

		try:
			ra = record['ra']
		except KeyError:
			ra = None

		try:
			dec = record['dec']
		except KeyError:
			dec = None

		try:
			ra_dec_error = record['ra_dec_error']
		except KeyError:
			ra_dec_error = None

		try:
			containment_probability = record['containment_probability']
		except KeyError:
			containment_probability = None

		try:
			luminosity_distance = record['luminosity_distance']
		except KeyError:
			luminosity_distance = None

		try:
			luminosity_distance_error = record['luminosity_distance_error']
		except KeyError:
			luminosity_distance_error = None

		server = smtplib.SMTP_SSL(Configuration.SMTP, Configuration.PORT)
		server.ehlo()
		server.login(Configuration.EMAIL, Configuration.PAS)
		msg = MIMEMultipart()
		msg['From'] = Configuration.EMAIL
		msg['To'] = ', '.join(Configuration.MAILING_LIST)
		msg['Subject'] = 'ALERT: Super-Kamioka Neutrino Detection Experiment (Super-Kamiokande)\n'
		body = 'Text body.'
		msg.attach(MIMEText(body, 'plain'))
		sms = msg.as_string()
		server.sendmail(Configuration.EMAIL, Configuration.MAILING_LIST, sms)
		server.quit()
		Utils.log('Alert message sent to mailing list recipients [' + topic + '].', 'info')

	@staticmethod
	def gcn_notices_swift_bat_guano(value, topic):

		record = json.loads(value)
		
		try:
			mission = record['mission']
		except KeyError:
			mission = None

		try:
			instrument = record['instrument']
		except KeyError:
			instrument = None

		try:
			messenger = record['messenger']
		except KeyError:
			messenger = None

		try:
			record_number = record['record_number']
		except KeyError:
			record_number = None

		try:
			alert_datetime = record['alert_datetime']
		except KeyError:
			alert_datetime = None

		try:
			alert_tense = record['alert_tense']
		except KeyError:
			alert_tense = None

		try:
			alert_type = record['alert_type']
		except KeyError:
			alert_type = None

		try:
			trigger_time = record['trigger_time']
		except KeyError:
			trigger_time = None

		try:
			follow_up_event = record['follow_up_event']
		except KeyError:
			follow_up_event = None

		try:
			follow_up_type = record['follow_up_type']
		except KeyError:
			follow_up_type = None

		try:
			data_archive_page = record['data_archive_page']
		except KeyError:
			data_archive_page = None

		try:
			alert_id = record['id']
		except KeyError:
			alert_id = None

		# 1 - guano.example.json
		if record_number == 1:

			try:
				rate_snr = record['rate_snr']
			except KeyError:
				rate_snr = None

			try:
				rate_duration = record['rate_duration']
			except KeyError:
				rate_duration = None

			try:
				rate_energy_range = record['rate_energy_range']
				rate_energy_range_min = rate_energy_range[0]
				rate_energy_range_max = rate_energy_range[1]
			except KeyError:
				rate_energy_range = None
				rate_energy_range_min = None
				rate_energy_range_max = None

			try:
				classification = record['classification']
				classification_grb = classification['GRB']
			except KeyError:
				classification = None

			try:
				far = record['far']
			except KeyError:
				far = None

		# 2 - guano.loc_map.example.json
		elif record_number == 2:

			try:
				healpix_file = record['healpix_file']
			except KeyError:
				healpix_file = None

			try:
				systematic_included = record['systematic_included']
			except KeyError:
				systematic_included = None

			try:
				rate_snr = record['rate_snr']
			except KeyError:
				rate_snr = None

			try:
				rate_duration = record['rate_duration']
			except KeyError:
				rate_duration = None

			try:
				rate_energy_range = record['rate_energy_range']
				rate_energy_range_min = rate_energy_range[0]
				rate_energy_range_max = rate_energy_range[1]
			except KeyError:
				rate_energy_range = None
				rate_energy_range_min = None
				rate_energy_range_max = None

			try:
				classification = record['classification']
				classification_grb = classification['GRB']
			except KeyError:
				classification = None

			try:
				far = record['far']
			except KeyError:
				far = None

		# 3 - guano.loc_arc_min.example.json
		elif record_number == 3:

			try:
				ra = record['ra']
			except KeyError:
				ra = None

			try:
				dec = record['dec']
			except KeyError:
				dec = None

			try:
				ra_dec_error = record['ra_dec_error']
			except KeyError:
				ra_dec_error = None

			try:
				containment_probability = record['containment_probability']
			except KeyError:
				containment_probability = None

			try:
				systematic_included = record['systematic_included']
			except KeyError:
				systematic_included = None

			try:
				rate_snr = record['rate_snr']
			except KeyError:
				rate_snr = None

			try:
				rate_duration = record['rate_duration']
			except KeyError:
				rate_duration = None

			try:
				rate_energy_range = record['rate_energy_range']
				rate_energy_range_min = rate_energy_range[0]
				rate_energy_range_max = rate_energy_range[1]
			except KeyError:
				rate_energy_range = None
				rate_energy_range_min = None
				rate_energy_range_max = None

			try:
				classification = record['classification']
				classification_grb = classification['GRB']
			except KeyError:
				classification = None

			try:
				far = record['far']
			except KeyError:
				far = None

		# 4 - guano.reraction.example.json
		elif record_number == 4:
			pass

		else:
			pass


	@staticmethod
	def igwn_gwalert(value, topic):

		record = json.loads(value)
		print(record)