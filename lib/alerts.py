from config import Configuration
from email.message import EmailMessage
import json
import requests
import smtplib

class Alerts:

	@staticmethod
	def filter_alert(value, topic):

		# GCN Kafka Alerts
		if topic == 'gcn.circulars':
			print('\nGCN Kafka Alert (' + str(topic) + ')')
			Alerts.gcn_circulars(value, topic)
			print()

		elif topic == 'gcn.heartbeat':
			pass
			#Alerts.gcn_heartbeat(value, topic)

		elif topic == 'gcn.notices.chime.frb':
			print('\nGCN Kafka Alert (' + str(topic) + ')')
			Alerts.gcn_notices_chime_frb(value, topic)
			print()

		elif topic == 'gcn.notices.dsa110.frb':
			print('GCN Kafka Alert (' + str(topic) + ')')
			Alerts.gcn_notices_dsa110_frb(value, topic)

		elif topic == 'gcn.notices.einstein_probe.wxt.alert':
			print('GCN Kafka Alert (' + str(topic) + ')')
			Alerts.gcn_notices_einstein_probe_wxt_alert(value, topic)

		elif topic == 'gcn.notices.icecube.gold_bronze_track_alerts':
			print('GCN Kafka Alert (' + str(topic) + ')')
			#Alerts.gcn_notices_icecube_gold_bronze_track_alerts(value, topic)

		elif topic == 'gcn.notices.icecube.lvk_nu_track_search':
			print('GCN Kafka Alert (' + str(topic) + ')')
			Alerts.gcn_notices_icecube_lvk_nu_track_search(value, topic)

		elif topic =='gcn.notices.superk.sn_alert':
			print('GCN Kafka Alert (' + str(topic) + ')')
			Alerts.gcn_notices_superk_sn_alert(value, topic)

		elif topic == 'gcn.notices.swift.bat.guano':
			print('GCN Kafka Alert (' + str(topic) + ')')
			Alerts.gcn_notices_swift_bat_guano(value, topic)

		elif topic == 'igwn.gwalert':
			print('GCN Kafka Alert (' + str(topic) + ')')
			Alerts.igwn_gwalert(value, topic)

		# Unknown Alerts
		else:
			print()
			print(value)
			print()

	@staticmethod
	def gcn_circulars(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN Circulars via GCN Kafka.

		:parameter value - The value of the message
		:parameter topic - The topic of the message
		:parameter alert - A toggle for broadcasting the alert

		:return - Nothing is returned
		'''

		record = json.loads(value)

		try:
			record_schema = record['$schema']
		except KeyError:
			record_schema = None

		try:
			record_event_id = record['eventId']
		except KeyError:
			record_event_id = None

		try:
			record_submitter = record['submitter']
		except KeyError:
			record_submitter = None

		try:
			record_submitted_how = record['submittedHow']
		except KeyError:
			record_submitted_how = None

		try:
			record_subject = record['subject']
		except KeyError:
			record_subject = None

		try:
			record_circular_id = record['circularId']
		except KeyError:
			record_circular_id = None

		try:
			record_format = record['format']
		except KeyError:
			record_format = None

		try:
			record_body = record['body']
		except KeyError:
			record_body = None

		try:
			record_created_on = record['createdOn']
		except KeyError:
			record_created_on = None

		print('\tSchema:', record_schema)
		print('\tEvent ID:', record_event_id)
		print('\tSubmitter:', record_submitter)
		print('\tSubmitted How:', record_submitted_how)
		print('\tSubject:', record_subject)
		print('\tCircular ID:', record_circular_id)
		print('\tFormat:', record_format)
		#print('\tBody:', record_body)
		print('\tCreated On:', record_created_on)

		if alert:
			body = str(record_circular_id) + '\n' + str(record_body)
			Alerts.send_alert(topic, body)

	@staticmethod
	def send_alert(topic, body):
		'''This function creates and sends an alert given a specific topic and body of text.

		:parameter topic - The topic string of the alert source
		:parameter body - The body of text describing the alert contents

		:return - Nothing is returned
		'''

		server = smtplib.SMTP_SSL(Configuration.SMTP, Configuration.PORT)
		server.ehlo()
		server.login(Configuration.EMAIL, Configuration.PAS)

		message = EmailMessage()
		message['Subject'] = 'TOROS Alert: ' + str(topic)
		message['From'] = Configuration.EMAIL
		message['To'] = ', '.join(Configuration.MAILING_LIST)
		message.set_content(body)

		server.send_message(message)
		server.quit()

	@staticmethod
	def gcn_heartbeat(value, topic):

		record = json.loads(value)

	@staticmethod
	def gcn_notices_chime_frb(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN Notices from the Canadian Hydrogen Intensity Mapping Experiment (CHIME) via GCN Kafka.

		:parameter value - The value of the message
		:parameter topic - The topic of the message
		:parameter alert - A toggle for broadcasting the alert

		:return - Nothing is returned
		'''

		record = json.loads(value)

		try:
			record_schema = record['$schema']
		except KeyError:
			record_schema = None

		try:
			record_alert_type = record['alert_type']
		except KeyError:
			record_alert_type = None

		# initial
		if record_alert_type == 'initial':
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			try:
				record_snr = record['snr']
			except KeyError:
				record_snr = None

			try:
				record_ra = record['ra']
			except KeyError:
				record_ra = None

			try:
				record_dec = record['dec']
			except KeyError:
				record_dec = None

			try:
				record_ra_dec_error = record['ra_dec_error']
			except KeyError:
				record_ra_dec_error = None

			try:
				record_dm = record['dm']
			except KeyError:
				record_dm = None

			try:
				record_dm_error = record['dm_error']
			except KeyError:
				record_dm_error = None

			try:
				record_dm_gal_ne_2001_max = record['dm_gal_ne_2001_max']
			except KeyError:
				record_dm_gal_ne_2001_max = None

			try:
				record_trigger_time_inf_freq = record['trigger_time_inf_freq']
			except KeyError:
				record_trigger_time_inf_freq = None

			try:
				record_trigger_time_inf_freq_error = record['trigger_time_inf_freq_error']
			except KeyError:
				record_trigger_time_inf_freq_error = None

			try:
				record_importance = record['importance']
			except KeyError:
				record_importance = None

			try:
				record_sampling_time = record['sampling_time']
			except KeyError:
				record_sampling_time = None

			try:
				record_spectral_band = record['spectral_band']
			except KeyError:
				record_spectral_band = None

			try:
				record_spectral_band_units = record['spectral_band_units']
			except KeyError:
				record_spectral_band_units = None

			try:
				record_npol = record['npol']
			except KeyError:
				record_npol = None

			try:
				record_tsys = record['tsys']
			except KeyError:
				record_tsys = None

			try:
				record_description = record['description']
			except KeyError:
				record_description = None

			# print the alert contents
			print('\tSchema:', record_schema)
			print('\tAlert type:', record_alert_type)
			print('\tID:', record_id)
			print('\tTrigger time:', record_trigger_time)
			print('\tTrigger time error (1-sigma):', record_trigger_time_error, 's')
			print('\tSNR:', record_snr)
			print('\tRA:', record_ra, 'deg')
			print('\tDec:', record_dec, 'deg')
			print('\tRA/Dec error:', record_ra_dec_error, 'deg')
			print('\tDM:', record_dm, 'pc/cm^3')
			print('\tDM error:', record_dm_error, 'pc/cm^3')
			print('\tDM (galactic, NE2001):', record_dm_gal_ne_2001_max, 'pc/cm^3')
			print('\tTrigger time (inf freq):', record_trigger_time_inf_freq)
			print('\tTrigger time (inf freq) error:', record_trigger_time_inf_freq_error, 's')
			print('\tImportance:', record_importance)
			print('\tSampling time:', record_sampling_time, 'ms')
			print('\tSpectral band:', record_spectral_band, record_spectral_band_units)
			print('\tNumber of polarizations:', record_npol)
			print('\tSystem temperature:', record_tsys, 'K')
			print('\tDescription:', record_description)

			# send an alert message
			if alert:
				body = str(record_id) + '\nRA: ' + str(record_ra) + ' deg\nDec: ' + str(record_dec) + ' deg'
				Alerts.send_alert(topic, body)

		# retraction
		elif record_alert_type == 'retraction':
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			try:
				record_description = record['description']
			except KeyError:
				record_description = None

			# print the alert contents
			print('\tSchema:', record_schema)
			print('\tAlert type:', record_alert_type)
			print('\tID:', record_id)
			print('\tTrigger time:', record_trigger_time)
			print('\tTrigger time error (1-sigma):', record_trigger_time_error, 's')
			print('\tDescription:', record_description)

			# send an alert message
			if alert:
				body = str(record_id) + '\nRetraction'
				Alerts.send_alert(topic, body)

		# subsequent
		elif record_alert_type == 'subsequent':
			try:
				record_known_source = record['known_source']
			except KeyError:
				record_known_source = None

			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			try:
				record_snr = record['snr']
			except KeyError:
				record_snr = None

			try:
				record_ra = record['ra']
			except KeyError:
				record_ra = None

			try:
				record_dec = record['dec']
			except KeyError:
				record_dec = None

			try:
				record_ra_dec_error = record['ra_dec_error']
			except KeyError:
				record_ra_dec_error = None

			try:
				record_dm = record['dm']
			except KeyError:
				record_dm = None

			try:
				record_dm_error = record['dm_error']
			except KeyError:
				record_dm_error = None

			try:
				record_dm_gal_ne_2001_max = record['dm_gal_ne_2001_max']
			except KeyError:
				record_dm_gal_ne_2001_max = None

			try:
				record_trigger_time_inf_freq = record['trigger_time_inf_freq']
			except KeyError:
				record_trigger_time_inf_freq = None

			try:
				record_trigger_time_inf_freq_error = record['trigger_time_inf_freq_error']
			except KeyError:
				record_trigger_time_inf_freq_error = None

			try:
				record_importance = record['importance']
			except KeyError:
				record_importance = None

			try:
				record_association_probability = record['association_probability']
			except KeyError:
				record_association_probability = None

			try:
				record_sampling_time = record['sampling_time']
			except KeyError:
				record_sampling_time = None

			try:
				record_spectral_band = record['spectral_band']
			except KeyError:
				record_spectral_band = None

			try:
				record_spectral_band_units = record['spectral_band_units']
			except KeyError:
				record_spectral_band_units = None

			try:
				record_npol = record['npol']
			except KeyError:
				record_npol = None

			try:
				record_tsys = record['tsys']
			except KeyError:
				record_tsys = None

			try:
				record_description = record['description']
			except KeyError:
				record_description = None

			# print the alert contents
			print('\tSchema:', record_schema)
			print('\tAlert type:', record_alert_type)
			print('\tID:', record_id)
			print('\tKnown source:', record_known_source)
			print('\tTrigger time:', record_trigger_time)
			print('\tTrigger time error (1-sigma):', record_trigger_time_error, 's')
			print('\tSNR:', record_snr)
			print('\tRA:', record_ra, 'deg')
			print('\tDec:', record_dec, 'deg')
			print('\tRA/Dec error:', record_ra_dec_error, 'deg')
			print('\tDM:', record_dm, 'pc/cm^3')
			print('\tDM error:', record_dm_error, 'pc/cm^3')
			print('\tDM (galactic, NE2001):', record_dm_gal_ne_2001_max, 'pc/cm^3')
			print('\tTrigger time (inf freq):', record_trigger_time_inf_freq)
			print('\tTrigger time (inf freq) error:', record_trigger_time_inf_freq_error, 's')
			print('\tImportance:', record_importance)
			print('\tAssociation probability:', record_association_probability)
			print('\tSampling time:', record_sampling_time, 'ms')
			print('\tSpectral band:', record_spectral_band, record_spectral_band_units)
			print('\tNumber of polarizations:', record_npol)
			print('\tSystem temperature:', record_tsys, 'K')
			print('\tDescription:', record_description)

			# send an alert message
			if alert:
				body = str(record_id) + '\nRA: ' + str(record_ra) + ' deg\nDec: ' + str(record_dec) + ' deg'
				Alerts.send_alert(topic, body)

		# update
		elif record_alert_type == 'update':
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			try:
				record_update_message = record['update_message']
			except KeyError:
				record_update_message = None

			try:
				record_description = record['description']
			except KeyError:
				record_description = None

			# print the alert contents
			print('\tSchema:', record_schema)
			print('\tAlert type:', record_alert_type)
			print('\tID:', record_id)
			print('\tTrigger time:', record_trigger_time)
			print('\tTrigger time error (1-sigma):', record_trigger_time_error, 's')
			print('\tUpdate message:', record_update_message)
			print('\tDescription:', record_description)

			# send an alert message
			if alert:
				body = str(record_id) + '\n' + str(record_update_message)
				Alerts.send_alert(topic, body)

	@staticmethod
	def gcn_notices_dsa110_frb(value, topic):

		record = json.loads(value)

		try:
			record_schema = record['$schema']
		except KeyError:
			record_schema = None

		try:
			record_alert_type = record['alert_type']
		except KeyError:
			record_alert_type = None

		if record_alert_type == 'initial':
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			try:
				record_snr = record['snr']
			except KeyError:
				record_snr = None

			try:
				record_dm = record['dm']
			except KeyError:
				record_dm = None

			try:
				record_event_duration = record['event_duration']
			except KeyError:
				record_event_duration = None

			try:
				record_ra = record['ra']
			except KeyError:
				record_ra = None

			try:
				record_dec = record['dec']
			except KeyError:
				record_dec = None

			try:
				record_ra_dec_error = record['ra_dec_error']
			except KeyError:
				record_ra_dec_error = None

			try:
				record_importance = record['importance']
			except KeyError:
				record_importance = None

		elif record_alert_type == 'retraction':
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			try:
				record_description = record['description']
			except KeyError:
				record_description = None

		elif record_alert_type == 'subsequent':
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			try:
				record_known_source = record['known_source']
			except KeyError:
				record_known_source = None

			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			try:
				record_snr = record['snr']
			except KeyError:
				record_snr = None

			try:
				record_dm = record['dm']
			except KeyError:
				record_dm = None

			try:
				record_event_duration = record['event_duration']
			except KeyError:
				record_event_duration = None

			try:
				record_ra = record['ra']
			except KeyError:
				record_ra = None

			try:
				record_dec = record['dec']
			except KeyError:
				record_dec = None

			try:
				record_ra_dec_error = record['ra_dec_error']
			except KeyError:
				record_ra_dec_error = None

			try:
				record_importance = record['importance']
			except KeyError:
				record_importance = None

		elif record_alert_type == 'update':
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			try:
				record_description = record['description']
			except KeyError:
				record_description = None

	@staticmethod
	def gcn_notices_einstein_probe_wxt_alert(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN Notices from the Einstein Probe Wide Field X-ray Telescope (EP-WXT) via GCN Kafka.

		:parameter value - The value of the message
		:parameter topic - The topic of the message
		:parameter alert - A toggle for broadcasting the alert

		:return - Nothing is returned
		'''

		record = json.loads(value)
		
		try:
			record_schema = record['$schema']
		except KeyError:
			record_schema = None

		try:
			record_instrument = record['instrument']
		except KeyError:
			record_instrument = None

		try:
			record_trigger_time = record['trigger_time']
		except KeyError:
			record_trigger_time = None

		try:
			record_id = record['id'][0]
		except KeyError:
			record_id = None

		try:
			record_ra = record['ra']
		except KeyError:
			record_ra = None

		try:
			record_dec = record['dec']
		except KeyError:
			record_dec = None

		try:
			record_ra_dec_error = record['ra_dec_error']
		except KeyError:
			record_ra_dec_error = None

		try:
			record_image_energy_range = record['image_energy_range']
		except KeyError:
			record_image_energy_range = None

		try:
			record_net_count_rate = record['net_count_rate']
		except KeyError:
			record_net_count_rate = None

		try:
			record_image_snr = record['image_snr']
		except KeyError:
			record_image_snr = None

		try:
			record_additional_info = record['additional_info']
		except KeyError:
			record_additional_info = None

		# print the alert contents
		print('\t$schema:', record_schema)
		print('\tInstrument:', record_instrument)
		print('\tTrigger time:', record_trigger_time)
		print('\tID:', record_id)
		print('\tRA:', record_ra)
		print('\tDec:', record_dec)
		print('\tRA/Dec error:', record_ra_dec_error)
		print('\tImage energy range:', record_image_energy_range)
		print('\tNet count rate:', record_net_count_rate)
		print('\tImage SNR:', record_image_snr)
		print('\tAdditional info:', record_additional_info)

		# send an alert message
		if alert:
			body = str(record_id) + '\n' + str(record_ra) + '\n' + str(record_dec)
			Alerts.send_alert(topic, body)

	@staticmethod
	def gcn_notices_icecube_lvk_nu_track_search(value, topic):

		record = json.loads(value)

		try:
			record_schema = record['$schema']
		except KeyError:
			record_schema = 'unknown: $schema'

		try:
			record_type = record['type']
		except KeyError:
			record_type = 'unknown: type'

		try:
			record_reference = record['reference']
			record_gcn_notices_lvk_alert = record_reference['gcn.notices.LVK.alert']
		except KeyError:
			record_reference = {'gcn.notices.LVK.alert':'unknown:'}
			record_gcn_notices_lvk_alert = record_reference['gcn.notices.LVK.alert']

		try:
			record_ref_id = record['ref_ID']
		except KeyError:
			record_ref_id = 'unknown: ref_ID'

		try:
			record_alert_datetime = record['alert_datetime']
		except KeyError:
			record_alert_datetime = 'unknown: alert_datetime'

		try:
			record_trigger_time = record['trigger_time']
		except KeyError:
			record_trigger_time = 'unknown: trigger_time'

		try:
			record_observation_start = record['observation_start']
		except KeyError:
			record_observation_start = 'unknown: observation_start'

		try:
			record_observation_stop = record['observation_stop']
		except KeyError:
			record_observation_stop = 'unknown: observation_stop'

		try:
			record_observation_livetime = record['observation_livetime']
		except KeyError:
			record_observation_livetime = -999.

		try:
			record_pval_generic = record['pval_generic']
		except KeyError:
			record_pval_generic = -999.

		try:
			record_pval_bayesian = record['pval_bayesian']
		except KeyError:
			record_pval_bayesian = -999.

		try:
			record_n_events_coincident = record['n_events_coincident']
			if record_n_events_coincident > 0:
				try:
					record_coincident_events = record['coincident_events']
					for i in range(record_n_events_coincident):
						record_coincident_event = record_coincident_events[i]
						try:
							record_event_dt = record_coincident_event['event_dt']
						except KeyError:
							record_event_dt = -999.
						try:
							record_localization = record_coincident_event['localization']
							try:
								record_ra = record_localization['ra']
							except KeyError:
								record_ra = -999.
							try:
								record_dec = record_localization['dec']
							except KeyError:
								record_dec = -999.
							try:
								record_ra_dec_error = record_localization['ra_dec_error']
							except KeyError:
								record_ra_dec_error = -999.
							try:
								record_containment_probability = record_localization['containment_probability']
							except KeyError:
								record_containment_probability = -999.
							try:
								record_systematic_included = record_localization['systematic_included']
							except KeyError:
								record_systematic_included = False
						except KeyError:
							record_localization = {}
						try:
							record_id = record_coincident_event['id'][0]
						except KeyError:
							record_id = 'unknown: id'
						try:
							record_event_pval_generic = record_coincident_event['event_pval_generic']
						except KeyError:
							record_event_pval_generic = -999.
						try:
							record_event_pval_bayesian = record_coincident_event['event_pval_bayesian']
						except KeyError:
							record_event_pval_bayesian = -999.
				except KeyError:
					record_coincident_events = None
		except KeyError:
			record_n_events_coincident = None

		try:
			record_most_probable_direction = record['most_probable_direction']
			record_most_probable_direction_ra = most_probable_direction['ra']
			record_most_probable_direction_dec = most_probable_direction['dec']
		except KeyError:
			record_most_probable_direction = {'ra': -999., 'dec': -999.}
			record_most_probable_direction_ra = record_most_probable_direction['ra']
			record_most_probable_direction_dec = record_most_probable_direction['dec']
		try:
			record_neutrino_flux_sensitivity_range = record['neutrino_flux_sensitivity_range']
			try:
				record_flux_sensitivity = record_neutrino_flux_sensitivity_range['flux_sensitivity']
				record_flux_sensitivity_min = record_flux_sensitivity[0]
				record_flux_sensitivity_max = record_flux_sensitivity[1]
			except KeyError:
				record_flux_sensitivity = [-999., -999.]
				record_flux_sensitivity_min = record_flux_sensitivity[0]
				record_flux_sensitivity_max = record_flux_sensitivity[1]

			try:
				record_sensitive_energy_range = record_neutrino_flux_sensitivity_range['sensitive_energy_range']
				record_sensitive_energy_range_min = record_sensitive_energy_range[0]
				record_sensitive_energy_range_max = record_sensitive_energy_range[1]
			except KeyError:
				record_sensitive_energy_range = [-999., -999.]
				record_sensitive_energy_range_min = record_sensitive_energy_range[0]
				record_sensitive_energy_range_max = record_sensitive_energy_range[1]

		except KeyError:
			record_neutrino_flux_sensitivity_range = {}

	@staticmethod
	def gcn_notices_superk_sn_alert(value, topic):

		record = json.loads(value)
		
		try:
			record_schema = record['$schema']
		except KeyError:
			record_schema = None

		try:
			record_mission = record['mission']
		except KeyError:
			record_mission = None

		try:
			record_messenger = record['messenger']
		except KeyError:
			record_messenger = None

		try:
			record_id = record['id']
		except KeyError:
			record_id = None

		try:
			record_record_number = record['record_number']
		except KeyError:
			record_record_number = None

		try:
			record_trigger_number = record['trigger_number']
		except KeyError:
			record_trigger_number = None

		try:
			record_alert_datetime = record['alert_datetime']
		except KeyError:
			record_alert_datetime = None

		try:
			record_alert_tense = record['current']
		except KeyError:
			record_alert_tense = None

		try:
			record_alert_type = record['alert_type']
		except KeyError:
			record_alert_type = None

		try:
			record_trigger_time = record['trigger_time']
		except KeyError:
			record_trigger_time = None

		try:
			record_processed_sample = record['processed_sample']
		except KeyError:
			record_processed_sample = None

		try:
			record_pipeline = record['pipeline']
		except KeyError:
			record_pipeline = None

		try:
			record_n_events = record['n_events']
		except KeyError:
			record_n_events = None

		try:
			record_n_ibd_events = record['n_ibd_events']
		except KeyError:
			record_n_ibd_events = None

		try:
			record_detection_interval = record['detection_interval']
		except KeyError:
			record_detection_interval = None

		try:
			record_rate_energy_range = record['rate_energy_range']
		except KeyError:
			record_rate_energy_range = None

		try:
			record_ra = record['ra']
		except KeyError:
			record_ra = None

		try:
			record_dec = record['dec']
		except KeyError:
			record_dec = None

		try:
			record_ra_dec_error = record['ra_dec_error']
		except KeyError:
			reord_ra_dec_error = None

		try:
			record_containment_probability = record['containment_probability']
		except KeyError:
			record_containment_probability = None

		try:
			record_luminosity_distance = record['luminosity_distance']
		except KeyError:
			record_luminosity_distance = None

		try:
			record_luminosity_distance_error = record['luminosity_distance_error']
		except KeyError:
			record_luminosity_distance_error = None

	@staticmethod
	def gcn_notices_swift_bat_guano(value, topic):

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