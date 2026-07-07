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
			Alerts.gcn_circulars(value, topic, False)
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
			print()

		elif topic == 'gcn.notices.einstein_probe.wxt.alert':
			print('GCN Kafka Alert (' + str(topic) + ')')
			Alerts.gcn_notices_einstein_probe_wxt_alert(value, topic)
			print()

		elif topic == 'gcn.notices.icecube.gold_bronze_track_alerts':
			print('GCN Kafka Alert (' + str(topic) + ')')
			Alerts.gcn_notices_icecube_gold_bronze_track_alerts(value, topic)
			print()

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

		The Canadian Hydrogen Intensity Mapping Experiment (CHIME) is a transit radio telescope located in Penticton, BC, Canada. It consists of four semi-cylindrical reflectors, each with 2000 square meters of collecting area and 256 dual-polarization antennas. The telescope observes between 400 MHz and 800 MHz and covers an instantaneous field of view of ~200 square degrees. CHIME houses several electronic backends, which are tailored for specific scientific goals, such as generating cosmological maps of hydrogen density, detecting fast radio bursts (FRBs), and observing and timing pulsars. In particular, the FRB backend operates at ~1 ms time resolution and ~0.4 MHz frequency resolution.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
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
			# trigger_time [string] - time of trigger [ISO 8601]
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			# trigger_time_error [number|array] - trigger time uncertainty with symmetric uncertainty; array if uncertainty is asymmetric [s, 1σ]
			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			# id [string|array] - instrument-specific trigger ID or alternate ID; array if multiple instrument-specific trigger IDs are present
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			# snr [number] - signal-to-noise ratio of the burst
			try:
				record_snr = record['snr']
			except KeyError:
				record_snr = None

			# ra [number] - ICRS right ascension; utilizes the J2000 epoch and an equatorial coordinate system [deg]
			try:
				record_ra = record['ra']
			except KeyError:
				record_ra = None

			# deg [number] - ICRS declination; utilizes the J2000 epoch and an equatorial coordinate system [deg]
			try:
				record_dec = record['dec']
			except KeyError:
				record_dec = None

			# ra_dec_error - ?
			try:
				record_ra_dec_error = record['ra_dec_error']
			except KeyError:
				record_ra_dec_error = None

			# dm [number] - dispersion measure (DM) of the burst; represents the integrated column density of free electrons along the line of sight [pc/cm^3]
			try:
				record_dm = record['dm']
			except KeyError:
				record_dm = None

			# dm_error [number|array] - uncertainty associated with the dispersion measure with symmetric uncertainty; array if uncertainty is asymmetric [pc/cm^3, 1σ]
			try:
				record_dm_error = record['dm_error']
			except KeyError:
				record_dm_error = None

			# dm_gal_ne_2001_max [number] - estimated contribution to the dispersion measure from the Galaxy using the NE2001 model [pc/cm^3]
			try:
				record_dm_gal_ne_2001_max = record['dm_gal_ne_2001_max']
			except KeyError:
				record_dm_gal_ne_2001_max = None

			# trigger_time_inf_freq [string] - time of the trigger at infinite frequency [ISO 8601]
			try:
				record_trigger_time_inf_freq = record['trigger_time_inf_freq']
			except KeyError:
				record_trigger_time_inf_freq = None

			# trigger_time_inf_freq_error [number] - error on the trigger time at infinite frequency [s]
			try:
				record_trigger_time_inf_freq_error = record['trigger_time_inf_freq_error']
			except KeyError:
				record_trigger_time_inf_freq_error = None

			# importance [number] - a machine learning score separating RFI (0) from an astrophysical signal (1)
			try:
				record_importance = record['importance']
			except KeyError:
				record_importance = None

			# sampling_time [number] - time resolution of the FRB search at the host observatory [ms]
			try:
				record_sampling_time = record['sampling_time']
			except KeyError:
				record_sampling_time = None

			# spectral_band [array] - observed spectral band, expressed in specific 'units' field
			try:
				record_spectral_band = record['spectral_band']
			except KeyError:
				record_spectral_band = None

			# spectral_band_units [string] - units for the spectral data; default unit is keV; options: keV, nm, MHz
			try:
				record_spectral_band_units = record['spectral_band_units']
			except KeyError:
				record_spectral_band_units = None

			# npol [number] - the number of polarizations of the real-time FRB search
			try:
				record_npol = record['npol']
			except KeyError:
				record_npol = None

			# tsys [number] - the system temperature of the real-time FRB search [K]
			try:
				record_tsys = record['tsys']
			except KeyError:
				record_tsys = None

			# description [string] - message description
			try:
				record_description = record['description']
			except KeyError:
				record_description = None

			# print the alert contents
			print('\tSchema:', record_schema)
			print('\tAlert type:', record_alert_type)
			print('\tID:', record_id)
			print('\tTrigger time:', record_trigger_time)
			print('\tTrigger time error (1σ):', record_trigger_time_error, 's')
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
				ln_01 = 'Schema: ' + str(record_schema) + '\n'
				ln_02 = 'Alert type: ' + str(record_alert_type) + '\n'
				ln_03 = 'ID: ' + str(record_id) + '\n'
				ln_04 = 'Trigger time: ' + str(record_trigger_time) + '\n'
				ln_05 = 'Trigger time error (1σ): ' + str(record_trigger_time_error) + ' s\n'
				ln_06 = 'SNR: ' + str(record_snr) + '\n'
				ln_07 = 'RA: ' + str(record_ra) + ' deg\n'
				ln_08 = 'Dec: ' + str(record_dec) + ' deg\n'
				ln_09 = 'RA/Dec error: ' + str(record_ra_dec_error) + ' deg\n'
				ln_10 = 'DM: ' + str(record_dm) + ' pc/cm^3\n'
				ln_11 = 'DM error: ' + str(record_dm_error) + ' pc/cm^3\n'
				ln_12 = 'DM (galactic, NE2001): ' + str(record_dm_gal_ne_2001_max) + ' pc/cm^3\n'
				ln_13 = 'Trigger time (inf freq): ' + str(record_trigger_time_inf_freq)
				ln_14 = 'Trigger time (inf freq) error: ' + str(record_trigger_time_inf_freq_error) + ' s\n'
				ln_15 = 'Importance: ' + str(record_importance) + '\n'
				ln_16 = 'Sampling time: ' + str(record_sampling_time) + ' ms\n'
				ln_17 = 'Spectral band: ' + str(record_spectral_band) + ' ' + str(record_spectral_band_units) + '\n'
				ln_18 = 'Number of polarizations: ' + str(record_npol) + '\n'
				ln_19 = 'System temperature: ' + str(record_tsys) + ' K\n'
				ln_20 = 'Description: ' + str(record_description) + '\n'
				body = ln_01 + ln_02 + ln_03 + ln_04 + ln_05 + ln_06 + ln_07 + ln_08 + ln_09 + ln_10 + ln_11 + ln_12 + ln_13 + ln_14 + ln_15 + ln_16 + ln_17 + ln_18 + ln_19 + ln_20
				Alerts.send_alert(topic, body)

		# retraction
		elif record_alert_type == 'retraction':
			# id [string|array] - instrument-specific trigger ID or alternate ID; array if multiple instrument-specific trigger IDs are present
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			# trigger_time [string] - time of trigger [ISO 8601]
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			# trigger_time_error [number|array] - trigger time uncertainty with symmetric uncertainty; array if uncertainty is asymmetric [s, 1σ]
			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			# description [string] - message description
			try:
				record_description = record['description']
			except KeyError:
				record_description = None

			# print the alert contents
			print('\tSchema:', record_schema)
			print('\tAlert type:', record_alert_type)
			print('\tID:', record_id)
			print('\tTrigger time:', record_trigger_time)
			print('\tTrigger time error (1σ):', record_trigger_time_error, 's')
			print('\tDescription:', record_description)

			# send an alert message
			if alert:
				ln_01 = 'Schema: ' + str(record_schema) + '\n'
				ln_02 = 'Alert type: ' + str(record_alert_type) + '\n'
				ln_03 = 'ID: ' + str(record_id) + '\n'
				ln_04 = 'Trigger time: ' + str(record_trigger_time) + '\n'
				ln_05 = 'Trigger time error (1σ): ' + str(record_trigger_time_error) + ' s\n'
				ln_06 = 'Description: ' + str(record_description) + '\n'
				body = ln_01 + ln_02 + ln_03 + ln_04 + ln_05 + ln_06
				Alerts.send_alert(topic, body)

		# subsequent
		elif record_alert_type == 'subsequent':
			# known_source [string] - Transient Name Server name of the known source (if a repeating FRB source)
			try:
				record_known_source = record['known_source']
			except KeyError:
				record_known_source = None

			# trigger_time [string] - time of trigger [ISO 8601]
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			# trigger_time_error [number|array] - trigger time uncertainty with symmetric uncertainty; array if uncertainty is asymmetric [s, 1σ]
			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			# id [string|array] - instrument-specific trigger ID or alternate ID; array if multiple instrument-specific trigger IDs are present
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			# snr [number] - signal-to-noise ratio of the burst
			try:
				record_snr = record['snr']
			except KeyError:
				record_snr = None

			# ra [number] - ICRS right ascension; utilizes the J2000 epoch and equatorial coordinate system [deg]
			try:
				record_ra = record['ra']
			except KeyError:
				record_ra = None

			# dec [number] - ICRS declination; utilizes the J2000 epoch and equatorial coordinate system [deg]
			try:
				record_dec = record['dec']
			except KeyError:
				record_dec = None

			# ra_dec_error - ?
			try:
				record_ra_dec_error = record['ra_dec_error']
			except KeyError:
				record_ra_dec_error = None

			# dm [number] - dispersion measure (DM) of the burst; represents the integrated column density of free electrons along the line of sight [pc/cm^3]
			try:
				record_dm = record['dm']
			except KeyError:
				record_dm = None

			# dm_error [number|array] - uncertainty associated with the dispersion measure with symmetric uncertainty; array if uncertainty is asymmetric [pc/cm^3, 1σ]
			try:
				record_dm_error = record['dm_error']
			except KeyError:
				record_dm_error = None

			# dm_gal_ne_2001_max [number] - estimated contribution to the dispersion measure from the Galaxy using the NE2001 model [pc/cm^3]
			try:
				record_dm_gal_ne_2001_max = record['dm_gal_ne_2001_max']
			except KeyError:
				record_dm_gal_ne_2001_max = None

			# trigger_time_inf_freq [string] - time of the trigger at infinite frequency [ISO 8601]
			try:
				record_trigger_time_inf_freq = record['trigger_time_inf_freq']
			except KeyError:
				record_trigger_time_inf_freq = None

			# trigger_time_inf_freq_error [number] - error on the trigger time at infinite frequency [s]
			try:
				record_trigger_time_inf_freq_error = record['trigger_time_inf_freq_error']
			except KeyError:
				record_trigger_time_inf_freq_error = None

			# importance [number] - a machine learning score separating RFI (0) from an astrophysical signal (1)
			try:
				record_importance = record['importance']
			except KeyError:
				record_importance = None

			# association_probability [number] - a score of known source association probability from poorly associated (0) to confidently associated (1)
			try:
				record_association_probability = record['association_probability']
			except KeyError:
				record_association_probability = None

			# sampling_time [number] - time resolution of the FRB search at the host observatory [ms]
			try:
				record_sampling_time = record['sampling_time']
			except KeyError:
				record_sampling_time = None

			# spectral_band [array] - observed spectral band, expressed in specific 'units' field
			try:
				record_spectral_band = record['spectral_band']
			except KeyError:
				record_spectral_band = None

			# spectral_band_units [string] - units for the spectral data; default unit is keV; options: keV, nm, MHz
			try:
				record_spectral_band_units = record['spectral_band_units']
			except KeyError:
				record_spectral_band_units = None

			# npol [number] - the number of polarizations of the real-time FRB search
			try:
				record_npol = record['npol']
			except KeyError:
				record_npol = None

			# tsys [number] - the system temperature of the real-time FRB search [K]
			try:
				record_tsys = record['tsys']
			except KeyError:
				record_tsys = None

			# description [string] - message description
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
			print('\tTrigger time error (1σ):', record_trigger_time_error, 's')
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
				ln_01 = 'Schema: ' + str(record_schema) + '\n'
				ln_02 = 'Alert type: ' + str(record_alert_type) + '\n'
				ln_03 = 'ID: ' + str(record_id) + '\n'
				ln_04 = 'Known source: ' + str(record_known_source) + '\n'
				ln_05 = 'Trigger time: ' + str(record_trigger_time) + '\n'
				ln_06 = 'Trigger time error (1σ): ' + str(record_trigger_time_error) + ' s\n'
				ln_07 = 'SNR: ' + str(record_snr) + '\n'
				ln_08 = 'RA: ' + str(record_ra) + ' deg\n'
				ln_09 = 'Dec: ' + str(record_dec) + ' deg\n'
				ln_10 = 'RA/Dec error: ' + str(record_ra_dec_error) + ' deg\n'
				ln_11 = 'DM: ' + str(record_dm) + ' pc/cm^3\n'
				ln_12 = 'DM error: ' + str(record_dm_error) + ' pc/cm^3\n'
				ln_13 = 'DM (galactic, NE2001): ' + str(record_dm_gal_ne_2001_max) + ' pc/cm^3\n'
				ln_14 = 'Trigger time (inf freq): ' + str(record_trigger_time_inf_freq) + '\n'
				ln_15 = 'Trigger time (inf freq) error: ' + str(record_trigger_time_inf_freq_error) + ' s\n'
				ln_16 = 'Importance: ' + str(record_importance) + '\n'
				ln_17 = 'Association probability: ' + str(record_association_probability) + '\n'
				ln_18 = 'Sampling time: ' + str(record_sampling_time) + ' ms\n'
				ln_19 = 'Spectral band: ' + str(record_spectral_band) + ' ' + str(record_spectral_band_units) + '\n'
				ln_20 = 'Number of polarizations: ' + str(record_npol) + '\n'
				ln_21 = 'System temperature: ' + str(record_tsys), ' K\n'
				ln_22 = 'Description: ' + str(record_description) + '\n'
				body = ln_01 + ln_02 + ln_03 + ln_04 + ln_05 + ln_06 + ln_07 + ln_08 + ln_09 + ln_10 + ln_11 + ln_12 + ln_13 + ln_14 + ln_15 + ln_16 + ln_17 + ln_18 + ln_19 + ln_20 + ln_21 + ln_22
				Alerts.send_alert(topic, body)

		# update
		elif record_alert_type == 'update':
			# id [string|array] - instrument-specific trigger ID or alternate ID; array if multiple instrument-specific trigger IDs are present
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			# trigger_time [string] - time of trigger [ISO 8601]
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			# trigger_time_error [number|array] - trigger time uncertainty with symmetric uncertainty; array if uncertainty is asymmetric [s, 1σ]
			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			# update_message [string] - message to be included in an update alert
			try:
				record_update_message = record['update_message']
			except KeyError:
				record_update_message = None

			# description [string] - message description
			try:
				record_description = record['description']
			except KeyError:
				record_description = None

			# print the alert contents
			print('\tSchema:', record_schema)
			print('\tAlert type:', record_alert_type)
			print('\tID:', record_id)
			print('\tTrigger time:', record_trigger_time)
			print('\tTrigger time error (1σ):', record_trigger_time_error, 's')
			print('\tUpdate message:', record_update_message)
			print('\tDescription:', record_description)

			# send an alert message
			if alert:
				ln_01 = 'Schema: ' + str(record_schema) + '\n'
				ln_02 = 'Alert type: ' + str(record_alert_type) + '\n'
				ln_03 = 'ID: ' + str(record_id) + '\n'
				ln_04 = 'Trigger time: ' + str(record_trigger_time) + '\n'
				ln_05 = 'Trigger time error (1σ): ' + str(record_trigger_time_error) + ' s\n'
				ln_06 = 'Update message: ' + str(record_update_message) + '\n'
				ln_07 = 'Description: ' + str(record_description) + '\n'
				body = ln_01 + ln_02 + ln_03 + ln_04 + ln_05 + ln_06 + ln_07
				Alerts.send_alert(topic, body)

	@staticmethod
	def gcn_notices_dsa110_frb(value, topic):
		'''This function processes GCN Notices from the Deep Synoptic Array-110 (DSA-110) via GCN Kafka.

		The Deep Synoptic Array-110 (DSA-110) is a radio interferometer purpose-built for fast radio burst (FRB) detection and direct localization. The array is located at the Owens Valley Radio Observatory (OVRO) comprised of 96 4.65-m dishes that continuously survey for FRBs at frequencies between 1280 and 1530 MHz. Over a three-year science program, the DSA-110 will deliver a sample of more than 300 FRBs, each localized to regions ~1 arcminute radius within 1 minute of detection. This is made possible by a suite of novel instrumentation, including precisely engineered antennas, ultra-low noise ambient temperature receivers, and a powerful real-time, autonomous data-analysis system.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
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
			# trigger time [string] - time of trigger [ISO 8601]
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			# id [string|array] - instrument-specific trigger ID or alternate ID; array if multiple instrument-specific trigger IDs are present
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			# snr [number] - signal-to-noise ratio of the burst
			try:
				record_snr = record['snr']
			except KeyError:
				record_snr = None

			# dm [number] - dispersion measure (DM) of the burst; represents the integrated column density of free electrons along the line of sight [pc/cm^3]
			try:
				record_dm = record['dm']
			except KeyError:
				record_dm = None

			# event_duration [number] - time duration of the event [ms]
			try:
				record_event_duration = record['event_duration']
			except KeyError:
				record_event_duration = None

			# ra [number] - ICRS right ascension; utilizes the J2000 epoch and an equatorial coordinate system [deg]
			try:
				record_ra = record['ra']
			except KeyError:
				record_ra = None

			# deg [number] - ICRS declination; utilizes the J2000 epoch and an equatorial coordinate system [deg]
			try:
				record_dec = record['dec']
			except KeyError:
				record_dec = None

			# ra_dec_error - ?
			try:
				record_ra_dec_error = record['ra_dec_error']
			except KeyError:
				record_ra_dec_error = None

			# importance [number] - a machine learning score separating RFI (0) from an astrophysical signal (1)
			try:
				record_importance = record['importance']
			except KeyError:
				record_importance = None

			# print the alert contents
			print('\tSchema:', record_schema)
			print('\tAlert type:', record_alert_type)
			print('\tID:', record_id)
			print('\tSNR:', record_snr)
			print('\tDM:', record_dm, 'pc/cm^3')
			print('\tEvent duration:', record_event_duration, 'ms')
			print('\tRA:', record_ra, 'deg')
			print('\tDec:', record_dec, 'deg')
			print('\tRA/Dec error:', record_ra_dec_error, 'deg')
			print('\tImportance:', record_importance)

			# send an alert message
			if alert:
				ln_01 = 'Schema: ' + str(record_schema) + '\n'
				ln_02 = 'Alert type: ' + str(record_alert_type) + '\n'
				ln_03 = 'ID: ' + str(record_id) + '\n'
				ln_04 = 'SNR: ' + str(record_snr) + '\n'
				ln_05 = 'DM: ' + str(record_dm) + ' pc/cm^3\n'
				ln_06 = 'Event duration: ' + str(record_event_duration) + ' ms\n'
				ln_07 = 'RA: ' + str(record_ra) + ' deg\n'
				ln_08 = 'Dec: ' + str(record_dec) + ' deg\n'
				ln_09 = 'RA/Dec error: ' + str(record_ra_dec_error) + ' deg\n'
				ln_10 = 'Importance: ' + str(record_importance) + '\n'
				body = ln_01 + ln_02 + ln_03 + ln_04 + ln_05 + ln_06 + ln_07 + ln_08 + ln_09 + ln_10
				Alerts.send_alert(topic, body)

		# retraction
		elif record_alert_type == 'retraction':
			# id [string|array] - instrument-specific trigger ID or alternate ID; array if multiple instrument-specific trigger IDs are present
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			# trigger_time [string] - time of trigger [ISO 8601]
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			# trigger_time_error [number|array] - trigger time uncertainty with symmetric uncertainty; array if uncertainty is asymmetric [s, 1σ]
			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			# description [string] - description of how the event was identified and what it represents
			try:
				record_description = record['description']
			except KeyError:
				record_description = None

			# print the alert contents
			print('\tSchema:', record_schema)
			print('\tAlert type:', record_alert_type)
			print('\tID:', record_id)
			print('\tTrigger time:', record_trigger_time)
			print('\tTrigger time error (1σ):', record_trigger_time_error, 's')
			print('\tDescription:', record_description)

			# send an alert message
			if alert:
				ln_01 = 'Schema: ' + str(record_schema) + '\n'
				ln_02 = 'Alert type: ' + str(record_alert_type) + '\n'
				ln_03 = 'ID: ' + str(record_id) + '\n'
				ln_04 = 'Trigger time: ' + str(record_trigger_time) + '\n'
				ln_05 = 'Trigger time error (1σ): ' + str(record_trigger_time_error) + ' s\n'
				ln_06 = 'Description: ' + str(record_description) + '\n'
				body = ln_01 + ln_02 + ln_03 + ln_04 + ln_05 + ln_06
				Alerts.send_alert(topic, body)

		# subsequent
		elif record_alert_type == 'subsequent':
			# trigger_time [string] - time of trigger [ISO 8601]
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			# trigger_time_error [number|array] - trigger time uncertainty with symmetric uncertainty; array if uncertainty is asymmetric [s, 1σ]
			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			# known_source [string] - name of a known FRB that this event has been associated to
			try:
				record_known_source = record['known_source']
			except KeyError:
				record_known_source = None

			# id [string|array] - instrument-specific trigger ID or alternate ID; array if multiple instrument-specific trigger IDs are present
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			# snr [number] - signal-to-noise ratio of the burst
			try:
				record_snr = record['snr']
			except KeyError:
				record_snr = None

			# dm [number] - dispersion measure (DM) of the burst; represents the integrated column density of free electrons along the line of sight [pc/cm^3]
			try:
				record_dm = record['dm']
			except KeyError:
				record_dm = None

			# event_duration [number] - time duration of the event [ms]
			try:
				record_event_duration = record['event_duration']
			except KeyError:
				record_event_duration = None

			# ra [number] - ICRS right ascension; utilizes the J2000 epoch and an equatorial coordinate system [deg]
			try:
				record_ra = record['ra']
			except KeyError:
				record_ra = None

			# dec [number] - ICRS declination; utilizes the J2000 epoch and an equatorial coordinate system [deg]
			try:
				record_dec = record['dec']
			except KeyError:
				record_dec = None

			# ra_dec_error - ?
			try:
				record_ra_dec_error = record['ra_dec_error']
			except KeyError:
				record_ra_dec_error = None

			# importance [number] - a machine learning score separating RFI (0) from an astrophysical signal (1)
			try:
				record_importance = record['importance']
			except KeyError:
				record_importance = None

			# print the alert contents
			print('\tSchema:', record_schema)
			print('\tAlert type:', record_alert_type)
			print('\tID:', record_id)
			print('\tKnown source:', record_known_source)
			print('\tTrigger time:', record_trigger_time)
			print('\tTrigger time error (1σ):', record_trigger_time_error, 's')
			print('\tSNR:', record_snr)
			print('\tDM:', record_dm, 'pc/cm^3')
			print('\tEvent duration:', record_event_duration, 'ms')
			print('\tRA:', record_ra, 'deg')
			print('\tDec:', record_dec, 'deg')
			print('\tRA/Dec error:', record_ra_dec_error, 'deg')
			print('\tImportance:', record_importance)

			# send an alert message
			if alert:
				ln_01 = 'Schema: ' + str(record_schema) + '\n'
				ln_02 = 'Alert type: ' + str(record_alert_type) + '\n'
				ln_03 = 'ID: ' + str(record_id) + '\n'
				ln_04 = 'Known source: ' + str(record_known_source) + '\n'
				ln_05 = 'Trigger time: ' + str(record_trigger_time) + '\n'
				ln_06 = 'Trigger time error (1σ): ' + str(record_trigger_time_error) + ' s\n'
				ln_07 = 'SNR: ' + str(record_snr) + '\n'
				ln_08 = 'DM: ' + str(record_dm) + ' pc/cm^3\n'
				ln_09 = 'Event duration: ' + str(record_event_duration) + ' ms\n'
				ln_10 = 'RA: ' + str(record_ra) + ' deg\n'
				ln_11 = 'Dec: ' + str(record_dec) + ' deg\n'
				ln_12 = 'RA/Dec error: ' + str(record_ra_dec_error) + ' deg\n'
				ln_13 = 'Imporance: ' + str(record_importance) + '\n'
				body = ln_01 + ln_02 + ln_03 + ln_04 + ln_05 + ln_06 + ln_07 + ln_08 + ln_09 + ln_10 + ln_11 + ln_12 + ln_13
				Alerts.send_alert(topic, body)

		# update
		elif record_alert_type == 'update':
			# id [string|array] - instrument-specific trigger ID or alternate ID; array if multiple instrument-specific trigger IDs are present
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			# trigger_time [string] - time of trigger [ISO 8601]
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			# trigger_time_error [number|array] - trigger time uncertainty with symmetric uncertainty; array if uncertainty is asymmetric [s, 1σ]
			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			# description [string] - description of how the event was identified and what it represents
			try:
				record_description = record['description']
			except KeyError:
				record_description = None

			# print the alert contents
			print('\tSchema:', record_schema)
			print('\tAlert type:', record_alert_type)
			print('\tID:', record_id)
			print('\tTrigger time:', record_trigger_time)
			print('\tTrigger time error (1σ):', record_trigger_time_error, 's')
			print('\tDescription:', record_description)

			# send an alert message
			if alert:
				ln_01 = 'Schema: ' + str(record_schema) + '\n'
				ln_02 = 'Alert type: ' + str(record_alert_type) + '\n'
				ln_03 = 'ID: ' + str(record_id) + '\n'
				ln_04 = 'Trigger time: ' + str(record_trigger_time) + '\n'
				ln_05 = 'Trigger time error (1σ): ' + str(record_trigger_time_error) + ' s\n'
				ln_06 = 'Description: ' + str(record_description) + '\n'
				body = ln_01 + ln_02 + ln_03 + ln_04 + ln_05 + ln_06
				Alerts.send_alert(topic, body)

	@staticmethod
	def gcn_notices_einstein_probe_wxt_alert(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN Notices from the Einstein Probe Wide Field X-ray Telescope (EP-WXT) via GCN Kafka.

		The Einstein Probe (EP) is a mission of the Chinese Academy of Sciences (CAS), in collaboration with the European Space Agency (ESA) and the Max Planck Institute for Extraterrestrial Physics (MPE), Germany, dedicated to time-domain high-energy astrophysics. Its primary goals are to discover high-energy transients and monitor variable objects. To achieve this, EP employs a very large instantaneous field-of-view (3600 square degrees), along with moderate spatial resolution (FWHM ~5 arcminutes) and energy resolution in the 0.5-5 keV band. EP has also the capability of performing fast and deep follow-up observations in the 0.5-10 keV energy band (effective area 2*300 cm^2 at 1 keV; half-power diameter, HPD, ~30 arcseconds), as well as quick downlink of transient alert messages.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		record = json.loads(value)
		
		# $schema - ?
		try:
			record_schema = record['$schema']
		except KeyError:
			record_schema = None

		# instrument [string] - name of the instrument reporting the event
		try:
			record_instrument = record['instrument']
		except KeyError:
			record_instrument = None

		# trigger_time [string] - time of trigger [ISO 8601]
		try:
			record_trigger_time = record['trigger_time']
		except KeyError:
			record_trigger_time = None

		# id [string|array] - instrument-specific trigger ID or alternate ID; array if multiple instrument-specific trigger IDs are present
		try:
			record_id = record['id'][0]
		except KeyError:
			record_id = None

		# ra [number] - ICRS right ascension; utilizes the J2000 epoch and an equatorial coordinate system [deg]
		try:
			record_ra = record['ra']
		except KeyError:
			record_ra = None

		# dec [number] - ICRS declination; utilizes the J2000 epoch and an equatorial coordinate system [deg]
		try:
			record_dec = record['dec']
		except KeyError:
			record_dec = None

		# ra_dec_error - ?
		try:
			record_ra_dec_error = record['ra_dec_error']
		except KeyError:
			record_ra_dec_error = None

		# image_energy_range [array] - low and high energy bounds used in image signal-to-noise ratio calculation [keV]
		try:
			record_image_energy_range = record['image_energy_range']
		except KeyError:
			record_image_energy_range = None

		# net_count_rate [number] - net count rate of the transient above the background over rate duration and rate energy range [counts/s]
		try:
			record_net_count_rate = record['net_count_rate']
		except KeyError:
			record_net_count_rate = None

		# image_snr [number] - image signal-to-noise ratio
		try:
			record_image_snr = record['image_snr']
		except KeyError:
			record_image_snr = None

		# additional_info [string] - additional information about the event
		try:
			record_additional_info = record['additional_info']
		except KeyError:
			record_additional_info = None

		# print the alert contents
		print('\tSchema:', record_schema)
		print('\tInstrument:', record_instrument)
		print('\tTrigger time:', record_trigger_time)
		print('\tID:', record_id)
		print('\tRA:', record_ra, 'deg')
		print('\tDec:', record_dec, 'deg')
		print('\tRA/Dec error:', record_ra_dec_error, 'deg')
		print('\tImage energy range:', record_image_energy_range, 'keV')
		print('\tNet count rate:', record_net_count_rate, 'counts/s')
		print('\tImage SNR:', record_image_snr)
		print('\tAdditional info:', record_additional_info)

		# send an alert message
		if alert:
			ln_01 = 'Schema: ' + str(record_schema) + '\n'
			ln_02 = 'Instrument: ' + str(record_instrument) + '\n'
			ln_03 = 'Trigger time: ' + str(record_trigger_time) + '\n'
			ln_04 = 'ID: ' + str(record_id) + '\n'
			ln_05 = 'RA: ' + str(record_ra) + ' deg\n'
			ln_06 = 'Dec: ' + str(record_dec) + ' deg\n'
			ln_07 = 'RA/Dec error: ' + str(record_ra_dec_error) + ' deg\n'
			ln_08 = 'Image energy range: ' + str(record_image_energy_range) + ' keV\n'
			ln_09 = 'Net count rate: ' + record_net_count_rate, ' counts/s\n'
			ln_10 = 'Image SNR: ' + str(record_image_snr) + '\n'
			ln_11 = 'Additional info: ' + str(record_additional_info) + '\n'
			body = ln_01 + ln_02 + ln_03 + ln_04 + ln_05 + ln_06 + ln_07 + ln_08 + ln_09 + ln_10 + ln_10 + ln_11
			Alerts.send_alert(topic, body)

	@staticmethod
	def gcn_notices_icecube_gold_bronze_track_alerts(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN Notices (high-energy single-track events) from the IceCube Neutrino Observatory via GCN Kafka.

		The IceCube Neutrino Observatory is a cubic-kilometer Cherenkov particle detector deployed in the Antarctic ice beneath the Amundsen-Scott South Pole Station. It consists of 86 strings of photo-detectors, extending to a depth of about 2,500 meters below the glacier's surface and instrumenting a cubic-kilometer of ice. The Digital Optical Module photo-detectors detect the light produced by relativistic charged particles produced by neutrino interactions in or near the instrumented volume of ice.

		IceCube is sensitive to neutrinos from all directions. As neutrinos pass through the ice, their interactions can leave track signatures (order-km in length) in the IceCube detector array when they produce secondary muons or compact signatures (cascades of ~m in extent) when they produce secondary electrons or hadrons. Track events can be reconstructed with an uncertainty of less than 1 degree, while cascade events have higher signal purity. IceCube triggers on signals for neutrino energies greater than ~10 GeV (10^10 eV) and can identify likely astrophysical neutrino events for neutrino energies greater than ~50 TeV (5 x 10^13 eV) and issue alerts to the community.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		record = json.loads(value)

		# $schema - GCN schema structure according to the nasa-gcn/gcn-schema GitHub project
		try:
			record_schema = record['$schema']
		except KeyError:
			record_schema = None

		# mission [string] - instrument issuing the alert (i.e. 'IceCube')
		try:
			record_mission = record['mission']
		except KeyError:
			record_mission = None

		# instrument [string] - detector configuration (e.g. 'IC86' refers to IceCube with 86 strings active, or full detector)
		try:
			record_instrument = record['instrument']
		except KeyError:
			record_instrument = None

		# messenger [string] - type of astrophysical messenger (i.e. 'Neutrino')
		try:
			record_messenger = record['messenger']
		except KeyError:
			record_messenger = None

		# pipeline [string] - alert type, either 'Bronze Track Alert' or 'Gold Track Alert'; bronze alerts have an average astrophysical probability of 30%; gold alerts have an average astrophysical probability of 50% (determined by the p_astro parameter)
		try:
			record_pipeline = record['pipeline']
		except KeyError:
			record_pipeline = None

		# record_number [number] - 0 for the first preliminary reconstruction; 1 for the revised reconstruction
		try:
			record_record_number = record['record_number']
		except KeyError:
			record_record_number = None

		# event_name [string] - name of the event (e.g. 'IceCube-YYMMDDA' where the letter increments if multiple alerts occur on the same day)
		try:
			record_event_name = record['event_name']
		except KeyError:
			record_event_name = None

		# id [string] - internal event ID, formatted as 'RUNID_EVENTID'
		try:
			record_id = record['id']
		except KeyError:
			record_id = None

		# alert_datetime [string] - timestamp of the GCN notice (UTC time)
		try:
			record_alert_datetime = record['alert_datetime']
		except KeyError:
			record_alert_datetime = None

		# alert_type [string] - 'initial' for preliminary reconstruction (Revision 0); 'update' for revised reconstruction (Revision 1); 'retraction' if the alert is retracted
		try:
			record_alert_type = record['alert_type']
		except KeyError:
			record_alert_type = None

		# alert_tense [string] - 'current' for real alerts; 'test' for tests; 'injections' for fake injections used as tests
		try:
			record_alert_tense = record['alert_tense']
		except KeyError:
			record_alert_tense = None

		# alert_topology [string] - event topology
		try:
			record_alert_topology = record['alert_topology']
		except KeyError:
			record_alert_topology = None

		# number_of_events [number] - number of events contributing to the alerts (1 being a single neutrino event)
		try:
			record_number_of_events = record['number_of_events']
		except KeyError:
			record_number_of_events = None

		# ra [number] - best-fit right ascension (J2000) [deg]
		try:
			record_ra = record['ra']
		except KeyError:
			record_ra = None

		# dec [number] - best-fit declination (J2000) [deg]
		try:
			record_dec = record['dec']
		except KeyError:
			record_dec = None

		# ra_dec_error [number] - combined positional uncertainty, calculated from the RA and Dec errors at 90% containment
		try:
			record_ra_dec_error = record['ra_dec_error']
		except KeyError:
			record_ra_dec_error = None

		# containment_probability [number] - confidence level for the positional error (e.g. '0.9' for 90%)
		try:
			record_containment_probability = record['containment_probability']
		except KeyError:
			record_containment_probability = None

		# systematic_included [boolean] - 'true' if systematics are included in the positional uncertainty, otherwise 'false'
		try:
			record_systematic_included = record['systematic_included']
		except KeyError:
			record_systematic_included = None

		# healpix_url [string] - URL to the multi-order probability map of the event; 'null' for preliminary reconstruction
		try:
			record_healpix_url = record['healpix_url']
		except KeyError:
			record_healpix_url = None

		# trigger_time [string] - time of the triggered neutrino event (UTC time)
		try:
			record_trigger_time = record['trigger_time']
		except KeyError:
			record_trigger_time = None

		# nu_energy [number] - estimated neutrino energy [TeV]
		try:
			record_nu_energy = record['nu_energy']
		except KeyError:
			record_nu_energy = None

		# p_astro [number] - estimated probability that the event is astrophysical (in the past, this was named 'signalness')
		try:
			record_p_astro = record['p_astro']
		except KeyError:
			record_p_astro = None

		# far [number] - false alarm rate (FAR), indicating the expected frequency of such events in IceCube [Hz]
		try:
			record_far = record['far']
		except KeyError:
			record_far = None

		# print the alert contents
		print('\tSchema:', record_schema)
		print('\tMission:', record_mission)
		print('\tInstrument:', record_instrument)
		print('\tMessenger:', record_messenger)
		print('\tPipeline:', record_pipeline)
		print('\tRecord number:', record_record_number)
		print('\tEvent name:', record_event_name)
		print('\tID:', record_id)
		print('\tAlert datetime:', record_alert_datetime)
		print('\tAlert type:', record_alert_type)
		print('\tAlert tense:', record_alert_tense)
		print('\tAlert topology:', record_alert_topology)
		print('\tNumber of events:', record_number_of_events)
		print('\tRA:', record_ra, 'deg')
		print('\tDec:', record_dec, 'deg')
		print('\tRA/Dec error:', record_ra_dec_error, 'deg')
		print('\tContainment probability:', record_containment_probability)
		print('\tSystematic included:', record_systematic_included)
		print('\tHEALPix URL:', record_healpix_url)
		print('\tTrigger time:', record_trigger_time)
		print('\tν energy:', record_nu_energy, 'TeV')
		print('\tAstrophysical probability:', record_p_astro)

		# send an alert message
		if alert:
			ln_01 = 'Schema: ' + str(record_schema) + '\n'
			ln_02 = 'Mission: ' + str(record_mission) + '\n'
			ln_03 = 'Instrument: ' + str(record_instrument) + '\n'
			ln_04 = 'Messenger: ' + str(record_messenger) + '\n'
			ln_05 = 'Pipeline: ' + str(record_pipeline) + '\n'
			ln_06 = 'Record number: ' + str(record_record_number) + '\n'
			ln_07 = 'Event name: ' + str(record_event_name) + '\n'
			ln_08 = 'ID: ' + str(record_id) + '\n'
			ln_09 = 'Alert datetime: ' + str(record_alert_datetime) + '\n'
			ln_10 = 'Alert type: ' + str(record_alert_type) + '\n'
			ln_11 = 'Alert tense: ' + str(record_alert_tense) + '\n'
			ln_12 = 'Alert topology: ' + str(record_alert_topology) + '\n'
			ln_13 = 'Number of events: ' + str(record_number_of_events) + '\n'
			ln_14 = 'RA: ' + str(record_ra) + ' deg\n'
			ln_15 = 'Dec: ' + str(record_dec) + ' deg\n'
			ln_16 = 'RA/Dec error: ' + str(record_ra_dec_error) + ' deg\n'
			ln_17 = 'Containment probability: ' + str(record_containment_probability) + '\n'
			ln_18 = 'Systematic included: ' + str(record_systematic_included) + '\n'
			ln_19 = 'HEALPix URL: ' + str(record_healpix_url) + '\n'
			ln_20 = 'Trigger time: ' + str(record_trigger_time) + '\n'
			ln_21 = 'ν energy: ' + str(record_nu_energy) + ' TeV\n'
			ln_22 = 'Astrophysical probability: ' + str(record_p_astro) + '\n'
			body = ln_01 + ln_02 + ln_03 + ln_04 + ln_05 + ln_06 + ln_07 + ln_08 + ln_09 + ln_10 + ln_11 + ln_12 + ln_13 + ln_14 + ln_15 + ln_16 + ln_17 + ln_18 + ln_19 + ln_20 + ln_21 + ln_22
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