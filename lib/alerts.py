from config import Configuration
from email.message import EmailMessage
import json
import requests
import smtplib

class Alerts:

	@staticmethod
	def chime_frb(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN Notices from the Canadian Hydrogen Intensity Mapping Experiment (CHIME) via GCN Kafka.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		record = json.loads(value)
		body = ''

		try:
			rec_schema = record['$schema']
			str_schema = '\tSchema: ' + str(rec_schema) + '\n'
		except KeyError:
			rec_schema = None
			str_schema = '\tSchema: none\n'
		body += str_schema
		print(str_schema)

		try:
			rec_alert_type = record['alert_type']
			str_alert_type = '\tAlert type: ' + str(rec_alert_type) + '\n'
		except KeyError:
			rec_alert_type = None
			str_alert_type = '\tAlert type: none\n'
		body += str_alert_type
		print(str_alert_type)

		# initial
		if rec_alert_type == 'initial':
			# trigger_time
			try:
				rec_trigger_time = record['trigger_time']
				str_trigger_time = '\tTrigger time: ' + str(rec_trigger_time) + '\n'
			except KeyError:
				rec_trigger_time = None
				str_trigger_time = '\tTrigger time: none\n'
			body += str_trigger_time
			print(str_trigger_time)

			# trigger_time_error
			try:
				rec_trigger_time_error = record['trigger_time_error']
				str_trigger_time_error = '\tTrigger time error: ' + str(rec_trigger_time_error) + ' s\n'
			except KeyError:
				rec_trigger_time_error = None
				str_trigger_time_error = '\tTrigger time error: none\n'
			body += str_trigger_time_error
			print(str_trigger_time_error)

			# id
			try:
				rec_id = record['id']
				str_id = '\tID: ' + str(rec_id) + '\n'
			except KeyError:
				rec_id = None
				str_id = '\tID: none\n'
			body += str_id
			print(str_id)

			# snr
			try:
				rec_snr = record['snr']
				str_snr = '\tSNR: ' + str(rec_snr) + '\n'
			except KeyError:
				rec_snr = None
				str_snr = '\tSNR: none\n'
			body += str_snr
			print(str_snr)

			# ra
			try:
				rec_ra = record['ra']
				str_ra = '\tRA: ' + str(rec_ra) + ' deg\n'
			except KeyError:
				rec_ra = None
				str_ra = '\tRA: none\n'
			body += str_ra
			print(str_ra)

			# dec
			try:
				rec_dec = record['dec']
				str_dec = '\tDec: ' + str(rec_dec) + ' deg\n'
			except KeyError:
				rec_dec = None
				str_dec = '\tDec: none\n'
			body += str_dec
			print(str_dec)

			# ra_dec_error
			try:
				rec_ra_dec_error = record['ra_dec_error']
				str_ra_dec_error = '\tRA/Dec error: ' + str(rec_ra_dec_error) + '\n'
			except KeyError:
				rec_ra_dec_error = None
				str_ra_dec_error = '\tRA/Dec error: none\n'
			body += str_ra_dec_error
			print(str_ra_dec_error)

			# dm
			try:
				rec_dm = record['dm']
				str_dm = '\tDM: ' + str(rec_dm) + ' pc/cm^3\n'
			except KeyError:
				rec_dm = None
				str_dm = '\tDM: none\n'
			body += str_dm
			print(str_dm)

			# dm_error
			try:
				rec_dm_error = record['dm_error']
				str_dm_error = '\tDM error: ' + str(rec_dm_error) + ' pc/cm^3\n'
			except KeyError:
				rec_dm_error = None
				str_dm_error = '\tDM error: none\n'
			body += str_dm_error
			print(str_dm_error)

			# dm_gal_ne_2001_max
			try:
				rec_dm_gal_ne_2001_max = record['dm_gal_ne_2001_max']
				str_dm_gal_ne_2001_max = '\tDM (galactic, NE2001): ' + str(rec_dm_gal_ne_2001_max) + ' pc/cm^3\n'
			except KeyError:
				rec_dm_gal_ne_2001_max = None
				str_dm_gal_ne_2001_max = '\tDM (galactic, NE2001): none\n'
			body += str_dm_gal_ne_2001_max
			print(str_dm_gal_ne_2001_max)

			# trigger_time_inf_freq
			try:
				rec_trigger_time_inf_freq = record['trigger_time_inf_freq']
				str_trigger_time_inf_freq = '\tTrigger time (inf freq): ' + str(rec_trigger_time_inf_freq) + '\n'
			except KeyError:
				rec_trigger_time_inf_freq = None
				str_trigger_time_inf_freq = '\tTrigger time (inf freq): none\n'
			body += str_trigger_time_inf_freq
			print(str_trigger_time_inf_freq)

			# trigger_time_inf_freq_error
			try:
				rec_trigger_time_inf_freq_error = record['trigger_time_inf_freq_error']
				str_trigger_time_inf_freq_error = '\tTrigger time (inf freq) error: ' + str(rec_trigger_time_inf_freq_error) + ' s\n'
			except KeyError:
				rec_trigger_time_inf_freq_error = None
				str_trigger_time_inf_freq_error = '\tTrigger time (inf freq) error: none\n'
			body += str_trigger_time_inf_freq_error
			print(str_trigger_time_inf_freq_error)

			# importance
			try:
				rec_importance = record['importance']
				str_importance = '\tImportance: ' + str(rec_importance) + '\n'
			except KeyError:
				rec_importance = None
				str_importance = '\tImportance: none\n'
			body += str_importance
			print(str_importance)

			# sampling_time
			try:
				rec_sampling_time = record['sampling_time']
				str_sampling_time = '\tSampling time: ' + str(rec_sampling_time) + ' ms\n'
			except KeyError:
				rec_sampling_time = None
				str_sampling_time = '\tSampling time: none\n'
			body += str_sampling_time
			print(str_sampling_time)

			# spectral_band
			try:
				rec_spectral_band = record['spectral_band']
				str_spectral_band = '\tSpectral band: ' + str(rec_spectral_band) + '\n'
			except KeyError:
				rec_spectral_band = None
				str_spectral_band = '\tSpectral band: none\n'
			body += str_spectral_band
			print(str_spectral_band)

			# spectral_band_units
			try:
				rec_spectral_band_units = record['spectral_band_units']
				str_spectral_band_units = '\tSpectral band units: ' + str(rec_spectral_band_units) + '\n'
			except KeyError:
				rec_spectral_band_units = None
				str_spectral_band_units = '\tSpectral band units: none\n'
			body += str_spectral_band_units
			print(str_spectral_band_units)

			# npol
			try:
				rec_npol = record['npol']
				str_npol = '\tNumber of polarizations: ' + str(rec_npol) + '\n'
			except KeyError:
				rec_npol = None
				str_npol = '\tNumber of polarizations: none\n'
			body += str_npol
			print(str_npol)

			# tsys
			try:
				rec_tsys = record['tsys']
				str_tsys = '\tSystem temperature: ' + str(rec_tsys) + ' K\n'
			except KeyError:
				rec_tsys = None
				str_tsys = '\tSystem temperature: none\n'
			body += str_tsys
			print(str_tsys)

			# description
			try:
				rec_description = record['description']
				str_description = '\tDescription: ' + str(rec_description) + '\n'
			except KeyError:
				rec_description = None
				str_description = '\tDescription: none\n'
			body += str_description
			print(str_description)

			# send an alert message
			if alert:
				Alerts.send_alert(topic, body)

		# retraction
		elif rec_alert_type == 'retraction':
			# id
			try:
				rec_id = record['id']
				str_id = '\tID: ' + str(rec_id) + '\n'
			except KeyError:
				rec_id = None
				str_id = '\tID: none\n'
			body += str_id
			print(str_id)

			# trigger_time
			try:
				rec_trigger_time = record['trigger_time']
				str_trigger_time = '\tTrigger time: ' + str(rec_trigger_time) + '\n'
			except KeyError:
				rec_trigger_time = None
				str_trigger_time = '\tTrigger time: none\n'
			body += str_trigger_time
			print(str_trigger_time)

			# trigger_time_error
			try:
				rec_trigger_time_error = record['trigger_time_error']
				str_trigger_time_error = '\tTrigger time error: ' + str(trigger_time_error) + ' s\n'
			except KeyError:
				rec_trigger_time_error = None
				str_trigger_time_error = '\tTrigger time error: none\n'
			body += str_trigger_time_error
			print(str_trigger_time_error)

			# description
			try:
				rec_description = record['description']
				str_description = '\tDescription: ' + str(rec_description) + '\n'
			except KeyError:
				rec_description = None
				str_description = '\tDescription: none\n'
			body += str_description
			print(str_description)

			# send an alert message
			if alert:
				Alerts.send_alert(topic, body)

		# subsequent
		elif rec_alert_type == 'subsequent':
			# known_source
			try:
				rec_known_source = record['known_source']
				str_known_source = '\tKnown source: ' + str(rec_known_source) + '\n'
			except KeyError:
				rec_known_source = None
				str_known_source = '\tKnown source: none\n'
			body += str_known_source
			print(str_known_source)

			# trigger_time
			try:
				rec_trigger_time = record['trigger_time']
				str_trigger_time = '\tTrigger time: ' + str(rec_trigger_time) + '\n'
			except KeyError:
				rec_trigger_time = None
				str_trigger_time = '\tTrigger time: none\n'
			body += str_trigger_time
			print(str_trigger_time)

			# trigger_time_error
			try:
				rec_trigger_time_error = record['trigger_time_error']
				str_trigger_time_error = '\tTrigger time error: ' + str(rec_trigger_time_error) + ' s\n'
			except KeyError:
				rec_trigger_time_error = None
				str_trigger_time_error = '\tTrigger time error: none\n'
			body += str_trigger_time_error
			print(str_trigger_time_error)

			# id
			try:
				rec_id = record['id']
				str_id = '\tID: ' + str(rec_id) + '\n'
			except KeyError:
				rec_id = None
				str_id = '\tID: none\n'
			body += str_id
			print(str_id)

			# snr
			try:
				rec_snr = record['snr']
				str_snr = '\tSNR: ' + str(rec_snr) + '\n'
			except KeyError:
				rec_snr = None
				str_snr = '\tSNR: none\n'
			body += str_snr
			print(str_snr)

			# ra
			try:
				rec_ra = record['ra']
				str_ra = '\tRA: ' + str(rec_ra) + ' deg\n'
			except KeyError:
				rec_ra = None
				str_ra = '\tRA: none\n'
			body += str_ra
			print(str_ra)

			# dec
			try:
				rec_dec = record['dec']
				str_dec = '\tDec: ' + str(rec_dec) + ' deg\n'
			except KeyError:
				rec_dec = None
				str_dec = '\tDec: none\n'
			body += str_dec
			print(str_dec)

			# ra_dec_error - ?
			try:
				rec_ra_dec_error = record['ra_dec_error']
				str_ra_dec_error = '\tRA/Dec error: ' + str(rec_ra_dec_error) + '\n'
			except KeyError:
				rec_ra_dec_error = None
				str_ra_dec_error = '\tRA/Dec error: none\n'
			body += str_ra_dec_error
			print(str_ra_dec_error)

			# dm
			try:
				rec_dm = record['dm']
				str_dm = '\tDM: ' + str(rec_dm) + ' pc/cm^3\n'
			except KeyError:
				rec_dm = None
				str_dm = '\tDM: none\n'
			body += str_dm
			print(str_dm)

			# dm_error
			try:
				rec_dm_error = record['dm_error']
				str_dm_error = '\tDM error: ' + str(rec_dm_error) + ' pc/cm^3\n'
			except KeyError:
				rec_dm_error = None
				str_dm_error = '\tDM error: none\n'
			body += str_dm_error
			print(str_dm_error)

			# dm_gal_ne_2001_max
			try:
				rec_dm_gal_ne_2001_max = record['dm_gal_ne_2001_max']
				str_dm_gal_ne_2001_max = '\tDM (galactic, NE2001): ' + str(rec_dm_gal_ne_2001_max) + ' pc/cm^3\n'
			except KeyError:
				rec_dm_gal_ne_2001_max = None
				str_dm_gal_ne_2001_max = '\tDM (galactic, NE2001): none\n'
			body += str_dm_gal_ne_2001_max
			print(str_dm_gal_ne_2001_max)

			# trigger_time_inf_freq
			try:
				rec_trigger_time_inf_freq = record['trigger_time_inf_freq']
				str_trigger_time_inf_freq = '\tTrigger time (inf freq): ' + str(rec_trigger_time_inf_freq) + '\n'
			except KeyError:
				rec_trigger_time_inf_freq = None
				str_trigger_time_inf_freq = '\tTrigger time (inf freq) error: none\n'
			body += str_trigger_time_inf_freq
			print(str_trigger_time_inf_freq)

			# trigger_time_inf_freq_error
			try:
				rec_trigger_time_inf_freq_error = record['trigger_time_inf_freq_error']
				str_trigger_time_inf_freq_error = '\tTrigger time (inf freq) error: ' + str(rec_trigger_time_inf_freq_error) + ' s\n'
			except KeyError:
				rec_trigger_time_inf_freq_error = None
				str_trigger_time_inf_freq_error = '\tTrigger time (inf freq) error: none\n'
			body += str_trigger_time_inf_freq_error
			print(str_trigger_time_inf_freq_error)

			# importance
			try:
				rec_importance = record['importance']
				str_importance = '\tImportance: ' + str(rec_importance) + '\n'
			except KeyError:
				rec_importance = None
				str_importance = '\tImportance: none\n'
			body += str_importance
			print(str_importance)

			# association_probability
			try:
				rec_association_probability = record['association_probability']
				str_association_probability = '\tAssociation probability: ' + str(rec_association_probability) + '\n'
			except KeyError:
				rec_association_probability = None
				str_association_probability = '\tAssociation probability: none\n'
			body += str_association_probability
			print(str_association_probability)

			# sampling_time
			try:
				rec_sampling_time = record['sampling_time']
				str_sampling_time = '\tSampling time: ' + str(rec_sampling_time) + ' ms\n'
			except KeyError:
				rec_sampling_time = None
				str_sampling_time = '\tSampling time: none\n'
			body += str_sampling_time
			print(str_sampling_time)

			# spectral_band
			try:
				rec_spectral_band = record['spectral_band']
				str_spectral_band = '\tSpectral band: ' + str(rec_spectral_band) + '\n'
			except KeyError:
				rec_spectral_band = None
				str_spectral_band = '\tSpectral band: none\n'
			body += str_spectral_band
			print(str_spectral_band)

			# spectral_band_units
			try:
				rec_spectral_band_units = record['spectral_band_units']
				str_spectral_band_units = '\tSpectral band units: ' + str(rec_spectral_band_units) + '\n'
			except KeyError:
				rec_spectral_band_units = None

			# npol
			try:
				record_npol = record['npol']
			except KeyError:
				record_npol = None

			# tsys
			try:
				record_tsys = record['tsys']
			except KeyError:
				record_tsys = None

			# description
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
			# id
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			# trigger_time
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			# trigger_time_error
			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			# update_message
			try:
				record_update_message = record['update_message']
			except KeyError:
				record_update_message = None

			# description
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
	def decode_classic_notice(value):
		'''This function decodes the value of a classic text GCN notice when received from GCN Kafka.

		:parameter value - the value of the message from the GCN Kafka consumer

		:return params - a dictionary of the parameters contained in the message value
		'''

		contents = value.decode('utf-8').split('\n')

		keys = []
		items = []

		for content in contents:
			if len(content.split(': ')) > 1:
				key = content.split(': ')[0].strip()
				keys.append(key)

				item = content.split(': ')[1].strip()
				items.append(item)

		params = dict(zip(keys, items))

		return params

	@staticmethod
	def dsa110_frb(value, topic):
		'''This function processes GCN Notices from the Deep Synoptic Array-110 (DSA-110) via GCN Kafka.

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
			# trigger time
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			# id
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			# snr
			try:
				record_snr = record['snr']
			except KeyError:
				record_snr = None

			# dm
			try:
				record_dm = record['dm']
			except KeyError:
				record_dm = None

			# event_duration
			try:
				record_event_duration = record['event_duration']
			except KeyError:
				record_event_duration = None

			# ra
			try:
				record_ra = record['ra']
			except KeyError:
				record_ra = None

			# deg
			try:
				record_dec = record['dec']
			except KeyError:
				record_dec = None

			# ra_dec_error
			try:
				record_ra_dec_error = record['ra_dec_error']
			except KeyError:
				record_ra_dec_error = None

			# importance
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
			# id
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			# trigger_time
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			# trigger_time_error
			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			# description
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
			# trigger_time
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			# trigger_time_error
			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			# known_source
			try:
				record_known_source = record['known_source']
			except KeyError:
				record_known_source = None

			# id
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			# snr
			try:
				record_snr = record['snr']
			except KeyError:
				record_snr = None

			# dm
			try:
				record_dm = record['dm']
			except KeyError:
				record_dm = None

			# event_duration
			try:
				record_event_duration = record['event_duration']
			except KeyError:
				record_event_duration = None

			# ra
			try:
				record_ra = record['ra']
			except KeyError:
				record_ra = None

			# dec
			try:
				record_dec = record['dec']
			except KeyError:
				record_dec = None

			# ra_dec_error
			try:
				record_ra_dec_error = record['ra_dec_error']
			except KeyError:
				record_ra_dec_error = None

			# importance
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
			# id
			try:
				record_id = record['id']
			except KeyError:
				record_id = None

			# trigger_time
			try:
				record_trigger_time = record['trigger_time']
			except KeyError:
				record_trigger_time = None

			# trigger_time_error
			try:
				record_trigger_time_error = record['trigger_time_error']
			except KeyError:
				record_trigger_time_error = None

			# description
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
	def einstein_probe_wxt_alert(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN Notices from the Einstein Probe Wide Field X-ray Telescope (EP-WXT) via GCN Kafka.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		record = json.loads(value)
		
		# $schema
		try:
			record_schema = record['$schema']
		except KeyError:
			record_schema = None

		# instrument
		try:
			record_instrument = record['instrument']
		except KeyError:
			record_instrument = None

		# trigger_time
		try:
			record_trigger_time = record['trigger_time']
		except KeyError:
			record_trigger_time = None

		# id
		try:
			record_id = record['id'][0]
		except KeyError:
			record_id = None

		# ra
		try:
			record_ra = record['ra']
		except KeyError:
			record_ra = None

		# dec
		try:
			record_dec = record['dec']
		except KeyError:
			record_dec = None

		# ra_dec_error
		try:
			record_ra_dec_error = record['ra_dec_error']
		except KeyError:
			record_ra_dec_error = None

		# image_energy_range
		try:
			record_image_energy_range = record['image_energy_range']
		except KeyError:
			record_image_energy_range = None

		# net_count_rate
		try:
			record_net_count_rate = record['net_count_rate']
		except KeyError:
			record_net_count_rate = None

		# image_snr
		try:
			record_image_snr = record['image_snr']
		except KeyError:
			record_image_snr = None

		# additional_info
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
			ln_09 = 'Net count rate: ' + str(record_net_count_rate), ' counts/s\n'
			ln_10 = 'Image SNR: ' + str(record_image_snr) + '\n'
			ln_11 = 'Additional info: ' + str(record_additional_info) + '\n'
			body = ln_01 + ln_02 + ln_03 + ln_04 + ln_05 + ln_06 + ln_07 + ln_08 + ln_09 + ln_10 + ln_10 + ln_11
			Alerts.send_alert(topic, body)

	@staticmethod
	def fermi_gbm_alert(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN classic text Fermi GBM Alert Notices (gcn.classic.text.FERMI_GBM_ALERT) from the Fermi Gamma-ray Space Telescope via GCN Kafka.

		GBM Alert Notices occur when the GBM first triggers. They do not contain a RA/Dec location of the burst; only a date, time, trigger criteria, trigger detection significance, and the algorithm used to make the detection.

		:parameter value - the value of the message
		;parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		params = Alerts.decode_classic_notice(value)

		rec_title = params['TITLE']
		rec_notice_date = params['NOTICE_DATE']
		rec_notice_type = params['NOTICE_TYPE']
		rec_record_num = params['RECORD_NUM']
		rec_trigger_num = params['TRIGGER_NUM']
		rec_grb_date = params['GRB_DATE']
		rec_grb_time = params['GRB_TIME']
		rec_trigger_signif = params['TRIGGER_SIGNIF']
		rec_trigger_dur = params['TRIGGER_DUR']
		rec_e_range = params['E_RANGE']
		rec_algorithm = params['ALGORITHM']
		rec_detectors = params['DETECTORS']
		rec_comments = params['COMMENTS']

	@staticmethod
	def fermi_gbm_fin_pos(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN classic text Fermi GBM Final Position Notices (gcn.classic.text.FERMI_GBM_FIN_POS) from the Fermi Gamma-ray Space Telescope via GCN Kafka.

		GBM Final Position notices occur when humans are involved in the analysis. The humans can make better selections on the data used. There can be 0, 1, or possibly more instances of this notice type per trigger. They will be issued for the 10-20% brightest of the GRBs. There may be others, less bright, that also get a final notice, but they will have a larger location uncertainty.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		params = Alerts.decode_classic_notice(value)

		rec_title = params['TITLE']
		rec_notice_date = params['NOTICE_DATE']
		rec_notice_type = params['NOTICE_TYPE']
		rec_record_num = params['RECORD_NUM']
		rec_trigger_num = params['TRIGGER_NUM']
		rec_grb_ra = params['GRB_RA']
		rec_grb_dec = params['GRB_DEC']
		rec_grb_error = params['GRB_ERROR']
		rec_grb_date = params['GRB_DATE']
		rec_grb_time = params['GRB_TIME']
		rec_grb_phi = params['GRB_PHI']
		rec_grb_theta = params['GRB_THETA']
		rec_e_range = params['E_RANGE']
		rec_loc_algorithm = params['LOC_ALGORITHM']
		rec_lc_url = params['LC_URL']
		rec_loc_url = params['LOC_URL']
		rec_sun_postn = params['SUN_POSTN']
		rec_sun_dist = params['SUN_DIST']
		rec_moon_postn = params['MOON_POSTN']
		rec_moon_dist = params['MOON_DIST']
		rec_moon_illum = params['MOON_ILLUM']
		rec_gal_coords = params['GAL_COORDS']
		rec_ecl_coords = params['ECL_COORDS']
		rec_comments = params['COMMENTS']

	@staticmethod
	def fermi_gbm_flt_pos(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN classic text Fermi GBM Flight Position Notices (gcn.classic.text.FERMI_GBM_FLT_POS) from the Fermi Gamma-ray Space Telescope via GCN Kafka.

		Contains the RA/Dec location for the burst detected by the GBM. The positions are calculated by the onboard flight software. They are issued 1-5 times per burst. The flight position notice comes second in the sequence of notices on the triggers/bursts. Since the position is based on the least-possible amount of data in the processing (just a small initial portion of the burst’s lightcurve), it has the lowest significance in the localization process. The uncertainty in the position is ~20 deg (1-sigma, radius) for the at-threshold bursts and ~10 deg for the bright bursts. (Both the ‘at-threshold’ and ‘bright’ refer to only the amount of photons in the time interval of the trigger-sampling interval of the burst lightcurve, not the total burst duration.) This notice contains detections of both GRBs and hard X-ray transients. There is both flight- and ground-software in place to correctly identify bursts from transients, but this identification is not 100% perfect.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		params = Alerts.decode_classic_notice(value)

		rec_title = params['TITLE']
		rec_notice_date = params['NOTICE_DATE']
		rec_notice_type = params['NOTICE_TYPE']
		rec_record_num = params['RECORD_NUM']
		rec_trigger_num = params['TRIGGER_NUM']
		rec_grb_ra = params['GRB_RA']
		rec_grb_dec = params['GRB_DEC']
		rec_grb_error = params['GRB_ERROR']
		rec_grb_inten = params['GRB_INTEN']
		rec_data_signif = params['DATA_SIGNIF']
		rec_integ_time = params['INTEG_TIME']
		rec_grb_date = params['GRB_DATE']
		rec_grb_time = params['GRB_TIME']
		rec_grb_phi = params['GRB_PHI']
		rec_grb_theta = params['GRB_THETA']
		rec_hard_ratio = params['HARD_RATIO']
		rec_loc_algorithm = params['LOC_ALGORITHM']
		rec_most_likely = params['MOST_LIKELY']
		rec_2nd_most_likely = params['2nd_MOST_LIKELY']
		rec_detectors = params['DETECTORS']
		rec_lc_url = params['LC_URL']
		rec_sun_postn = params['SUN_POSTN']
		rec_sun_dist = params['SUN_DIST']
		rec_moon_postn = params['MOON_POSTN']
		rec_moon_dist = params['MOON_DIST']
		rec_moon_illum = params['MOON_ILLUM']
		rec_gal_coords = params['GAL_COORDS']
		rec_ecl_coords = params['ECL_COORDS']
		rec_comments = params['COMMENTS']

	@staticmethod
	def fermi_gbm_gnd_pos(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN classic text Fermi GBM Ground Position Notices (gcn.classic.text.FERMI_GBM_GND_POS) from the Fermi Gamma-ray Space Telescope via GCN Kafka.

		Contains the RA/Dec location for the burst detected by the GBM. The positions are calculated by automated ground software as soon as they are received on the ground. More sophisticated algorithms can be applied to the data to improve the location accuracy. More data from the ongoing progression of the burst lightcurve are used in this calculation. There can be 0, 1, or more instances of this notice type per trigger.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		rec_title = params['TITLE']
		rec_notice_date = params['NOTICE_DATE']
		rec_notice_type = params['NOTICE_TYPE']
		rec_record_num = params['RECORD_NUM']
		rec_trigger_num = params['TRIGGER_NUM']
		rec_grb_ra = params['GRB_RA']
		rec_grb_dec = params['GRB_DEC']
		rec_grb_error = params['GRB_ERROR']
		rec_data_signif = params['DATA_SIGNIF']
		rec_data_interval = params['DATA_INTERVAL']
		rec_grb_date = params['GRB_DATE']
		rec_grb_time = params['GRB_TIME']
		rec_grb_phi = params['GRB_PHI']
		rec_grb_theta = params['GRB_THETA'] # ?
		rec_grb_zenith = params['GRB_ZENITH'] # ?
		rec_e_range = params['E_RANGE']
		rec_loc_algorithm = params['LOC_ALGORITHM']
		rec_lc_url = params['LC_URL']
		rec_sun_postn = params['SUN_POSTN']
		rec_sun_dist = params['SUN_DIST']
		rec_moon_postn = params['MOON_POSTN']
		rec_moon_dist = params['MOON_DIST']
		rec_moon_illum = params['MOON_ILLUM']
		rec_gal_coords = params['GAL_COORDS']
		rec_ecl_coords = params['ECL_COORDS']
		rec_comments = params['COMMENTS']

	@staticmethod
	def fermi_gbm_pos_test(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN classic text Fermi GBM Position Test Notices (gcn.classic.text.FERMI_GBM_POS_TEST) from the Fermi Gamma-ray Space Telescope via GCN Kafka.

		Identical in format and content to the GBM Flight Position Notice, except that it contains a computer generated RA/Dec location and all the other fields. This “test” notice is generated by the GCN/TAN computer every ~3.6 hrs, allowing the receiving site to “practice” on the GBM Flight Position Notice.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		params = Alerts.decode_classic_notice(value)

		rec_title = params['TITLE']
		rec_notice_date = params['NOTICE_DATE']
		rec_notice_type = params['NOTICE_TYPE']
		rec_record_num = params['RECORD_NUM']
		rec_trigger_num = params['TRIGGER_NUM']
		rec_grb_ra = params['GRB_RA']
		rec_grb_dec = params['GRB_DEC']
		rec_grb_error = params['GRB_ERROR']
		rec_grb_inten = params['GRB_INTEN']
		rec_data_signif = params['DATA_SIGNIF']
		rec_integ_time = params['INTEG_TIME']
		rec_grb_date = params['GRB_DATE']
		rec_grb_time = params['GRB_TIME']
		rec_grb_phi = params['GRB_PHI']
		rec_grb_theta = params['GRB_THETA']
		rec_hard_ratio = params['HARD_RATIO']
		rec_loc_algorithm = params['LOC_ALGORITHM']
		rec_most_likely = params['MOST_LIKELY']
		rec_2nd_most_likely = params['2nd_MOST_LIKELY']
		rec_detectors = params['DETECTORS']
		rec_lc_url = params['LC_URL'] # KeyError
		rec_sun_postn = params['SUN_POSTN']
		rec_sun_dist = params['SUN_DIST']
		rec_moon_postn = params['MOON_POSTN']
		rec_moon_dist = params['MOON_DIST']
		rec_moon_illum = params['MOON_ILLUM']
		rec_gal_coords = params['GAL_COORDS']
		rec_ecl_coords = params['ECL_COORDS']
		rec_comments = params['COMMENTS']

	@staticmethod
	def fermi_gbm_subthresh(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN classic text Fermi GBM Subthreshold Notices (gcn.classic.text.FERMI_GBM_SUBTHRESH) from the Fermi Gamma-ray Space Telescope via GCN Kafka.

		This is a separate stream of transients. These are produced from a ground pipeline processing the data looking for transients that were below the onboard s/w trigger threshold level. These "subthreshold" events are typically used for coincidence searches with other data streams.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		params = Alerts.decode_classic_notice(value)

		rec_title = params['TITLE']
		rec_notice_date = params['NOTICE_DATE']
		rec_notice_type = params['NOTICE_TYPE']
		trans_num = params['TRANS_NUM']
		rec_full_id_num = params['FULL_ID_NUM']
		rec_trans_ra = params['TRANS_RA']
		rec_trans_dec = params['TRANS_DEC']
		rec_trans_error = params['TRANS_ERROR']
		rec_trans_duration = params['TRANS_DURATION']
		rec_trans_date = params['TRANS_DATE']
		rec_trans_time = params['TRANS_TIME']
		rec_earth_angle = params['EARTH_ANGLE']
		rec_spectral_class = params['SPECTRAL_CLASS']
		rec_type_class = params['TYPE_CLASS']
		rec_reliability = params['RELIABILITY']
		rec_healpix_url = params['HEALPIX_URL']
		rec_map_url = params['MAP_URL']
		rec_lc_url = params['LC_URL']
		rec_sun_postn = params['SUN_POSTN']
		rec_sun_dist = params['SUN_DIST']
		rec_moon_postn = params['MOON_POSTN']
		rec_moon_dist = params['MOON_DIST']
		rec_moon_illum = params['MOON_ILLUM']
		rec_gal_coords = params['GAL_COORDS']
		rec_ecl_coords = params['ECL_COORDS']
		rec_comments = params['COMMENTS']

	@staticmethod
	def fermi_lat_monitor(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN classic text Fermi LAT Monitor Notices (gcn.classic.text.FERMI_LAT_MONITOR) from the Fermi Gamma-ray Space Telescope via GCN Kafka.

		Reports the detection of flares from previously known sources. The position comes from the LAT ground-processing software plus possible human inspection and refinement of the results. LAT monitors 168 sources during its all-sky scanning mode. When one of these sources has a flare with a flux increase typically >5 sigma, a LAT Monitor Notice is issued. No location error is given because these are flares from known sources.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		params = Alerts.decode_classic_notice(value)

		rec_title = params['TITLE']
		rec_notice_date = params['NOTICE_DATE']
		rec_notice_type = params['NOTICE_TYPE']
		rec_source = params['SOURCE']
		rec_ref_num = params['REF_NUM']
		rec_ra = params['RA']
		rec_dec = params['DEC']
		rec_curr_flux = params['CURR_FLUX']
		rec_base_flux = params['BASE_FLUX']
		rec_time_scale = params['TIME_SCALE']
		rec_outburst_date = params['OUTBURST_DATE']
		rec_outburst_time = params['OUTBURST_TIME']
		rec_soln_status = params['SOLN_STATUS']
		rec_lc_url = params['LC_URL']
		rec_sun_postn = params['SUN_POSTN']
		rec_sun_dist = params['SUN_DIST']
		rec_moon_postn = params['MOON_POSTN']
		rec_moon_dist = params['MOON_DIST']
		rec_moon_illum = params['MOON_ILLUM']
		rec_gal_coords = params['GAL_COORDS']
		rec_ecl_coords = params['ECL_COORDS']
		rec_comments = params['COMMENTS']

	@staticmethod
	def fermi_lat_offline(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN classic text Fermi LAT Position Offline Notices (gcn.classic.text.FERMI_LAT_OFFLINE) from the Fermi Gamma-ray Space Telescope via GCN Kafka.

		Contains an RA/Dec location for a burst detected by LAT. The position comes from the LAT ground processing software plus human inspection and refinement of the results. The position uncertainty is 10-60 arcmin (depending on the number of photons detected and their energies). They are issued only once (rarely twice) per ground-detected burst.
		
		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		params = Alerts.decode_classic_notice(value)

		rec_title = params['TITLE']
		rec_notice_date = params['NOTICE_DATE']
		rec_notice_type = params['NOTICE_TYPE']
		rec_trigger_num = params['TRIGGER_NUM']
		rec_grb_ra = params['GRB_RA']
		rec_grb_dec = params['GRB_DEC']
		rec_grb_error = params['GRB_ERROR']
		rec_grb_date = params['GRB_DATE']
		rec_grb_time = params['GRB_TIME']
		rec_trigger_id = params['TRIGGER_ID']
		rec_sun_postn = params['SUN_POSTN']
		rec_sun_dist = params['SUN_DIST']
		rec_moon_postn = params['MOON_POSTN']
		rec_moon_dist = params['MOON_DIST']
		rec_moon_illum = params['MOON_ILLUM']
		rec_gal_coords = params['GAL_COORDS']
		rec_ecl_coords = params['ECL_COORDS']
		rec_comments = params['COMMENTS']

	@staticmethod
	def fermi_lat_pos_diag(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN classic text Fermi LAT Position Diagnostic Notices (gcn.classic.text.FERMI_LAT_POS_DIAG) from the Fermi Gamma-ray Space Telescope via GCN Kafka.

		Contains an RA/Dec location for the burst/afterglow detected by LAT. They are issued only once per burst. It has the longest integration of source and background photons. The position comes from the LAT flight software. The position uncertainty is 60-120 arcmin (depending on the number of photons detected and their energies).

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		params = Alerts.decode_classic_notice(value)

		rec_title = params['TITLE']
		rec_notice_date = params['NOTICE_DATE']
		rec_notice_type = params['NOTICE_TYPE']
		rec_record_num = params['RECORD_NUM']
		rec_trigger_num = params['TRIGGER_NUM']
		rec_grb_ra = params['GRB_RA']
		rec_grb_dec = params['GRB_DEC']
		rec_grb_error = params['GRB_ERROR']
		rec_grb_inten_tot = params['GRB_INTEN_TOT']
		rec_grb_inten1 = params['GRB_INTEN1']
		rec_grb_inten2 = params['GRB_INTEN2']
		rec_grb_inten3 = params['GRB_INTEN3']
		rec_grb_inten4 = params['GRB_INTEN4']
		rec_integ_dur = params['GRB_INTEG_DUR']
		rec_trigger_index = params['TRIGGER_INDEX']
		rec_grb_date = params['GRB_DATE']
		rec_grb_time = params['GRB_TIME']
		rec_soln_status = params['SOLN_STATUS']
		rec_temp_test_stat = params['TEMP_TEST_STAT']
		rec_image_test_stat = params['IMAGE_TEST_STAT']
		rec_sun_postn = params['SUN_POSTN']
		rec_sun_dist = params['SUN_DIST']
		rec_moon_postn = params['MOON_POSTN']
		rec_moon_dist = params['MOON_DIST']
		rec_moon_illum = params['MOON_ILLUM']
		rec_gal_coords = params['GAL_COORDS']
		rec_ecl_coords = params['ECL_COORDS']
		rec_comments = params['COMMENTS']

	@staticmethod
	def fermi_lat_pos_ini(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN classic text Fermi LAT Position Initial Notices (gcn.classic.text.FERMI_LAT_POS_INI) from the Fermi Gamma-ray Space Telescope via GCN Kafka.

		Contains an RA/Dec location for the burst/afterglow detected by LAT. It starts a sequence of LAT notices. It is issued only once per burst. It has the shortest integration of source and background photons. The position comes from the LAT flight software scanning the time and spatial domains. The position uncertainty is 60-120 arcmin (depending on the number of photons detected and their energies).

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		params = Alerts.decode_classic_notice(value)

		rec_title = params['TITLE']
		rec_notice_date = params['NOTICE_DATE']
		rec_notice_type = params['NOTICE_TYPE']
		rec_record_num = params['RECORD_NUM']
		rec_trigger_num = params['TRIGGER_NUM']
		rec_grb_ra = params['GRB_RA']
		rec_grb_dec = params['GRB_DEC']
		rec_grb_error = params['GRB_ERROR']
		rec_grb_inten_tot = params['GRB_INTEN_TOT']
		rec_grb_inten1 = params['GRB_INTEN1']
		rec_grb_inten2 = params['GRB_INTEN2']
		rec_grb_inten3 = params['GRB_INTEN3']
		rec_grb_inten4 = params['GRB_INTEN4']
		rec_integ_dur = params['INTEG_DUR']
		rec_trigger_index = params['TRIGGER_INDEX']
		rec_grb_date = params['GRB_DATE']
		rec_grb_time = params['GRB_TIME']
		rec_soln_status = params['SOLN_STATUS']
		rec_temp_test_stat = params['TEMP_TEST_STAT']
		rec_image_test_stat = params['IMAGE_TEST_STAT']
		rec_sun_postn = params['SUN_POSTN']
		rec_sun_dist = params['SUN_DIST']
		rec_moon_postn = params['MOON_POSTN']
		rec_moon_dist = params['MOON_DIST']
		rec_moon_illum = params['MOON_ILLUM']
		rec_gal_coords = params['GAL_COORDS']
		rec_ecl_coords = params['ECL_COORDS']
		rec_comments = params['COMMENTS']

	@staticmethod
	def fermi_lat_pos_test(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN classic text Fermi LAT Position Test Notices (gcn.classic.text.FERMI_LAT_POS_TEST) from the Fermi Gamma-ray Space Telescope via GCN Kafka.

		Identical in format and content to the LAT Position Update Notice, except that it contains a computer generated RA/Dec location and all the other fields. This “test” notice is generated by the GCN/TAN computer every ~3.6 hrs, allowing the receiving site to “practice” on the LAT Position Update Notice.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		params = Alerts.decode_classic_notice(value)

		rec_title = params['TITLE']
		rec_notice_date = params['NOTICE_DATE']
		rec_notice_type = params['NOTICE_TYPE']
		rec_record_num = params['RECORD_NUM']
		rec_trigger_num = params['TRIGGER_NUM']
		rec_grb_ra = params['GRB_RA']
		rec_grb_dec = params['GRB_DEC']
		rec_grb_error = params['GRB_ERROR']
		rec_grb_inten_tot = params['GRB_INTEN_TOT']
		rec_grb_inten1 = params['GRB_INTEN1']
		rec_grb_inten2 = params['GRB_INTEN2']
		rec_grb_inten3 = params['GRB_INTEN3']
		rec_grb_inten4 = params['GRB_INTEN4']
		rec_integ_dur = params['INTEG_DUR']
		rec_trigger_index = params['TRIGGER_INDEX'] #KeyError
		rec_grb_date = params['GRB_DATE']
		rec_grb_time = params['GRB_TIME']
		rec_soln_status = params['SOLN_STATUS']
		rec_temp_test_stat = params['TEMP_TEST_STAT']
		rec_image_test_stat = params['IMAGE_TEST_STAT']
		rec_sun_postn = params['SUN_POSTN']
		rec_sun_dist = params['SUN_DIST']
		rec_moon_postn = params['MOON_POSTN']
		rec_moon_dist = params['MOON_DIST']
		rec_moon_illum = params['MOON_ILLUM']
		rec_gal_coords = params['GAL_COORDS']
		rec_ecl_coords = params['ECL_COORDS']
		rec_comments = params['COMMENTS']

	@staticmethod
	def fermi_lat_pos_upd(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN classic text Fermi LAT Position Update Notices (gcn.classic.text.FERMI_LAT_POS_UPD) from the Fermi Gamma-ray Space Telescope via GCN Kafka.

		Contains an RA/Dec location for the burst/afterglow detected by LAT. They are issued six times per burst, but only the second issuance is available to the public. The choice of the second issuance is a trade between accumulating source photons without accumulating too many background photons. The position in this notice type comes from the LAT flight software. The position uncertainty is 60-120 arcmin (depending on the number of photons detected and their energies).

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		params = Alerts.decode_classic_notice(value)

		rec_title = params['TITLE']
		rec_notice_date = params['NOTICE_DATE']
		rec_notice_type = params['NOTICE_TYPE']
		rec_record_num = params['RECORD_NUM']
		rec_trigger_num = params['TRIGGER_NUM']
		rec_grb_ra = params['GRB_RA']
		rec_grb_dec = params['GRB_DEC']
		rec_grb_error = params['GRB_ERROR']
		rec_grb_inten_tot = params['GRB_INTEN_TOT']
		rec_grb_inten1 = params['GRB_INTEN1']
		rec_grb_inten2 = params['GRB_INTEN2']
		rec_grb_inten3 = params['GRB_INTEN3']
		rec_grb_inten4 = params['GRB_INTEN4']
		rec_integ_dur = params['INTEG_DUR']
		rec_trigger_index = params['TRIGGER_INDEX']
		rec_grb_date = params['GRB_DATE']
		rec_grb_time = params['GRB_TIME']
		rec_soln_status = params['SOLN_STATUS']
		rec_temp_test_stat = params['TEMP_TEST_STAT']
		rec_image_test_stat = params['IMAGE_TEST_STAT']
		rec_sun_postn = params['SUN_POSTN']
		rec_sun_dist = params['SUN_DIST']
		rec_moon_postn = params['MOON_POSTN']
		rec_moon_dist = params['MOON_DIST']
		rec_moon_illum = params['MOON_ILLUM']
		rec_gal_coords = params['GAL_COORDS']
		rec_ecl_coords = params['ECL_COORDS']
		rec_comments = params['COMMENTS']

	@staticmethod
	def fermi_pointdir(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN classic text Fermi Pointing Direction Notices (gcn.classic.text.FERMI_POINTDIR) from the Fermi Gamma-ray Space Telescope via GCN Kafka.

		Gives the current point direction of the spacecraft and the pointing direction at 2-minute intervals for the next 60 minutes. This extrapolation into the future is based on the s/c ephemeris published by the Fermi MOC. This notice type is issued every 60 minutes. The Fermi spacecraft is zenith-pointing oriented, but to increase the sky coverage, the spacecraft nods (slews) +/- ~30 deg north-n-south once per orbit. These ‘nods’ show up as discontinuities in the RA/Dec samples every 96 minutes. On rare occasions the actual spacecraft pointing direction will deviate from the planned pointing direction because of a TOO or when the spacecraft repoints in response to a GBM burst detection.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		params = Alerts.decode_classic_notice(value)

		rec_title = params['TITLE']
		rec_notice_date = params['NOTICE_DATE']
		rec_notice_type = params['NOTICE_TYPE']
		rec_curr_point_ra = params['CURR_POINT_RA']
		rec_curr_point_dec = params['CURR_POINT_DEC']
		rec_curr_date = params['CURR_DATE']
		rec_curr_time = params['CURR_TIME']
		rec_delta_time = params['DELTA_TIME']
		rec_sun_postn = params['SUN_POSTN']
		rec_sun_dist = params['SUN_DIST']
		rec_moon_postn = params['MOON_POSTN']
		rec_moon_dist = params['MOON_DIST']
		rec_moon_illum = params['MOON_ILLUM']
		rec_gal_coords = params['GAL_COORDS']
		rec_ecl_coords = params['ECL_COORDS']
		rec_future_ra_dec = params['FUTURE_RA_DEC']
		rec_comments = params['COMMENTS']

	@staticmethod
	def filter_notice(value, topic):

		print('\nGCN Notice [' + str(topic) + ']')

		if topic == 'gcn.circulars':
			Alerts.gcn_circulars(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_GBM_ALERT':
			Alerts.fermi_gbm_alert(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_GBM_FIN_POS':
			Alerts.fermi_gbm_fin_pos(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_GBM_FLT_POS':
			Alerts.fermi_gbm_flt_pos(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_GBM_GND_POS':
			Alerts.fermi_gbm_gnd_pos(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_GBM_POS_TEST':
			Alerts.fermi_gbm_pos_test(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_GBM_SUBTHRESH':
			Alerts.fermi_gbm_subthresh(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_LAT_MONITOR':
			Alerts.fermi_lat_monitor(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_LAT_POS_DIAG':
			Alerts.fermi_lat_pos_diag(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_LAT_POS_INI':
			Alerts.fermi_lat_pos_ini(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_LAT_OFFLINE':
			Alerts.fermi_lat_pos_offline(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_LAT_POS_TEST':
			Alerts.fermi_lat_pos_test(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_LAT_POS_UPD':
			Alerts.fermi_lat_pos_upd(value, topic, False)

		elif topic == 'gcn.classic.text.FERMI_POINTDIR':
			Alerts.fermi_pointdir(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_ACTUAL_POINTDIR':
			Alerts.swift_actual_pointdir(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_BAT_GRB_LC':
			Alerts.swift_bat_grb_lc(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_BAT_GRB_POS_ACK':
			Alerts.swift_bat_grb_pos_ack(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_BAT_GRB_POS_TEST':
			Alerts.swift_bat_grb_pos_test(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_BAT_QL_POS':
			Alerts.swift_bat_grb_ql_pos(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_BAT_SCALEDMAP':
			Alerts.swift_bat_scaledmap(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_BAT_TRANS':
			Alerts.swift_bat_trans(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_FOM_OBS':
			Alerts.swift_fom_obs(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_POINTDIR':
			Alerts.swift_pointdir(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_SC_SLEW':
			Alerts.swift_sc_slew(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_TOO_FOM':
			Alerts.swift_too_fom(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_TOO_SC_SLEW':
			Alerts.swift_too_sc_slew(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_DBURST':
			Alerts.swift_uvot_dburst(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_DBURST_PROC':
			Alerts.swift_uvot_dburst_proc(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_EMERGENCY':
			Alerts.swift_uvot_emergency(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_FCHART':
			Alerts.swift_uvot_fchart(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_FCHART_PROC':
			Alerts.swift_uvot_fchart_proc(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_POS':
			Alerts.swift_uvot_pos(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_POS_NACK':
			Alerts.swift_uvot_pos_nack(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_CENTROID':
			Alerts.swift_xrt_centroid(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_EMERGENCY':
			Alerts.swift_xrt_emergency(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_IMAGE':
			Alerts.swift_xrt_image(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_IMAGE_PROC':
			Alerts.swift_xrt_image_proc(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_LC':
			Alerts.swift_xrt_lc(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_POSITION':
			Alerts.swift_xrt_position(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_SPECTRUM':
			Alerts.swift_xrt_spectrum(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_SPECTRUM_PROC':
			Alerts.swift_xrt_spectrum_proc(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_SPER':
			Alerts.swift_xrt_sper(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_SPER_PROC':
			Alerts.swift_xrt_sper_proc(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_THRESHPIX':
			Alerts.swift_xrt_threshpix(value, topic, False)

		elif topic == 'gcn.classic.text.SWIFT_XRT_THRESHPIX_PROC':
			Alerts.swift_xrt_threshpix_proc(value, topic, False)

		elif topic == 'gcn.notices.chime.frb':
			Alerts.chime_frb(value, topic, False)

		elif topic == 'gcn.notices.dsa110.frb':
			Alerts.dsa110_frb(value, topic, False)

		elif topic == 'gcn.notices.einstein_probe.wxt.alert':
			Alerts.einstein_probe_wxt_alert(value, topic, False)

		elif topic == 'gcn.notices.icecube.gold_bronze_track_alerts':
			Alerts.icecube_gold_bronze_track_alerts(value, topic, False)

		elif topic == 'gcn.notices.icecube.lvk_nu_track_search':
			Alerts.icecube_lvk_nu_track_search(value, topic, False)

		elif topic == 'gcn.notices.maxi.gsc.known':
			Alerts.maxi_gsc_known(value, topic, False)

		elif topic == 'gcn.notices.maxi.gsc.unknown':
			Alerts.maxi_gsc_unknown(value, topic, False)

		elif topic == 'gcn.notices.superk.sn_alert':
			Alerts.superk_sn_alert(value, topic, False)

		elif topic == 'gcn.notices.svom.voevent.eclairs':
			Alerts.svom_eclairs(value, topic, False)

		elif topic == 'gcn.notices.swift.bat.guano':
			Alerts.swift_bat_guano(value, topic, False)

		elif topic == 'igwn.gwalert':
			Alerts.igwn_gwalert(value, topic, False)

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
	def get_notice_topics():
		'''This function reads a file listing available GCN topics and returns them as a list.

		:return topic_list - a list of available GCN topics
		:return n_topics - the number of available GCN topics
		'''

		topic_path = Configuration.PLOUTON_DIRECTORY + 'alerts/gcn_topics' + Configuration.TABLE_EXTENSION

		topic_file = open(topic_path, 'r')
		topic_list = []

		for ln in topic_file:
			if ln[0] == '#':
				pass
			else:
				ln = ln.rstrip('\n')
				topic_list.append(ln)

		n_topics = len(topic_list)

		return topic_list, n_topics

	@staticmethod
	def icecube_gold_bronze_track_alerts(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN Notices (high-energy single-track events) from the IceCube Neutrino Observatory via GCN Kafka.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		record = json.loads(value)

		# $schema
		try:
			record_schema = record['$schema']
		except KeyError:
			record_schema = None

		# mission
		try:
			record_mission = record['mission']
		except KeyError:
			record_mission = None

		# instrument
		try:
			record_instrument = record['instrument']
		except KeyError:
			record_instrument = None

		# messenger
		try:
			record_messenger = record['messenger']
		except KeyError:
			record_messenger = None

		# pipeline
		try:
			record_pipeline = record['pipeline']
		except KeyError:
			record_pipeline = None

		# record_number
		try:
			record_record_number = record['record_number']
		except KeyError:
			record_record_number = None

		# event_name
		try:
			record_event_name = record['event_name']
		except KeyError:
			record_event_name = None

		# id
		try:
			record_id = record['id']
		except KeyError:
			record_id = None

		# alert_datetime
		try:
			record_alert_datetime = record['alert_datetime']
		except KeyError:
			record_alert_datetime = None

		# alert_type
		try:
			record_alert_type = record['alert_type']
		except KeyError:
			record_alert_type = None

		# alert_tense
		try:
			record_alert_tense = record['alert_tense']
		except KeyError:
			record_alert_tense = None

		# alert_topology
		try:
			record_alert_topology = record['alert_topology']
		except KeyError:
			record_alert_topology = None

		# number_of_events
		try:
			record_number_of_events = record['number_of_events']
		except KeyError:
			record_number_of_events = None

		# ra
		try:
			record_ra = record['ra']
		except KeyError:
			record_ra = None

		# dec
		try:
			record_dec = record['dec']
		except KeyError:
			record_dec = None

		# ra_dec_error
		try:
			record_ra_dec_error = record['ra_dec_error']
		except KeyError:
			record_ra_dec_error = None

		# containment_probability
		try:
			record_containment_probability = record['containment_probability']
		except KeyError:
			record_containment_probability = None

		# systematic_included
		try:
			record_systematic_included = record['systematic_included']
		except KeyError:
			record_systematic_included = None

		# healpix_url
		try:
			record_healpix_url = record['healpix_url']
		except KeyError:
			record_healpix_url = None

		# trigger_time
		try:
			record_trigger_time = record['trigger_time']
		except KeyError:
			record_trigger_time = None

		# nu_energy
		try:
			record_nu_energy = record['nu_energy']
		except KeyError:
			record_nu_energy = None

		# p_astro
		try:
			record_p_astro = record['p_astro']
		except KeyError:
			record_p_astro = None

		# far
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
	def icecube_lvk_nu_track_search(value, topic, alert=Configuration.ALERT):
		'''This function processes GCN Notices (coincident GW track events) from the IceCube Neutrino Observatory via GCN Kafka.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		record = json.loads(value)

		# $schema [string] - GCN schema structure according to the nasa-gcn/gcn-schema GitHub project
		try:
			record_schema = record['$schema']
		except KeyError:
			record_schema = None

		# type [string] - alert type (i.e. 'IceCube LVK Alert Nu Track Search')
		try:
			record_type = record['type']
		except KeyError:
			record_type = None

		# reference [string] - specifies the LVK alert GCN being followed up including the GW map revision used in the analysis
		try:
			record_reference = record['reference']
		except KeyError:
			record_reference = None

		# ref_ID [string] - the LVK event name being followed up
		try:
			record_ref_id = record['ref_ID']
		except KeyError:
			record_ref_id = None

		# alert_datetime [string] - timestamp of the GCN notice (UTC time)
		try:
			record_alert_datetime = record['alert_datetime']
		except KeyError:
			record_alert_datetime = None

		# trigger_time [string] - timestamp of the merger event (UTC time)
		try:
			record_trigger_time = record['trigger_time']
		except KeyError:
			record_trigger_time = None

		# observation_start [string] - beginning of the analysis time window (UTC time)
		try:
			record_observation_start = record['observation_start']
		except KeyError:
			record_observation_start = None

		# observation_stop [string] - end of the analysis time window (UTC time)
		try:
			record_observation_stop = record['observation_stop']
		except KeyError:
			record_observation_stop = None

		# observation_livetime [number] - duration of the analysis [s]
		try:
			record_observation_livetime = record['observation_livetime']
		except KeyError:
			record_observation_livetime = None

		# pval_generic [number] - the overall search results p-values for the generic search
		try:
			record_pval_generic = record['pval_generic']
		except KeyError:
			record_pval_generic = None

		# pval_bayesian [number] - the overall search results p-values for the Bayesian search
		try:
			record_pval_bayesian = record['pval_bayesian']
		except KeyError:
			record_pval_bayesian = None

		# n_events_coincident [number] - the number of events in the time window with a per-event p-value less than 10%
		try:
			record_n_events_coincident = record['n_events_coicident']
			if record_n_events_coincident > 0:
				for i in range(record_n_events_coincident):
					# event_dt [number] - relative timing of neutrino to GW alert [s]
					try:
						record_event_dt = record['coincident_events'][i]['event_dt']
					except KeyError:
						record_event_dt = None

					# ra [number] - reconstructed neutrino event right ascension (J2000) [deg]
					# dec [number] - reconstructed neutrino event declination (J2000) [deg]
					# ra_dec_error [number] - circular detection error at 90% confidence [deg]
					try:
						record_localization = record['coincident_events'][i]['localization']
						record_localization_ra = record_localization['ra']
						record_localization_dec = record_localization['dec']
						record_ra_dec_error = record_localization['ra_dec_error']
						record_containment_probability = record_localization['containment_probability']
						record_systematic_included = record_localization['systematic_included']
					except KeyError:
						record_localization = None
						record_localization_ra = None
						record_localization_dec = None
						record_ra_dec_error = None
						record_containment_probability = None
						record_systematic_included = None

					# id [string] - unique identifier label for the coincident neutrino event, formatted as 'RUNID_EVENTID'
					try:
						record_id = record['coincident_events'][i]['id']
					except KeyError:
						record_id = None

					# event_pval_generic [number] - per-event p-value for the generic search
					try:
						record_event_pval_generic = record['coincident_events'][i]['event_pval_generic']
					except KeyError:
						record_event_pval_generic = None

					# event_pval_bayesian [number] - per-event p-value for the Bayesian search
					try:
						record_event_pval_bayesian = record['coincident_events'][i]['event_pval_bayesian']
					except KeyError:
						record_event_pval_bayesian = None

		except KeyError:
			record_n_events_coincident = None

		# most_probable_direction [number] - the most likely direction from the generic neutrino source coincidence search (J2000) [deg]
		try:
			record_most_probable_direction = record['most_probable_direction']
			record_most_probable_direction_ra = record['most_probable_direction']['ra']
			record_most_probable_direction_dec = record['most_probable_direction']['dec']
		except KeyError:
			record_most_probable_direction = None
			record_most_probable_direction_ra = None
			record_most_probable_direction_dec = None

		# flux_sensitivity [number] - time integrated flux sensitivity range, in the format (min, max) sensitivity [GeV cm^-2]
		# sensitive_energy_range [number] - energy sensitivity range, in the format (lower, upper) energy [GeV]
		try:
			record_neutrino_flux_sensitivity_range = record['neutrino_flux_sensitivity_range']
			record_flux_sensitivity = record['neutrino_flux_sensitivity_range']['flux_sensitivity']
			record_sensitive_energy_range = record['neutrino_flux_sensitivity_range']['sensitive_energy_range']
		except KeyError:
			record_neutrino_flux_sensitivity_range = None
			record_flux_sensitivity = None
			record_sensitive_energy_range = None

		# print the alert contents
		print('\tSchema:', record_schema)
		print('\tType:', record_type)
		print('\tReference:', record_reference)
		print('\tRef ID:', record_ref_id)
		print('\tAlert datetime:', record_alert_datetime)
		print('\tTrigger time:', record_trigger_time)
		print('\tObservation start:', record_observation_start)
		print('\tObservation stop:', record_observation_stop)
		print('\tObservation livetime:', record_observation_livetime, 's')
		print('\tp-value (generic):', record_pval_generic)
		print('\tp-value (Bayesian):', record_pval_bayesian)
		print('\tNumber of coincident events:', record_n_events_coincident)

	@staticmethod
	def igwn_gwalert(value, topic, alert=Configuration.ALERT):

		record = json.loads(value)

	@staticmethod
	def maxi_gsc_known(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def maxi_gsc_unknown(value, topic, alert=Configuration.ALERT):
		return

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
	def superk_sn_alert(value, topic, alert=Configuration.ALERT):

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
	def svom_eclairs(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def svom_grm(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def svom_mxt(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_actual_pointdir(value, topic, alert=Configuration.ALERT):
		'''This function processes Swift Actual Pointing Direction Notices from the Neil Gehrels Swift Observatory via GCN Classic over Kafka.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		record = Alerts.decode_classic_notice(value)

		rec_title = record['TITLE']
		str_title = '\tTitle: ' + str(rec_title)
		print(str_title)

		rec_notice_date = record['NOTICE_DATE']
		str_notice_date = '\tNotice date: ' + str(rec_notice_date)
		print(str_notice_date)

		rec_notice_type = record['NOTICE_TYPE']
		str_notice_type = '\tNotice type: ' + str(rec_notice_type)
		print(str_notice_type)

		rec_curr_point_ra = record['CURR_POINT_RA']
		str_curr_point_ra = '\tCurrent point (RA): ' + str(rec_curr_point_ra)
		print(str_curr_point_ra)

		rec_curr_point_dec = record['CURR_POINT_DEC']
		str_curr_point_dec = '\tCurrent point (Dec): ' + str(rec_curr_point_dec)
		print(str_curr_point_dec)

		rec_curr_point_roll = record['CURR_POINT_ROLL']
		str_curr_point_roll = '\tCurrent point (Roll): ' + str(rec_curr_point_roll)
		print(str_curr_point_roll)

		rec_slew_time = record['SLEW_TIME']
		str_slew_time = '\tSlew time: ' + str(rec_slew_time)
		print(str_slew_time)

		rec_slew_date = record['SLEW_DATE']
		str_slew_date = '\tSlew date: ' + str(rec_slew_date)
		print(str_slew_date)

		rec_tgt_num = record['TGT_NUM']
		str_tgt_num = '\tTarget number: ' + str(rec_tgt_num)
		print(str_tgt_num)

		rec_sun_postn = record['SUN_POSTN']
		str_sun_postn = '\tSun position: ' + str(rec_sun_postn)
		print(str_sun_postn)
		
		rec_sun_dist = record['SUN_DIST']
		str_sun_dist = '\tSun distance: ' + str(rec_sun_dist)
		print(str_sun_dist)
		
		rec_moon_postn = record['MOON_POSTN']
		str_moon_postn = '\tMoon position: ' + str(rec_moon_postn)
		print(str_moon_postn)
		
		rec_moon_dist = record['MOON_DIST']
		str_moon_dist = '\tMoon distance: ' + str(rec_moon_dist)
		print(str_moon_dist)
		
		rec_moon_illum = record['MOON_ILLUM']
		str_moon_illum = '\tMoon illumination: ' + str(rec_moon_illum)
		print(str_moon_illum)

		rec_gal_coords = record['GAL_COORDS']
		str_gal_coords = '\tGalactic coordinates: ' + str(rec_gal_coords)
		print(str_gal_coords)

		rec_ecl_coords = record['ECL_COORDS']
		str_ecl_coords = '\tEcliptic coordinates: ' + str(rec_ecl_coords)
		print(str_ecl_coords)

		rec_comments = record['COMMENTS']
		str_comments = '\tComments: ' + str(rec_comments)
		print(str_comments)

	@staticmethod
	def swift_bat_grb_lc(value, topic, alert=Configuration.ALERT):
		'''This function processes Swift BAT Lightcurve Notices from the Neil Gehrels Swift Observatory via GCN Classic over Kafka.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		record = Alerts.decode_classic_notice(value)

		rec_title = record['TITLE']
		str_title = '\tTitle: ' + str(rec_title)
		print(str_title)

		rec_notice_date = record['NOTICE_DATE']
		str_notice_date = '\tNotice date: ' + str(rec_notice_date)
		print(str_notice_date)

		rec_notice_type = record['NOTICE_TYPE']
		str_notice_type = '\tNotice type: ' + str(rec_notice_type)
		print(str_notice_type)

		rec_trigger_num = record['TRIGGER_NUM']
		str_trigger_num = '\tTrigger number: ' + str(rec_trigger_num)
		print(str_trigger_num)

		rec_grb_ra = record['GRB_RA']
		str_grb_ra = '\tGRB RA: ' + str(rec_grb_ra)
		print(str_grb_ra)

		rec_grb_dec = record['GRB_DEC']
		str_grb_dec = '\tGRB Dec: ' + str(rec_grb_dec)
		print(str_grb_dec)

		rec_grb_date = record['GRB_DATE']
		str_grb_date = '\tGRB date: ' + str(rec_grb_date)
		print(str_grb_date)

		rec_grb_time = record['GRB_TIME']
		str_grb_time = '\tGRB time: ' + str(rec_grb_time)
		print(str_grb_time)

		rec_trigger_index = record['TRIGGER_INDEX']
		str_trigger_index = '\tTrigger index: ' + str(rec_trigger_index)
		print(str_trigger_index)

		rec_grb_phi = record['GRB_PHI']
		str_grb_phi = '\tGRB phi: ' + str(rec_grb_phi)
		print(str_grb_phi)

		rec_grb_theta = record['GRB_THETA']
		str_grb_theta = '\tGRB theta: ' + str(rec_grb_theta)
		print(str_grb_theta)

		rec_delta_time = record['DELTA_TIME']
		str_delta_time = '\tDelta time: ' + str(rec_delta_time)
		print(str_delta_time)

		rec_lc_url = record['LC_URL']
		str_lc_url = '\tLightcurve URL: ' + str(rec_lc_url)
		print(str_lc_url)

		rec_sun_postn = record['SUN_POSTN']
		str_sun_postn = '\tSun position: ' + str(rec_sun_postn)
		print(str_sun_postn)
		
		rec_sun_dist = record['SUN_DIST']
		str_sun_dist = '\tSun distance: ' + str(rec_sun_dist)
		print(str_sun_dist)
		
		rec_moon_postn = record['MOON_POSTN']
		str_moon_postn = '\tMoon position: ' + str(rec_moon_postn)
		print(str_moon_postn)
		
		rec_moon_dist = record['MOON_DIST']
		str_moon_dist = '\tMoon distance: ' + str(rec_moon_dist)
		print(str_moon_dist)
		
		rec_moon_illum = record['MOON_ILLUM']
		str_moon_illum = '\tMoon illumination: ' + str(rec_moon_illum)
		print(str_moon_illum)
		
		rec_gal_coords = record['GAL_COORDS']
		str_gal_coords = '\tGalactic coordinates: ' + str(rec_gal_coords)
		print(str_gal_coords)
		
		rec_ecl_coords = record['ECL_COORDS']
		str_ecl_coords = '\tEcliptic coordinates: ' + str(rec_ecl_coords)
		print(str_ecl_coords)
		
		rec_comments = record['COMMENTS']
		str_comments = '\tComments: ' + str(rec_comments)
		print(str_comments)

	@staticmethod
	def swift_bat_grb_pos_ack(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_bat_grb_pos_test(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_bat_guano(value, topic, alert=Configuration.ALERT):

		record = json.loads(value)

		rec_schema = record['$schema']
		rec_mission = record['mission']
		rec_instrument = record['instrument']
		rec_messenger = record['messenger']
		rec_record_number = record['record_number']
		rec_alert_datetime = record['alert_datetime']
		rec_alert_tense = record['alert_tense']
		rec_trigger_time = record['trigger_time']
		rec_follow_up_event = record['follow_up_event']
		rec_follow_up_type = record['follow_up_type']
		rec_data_archive_page = record['data_archive_page']
		rec_id = record['id']

		# 1 - guano.example.json
		if record_number == 1:
			rec_rate_snr = record['rate_snr']
			rec_rate_duration = record['rate_duration']
			rec_rate_energy_range = record['rate_energy_range']
			rec_rate_energy_range_min = record['rate_energy_range'][0]
			rec_rate_energy_range_max = record['rate_energy_range'][1]
			rec_classification = record['classification']
			rec_classification_grb = record['classification']['GRB']
			rec_far = record['far']

		# 2 - guano.loc_map.example.json
		elif record_number == 2:
			rec_healpix_file = record['healpix_file']
			rec_systematic_included = record['systematic_included']
			rec_rate_snr = record['rate_snr']
			rec_rate_duration = record['rate_duration']
			rec_rate_energy_range = record['rate_energy_range']
			rec_rate_energy_range_min = record['rate_energy_range'][0]
			rec_rate_energy_range_max = record['rate_energy_range'][1]
			rec_classification = record['classification']
			rec_classification_grb = record['classification']['GRB']
			rec_far = record['far']

		# 3 - guano.loc_arc_min.example.json
		elif record_number == 3:
			rec_ra = record['ra']
			rec_dec = record['dec']
			rec_ra_dec_error = record['ra_dec_error']
			rec_containment_probability = record['containment_probability']
			rec_systematic_included = record['systematic_included']
			rec_rate_snr = record['rate_snr']
			rec_rate_duration = record['rate_duration']
			rec_rate_energy_range = record['rate_energy_range']
			rec_rate_energy_range_min = record['rate_energy_range'][0]
			rec_rate_energy_range_max = record['rate_energy_range'][1]
			rec_classification = record['classification']
			rec_classification_grb = record['classification']['GRB']
			rec_far = record['far']

		# 4 - guano.reraction.example.json
		elif record_number == 4:
			pass

		else:
			pass

	@staticmethod
	def swift_bat_ql_pos(value, topic, alert=Configuration.ALERT):
		'''This function processes Swift BAT QuickLook Position Notices from the Neil Gehrels Swift Observatory via GCN Classic over Kafka.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''
		
		record = Alerts.decode_classic_notice(value)

		rec_title = record['TITLE']
		str_title = '\tTitle: ' + str(rec_title)
		print(str_title)

		rec_notice_date = record['NOTICE_DATE']
		str_notice_date = '\tNotice date: ' + str(rec_notice_date)
		print(str_notice_date)

		rec_notice_type = record['NOTICE_TYPE']
		str_notice_type = '\tNotice type: ' + str(rec_notice_type)
		print(str_notice_type)

		rec_trigger_num = record['TRIGGER_NUM']
		str_trigger_num = '\tTrigger number: ' + str(rec_trigger_num)
		print(str_trigger_num)

		rec_grb_ra = record['GRB_RA']
		str_grb_ra = '\tGRB RA: ' + str(rec_grb_ra)
		print(str_grb_ra)

		rec_grb_dec = record['GRB_DEC']
		str_grb_dec = '\tGRB Dec: ' + str(rec_grb_dec)
		print(str_grb_dec)

		rec_grb_error = record['GRB_ERROR']
		str_grb_error = '\tGRB error: ' + str(rec_grb_error)
		print(str_grb_error)

		rec_grb_date = record['GRB_DATE']
		str_grb_date = '\tGRB date: ' + str(rec_grb_date)
		print(str_grb_date)

		rec_grb_time = record['GRB_TIME']
		str_grb_time = '\tGRB time: ' + str(rec_grb_time)
		print(str_grb_time)

		rec_rate_signif = record['RATE_SIGNIF']
		str_rate_signif = '\tRate significance: ' + str(rec_rate_signif)
		print(str_rate_signif)

		rec_soln_status = record['SOLN_STATUS']
		str_soln_status = '\tSoln status: ' + str(rec_soln_status)
		print(str_soln_status)

		rec_flags = record['FLAGS']
		str_flags = '\tFlags: ' + str(rec_flags)
		print(str_flags)

		rec_merit = record['MERIT']
		str_merit = '\tMerit: ' + str(rec_merit)
		print(str_merit)

		rec_sun_postn = record['SUN_POSTN']
		str_sun_postn = '\tSun position: ' + str(rec_sun_postn)
		print(str_sun_postn)
		
		rec_sun_dist = record['SUN_DIST']
		str_sun_dist = '\tSun distance: ' + str(rec_sun_dist)
		print(str_sun_dist)
		
		rec_moon_postn = record['MOON_POSTN']
		str_moon_postn = '\tMoon position: ' + str(rec_moon_postn)
		print(str_moon_postn)
		
		rec_moon_dist = record['MOON_DIST']
		str_moon_dist = '\tMoon distance: ' + str(rec_moon_dist)
		print(str_moon_dist)
		
		rec_moon_illum = record['MOON_ILLUM']
		str_moon_illum = '\tMoon illumination: ' + str(rec_moon_illum)
		print(str_moon_illum)
		
		rec_gal_coords = record['GAL_COORDS']
		str_gal_coords = '\tGalactic coordinates: ' + str(rec_gal_coords)
		print(str_gal_coords)
		
		rec_ecl_coords = record['ECL_COORDS']
		str_ecl_coords = '\tEcliptic coordinates: ' + str(rec_ecl_coords)
		print(str_ecl_coords)
		
		rec_comments = record['COMMENTS']
		str_comments = '\tComments: ' + str(rec_comments)
		print(str_comments)

	@staticmethod
	def swift_bat_scaledmap(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_bat_trans(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_fom_obs(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_pointdir(value, topic, alert=Configuration.ALERT):
		'''This function processes Swift Pointing Direction Notices from the Neil Gehrels Swift Observatory via GCN Classic over Kafka.

		:parameter value - the value of the message
		:parameter topic - the topic of the message
		:parameter alert - a toggle for broadcasting the alert

		:return - nothing is returned
		'''

		record = Alerts.decode_classic_notice(value)

		rec_title = record['TITLE']
		str_title = '\tTitle: ' + str(rec_title)
		print(str_title)

		rec_notice_date = record['NOTICE_DATE']
		str_notice_date = '\tNotice date: ' + str(rec_notice_date)
		print(str_notice_date)

		rec_notice_type = record['NOTICE_TYPE']
		str_notice_type = '\tNotice type: ' + str(rec_notice_type)
		print(str_notice_type)
		
		rec_next_point_ra = record['NEXT_POINT_RA']
		str_next_point_ra = '\tNext point (RA): ' + str(rec_next_point_ra)
		print(str_next_point_ra)
		
		rec_next_point_dec = record['NEXT_POINT_DEC']
		str_next_point_dec = '\tNext point (Dec): ' + str(rec_next_point_dec)
		print(str_next_point_dec)
		
		rec_next_point_roll = record['NEXT_POINT_ROLL']
		str_next_point_roll = '\tNext point (Roll): ' + str(rec_next_point_roll)
		print(str_next_point_roll)
		
		rec_slew_time = record['SLEW_TIME']
		str_slew_time = '\tSlew time: ' + str(rec_slew_time)
		print(str_slew_time)
		
		rec_slew_date = record['SLEW_DATE']
		str_slew_date = '\tSlew date: ' + str(rec_slew_date)
		print(str_slew_date)
		
		rec_obs_time = record['OBS_TIME']
		str_obs_time = '\tObservation time: ' + str(rec_obs_time)
		print(str_obs_time)
		
		rec_tgt_name = record['TGT_NAME']
		str_tgt_name = '\tTarget name: ' + str(rec_tgt_name)
		print(str_tgt_name)
		
		rec_tgt_num = record['TGT_NUM']
		str_tgt_num = '\tTarget number: ' + str(rec_tgt_num)
		print(str_tgt_num)
		
		rec_merit = record['MERIT']
		str_merit = '\tMerit: ' + str(rec_merit)
		print(str_merit)
		
		rec_inst_modes = record['INST_MODES']
		str_inst_modes = '\tInstrument modes: ' + str(rec_inst_modes)
		print(str_inst_modes)
		
		rec_sun_postn = record['SUN_POSTN']
		str_sun_postn = '\tSun position: ' + str(rec_sun_postn)
		print(str_sun_postn)
		
		rec_sun_dist = record['SUN_DIST']
		str_sun_dist = '\tSun distance: ' + str(rec_sun_dist)
		print(str_sun_dist)
		
		rec_moon_postn = record['MOON_POSTN']
		str_moon_postn = '\tMoon position: ' + str(rec_moon_postn)
		print(str_moon_postn)
		
		rec_moon_dist = record['MOON_DIST']
		str_moon_dist = '\tMoon distance: ' + str(rec_moon_dist)
		print(str_moon_dist)
		
		rec_moon_illum = record['MOON_ILLUM']
		str_moon_illum = '\tMoon illumination: ' + str(rec_moon_illum)
		print(str_moon_illum)
		
		rec_gal_coords = record['GAL_COORDS']
		str_gal_coords = '\tGalactic coordinates: ' + str(rec_gal_coords)
		print(str_gal_coords)
		
		rec_ecl_coords = record['ECL_COORDS']
		str_ecl_coords = '\tEcliptic coordinates: ' + str(rec_ecl_coords)
		print(str_ecl_coords)
		
		rec_comments = record['COMMENTS']
		str_comments = '\tComments: ' + str(rec_comments)
		print(str_comments)

	@staticmethod
	def swift_sc_slew(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_too_fom(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_too_sc_slew(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_uvot_dburst(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_uvot_dburst_proc(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_uvot_emergency(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_uvot_fchart(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_uvot_fchart_proc(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_uvot_pos(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_uvot_pos_nack(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_xrt_centroid(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_xrt_emergency(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_xrt_image(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_xrt_image_proc(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_xrt_lc(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_xrt_position(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_xrt_spectrum(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_xrt_spectrum_proc(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_xrt_sper(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_xrt_sper_proc(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_xrt_threshpix(value, topic, alert=Configuration.ALERT):
		return

	@staticmethod
	def swift_xrt_threshpix_proc(value, topic, alert=Configuration.ALERT):
		return