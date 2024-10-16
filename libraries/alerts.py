from config import Configuration
from libraries.utils import Utils

from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import requests
import smtplib

import logging
logging.getLogger('requests').setLevel(logging.WARNING)
logging.getLogger('urllib3').setLevel(logging.WARNING)

class Alerts:

	@staticmethod
	def alert_team(alert_type='test', event_name='test'):
		''' This function sends text messages and emails to the team when an alert is detected.

		:parameter alert_type - The example alert type for the text to send; the default is to simply send a test alert.

		:return - Nothing is returned, howerver a message is sent to devices.

		'''

		# log in to the email server using the appropriate credentials
		server = smtplib.SMTP_SSL(Configuration.SMTP, Configuration.PORT)
		server.ehlo()
		server.login(Configuration.EMAIL, Configuration.PAS)

		# set up the message based on the type of alert
		msg = MIMEMultipart()
		msg['From'] = Configuration.EMAIL
		msg['To'] = ', '.join(Configuration.MAILING_LIST)

		# filter alert types
		if alert_type == 'PRELIMINARY':
			msg['Subject'] = 'Preliminary LVK-GW alert: ' + event_name + '\n'
			body = 'Ole! A new GW-event has been detected by LVK. Check queue and results.\n'

		elif alert_type == 'INITIAL':
			msg['Subject'] = 'Initial LVK-GW alert: ' + event_name + '\n'
			body = 'Ole! A new GW-event has been detected by LVK. Check queue and results.\n'
		
		elif alert_type == 'RETRACTION':
			msg['Subject'] = 'Event: ' + event_name + ' retracted\n'
			body = 'The previous alert has been retracted.\n'
		
		elif alert_type == 'OBSERVER':
			msg['Subject'] = 'Observer not present\n'
			body = 'A new alert has been detected, but the observer may not be present.\n'
		
		else:
			msg['Subject'] = 'Alert Test for event: ' + event_name + '\n'
			body = 'This is a test of the alert message system.\n'

		# convert message to text and send
		msg.attach(MIMEText(body, 'plain'))
		sms = msg.as_string()
		server.sendmail(Configuration.EMAIL, Configuration.MAILING_LIST, sms)

		# clean up
		server.quit()

		# log that the message was sent
		Utils.log('The following alert message was sent to the team: ' + body, 'info')

		return