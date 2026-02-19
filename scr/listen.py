from config import Configuration
from lib.routing import Router
from lib.utilities import Utils
from gcn_kafka import Consumer
import os

os.system('clear')
Utils.log('Beginning listener for astronomical alerts.', 'info')

consumer = Consumer(client_id=Configuration.CLIENT_ID, client_secret=Configuration.CLIENT_SECRET)
consumer.subscribe(Configuration.AVAILABLE_TOPICS)

while True:
	for message in consumer.consume(timeout=1):
		message_error = message.error()
		message_offset = message.offset()
		message_topic = message.topic()
		message_value = message.value()

		if message_error:
			Utils.log(message_error, 'debug')
			Utils.log(message_offset, 'debug')
			Utils.log(message_topic, 'debug')
			Utils.log(message_value, 'debug')

		else:
			Router.filter_alert(message_value, message_topic)