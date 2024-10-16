from config import Configuration
from libraries.router import Router
from libraries.utils import Utils

from gcn_kafka import Consumer

import logging
logging.getLogger('gcn').setLevel(logging.WARNING)

consumer = Consumer(client_id=Configuration.CLIENT_ID,
					client_secret=Configuration.CLIENT_SECRET)

consumer.subscribe(Configuration.AVAILABLE_TOPICS)

Utils.setup_hades(Configuration.MAIN_DIR)

while True:

	for message in consumer.consume(timeout=1):

		message_error = message.error()
		message_topic = message.topic()
		message_value = message.value()

		if message_error:
			Utils.log(message_error, 'info')
			Utils.log(message_topic, 'info')
			Utils.log(message_value, 'info')

		else:
			Router.route_alert(message_value, message_topic)