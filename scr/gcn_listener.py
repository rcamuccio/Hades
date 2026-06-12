from config import Configuration
from lib.alerts import Alerts

from gcn_kafka import Consumer
import os

os.system('clear')
print('[\033[1m' + 'ᾍδης ζῇ' + '\033[0m] - Running GCN listener\n')

consumer = Consumer(client_id=Configuration.CLIENT_ID, client_secret=Configuration.CLIENT_SECRET)
consumer.subscribe(Configuration.AVAILABLE_TOPICS)

while True:
	for message in consumer.consume(timeout=1):
		message_error = message.error()
		message_offset = message.offset()
		message_topic = message.topic()
		message_value = message.value()

		if message_error:
			print()
			print(message_error)
			print(message_offset)
			print(message_topic)
			print(message_value)
			print()

		else:
			Alerts.filter_alert(message_value, message_topic)