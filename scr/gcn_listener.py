from config import Configuration
from lib.alerts import Alerts
from gcn_kafka import Consumer
import os
import sys

os.system('clear')
print('[\033[1m' + 'ᾍδης ζῇ' + '\033[0m] - Running GCN listener\n')

# create a consumer object
consumer = Consumer(client_id=Configuration.CLIENT_ID, client_secret=Configuration.CLIENT_SECRET, **{'log_level': 0})

# get the available notice topics
topic_list, n_topics = Alerts.get_notice_topics()
consumer.subscribe(topic_list)

while True:
	try:
		for msg in consumer.consume(timeout=1):
			msg_err = msg.error()
			msg_off = msg.offset()
			msg_top = msg.topic()
			msg_val = msg.value()

			if msg_err:
				print(msg_err)
				print(msg_off)
				print(msg_top)
				print(msg_val)

			else:
				Alerts.filter_notice(msg_val, msg_top)

	# handle manual interruptions
	except KeyboardInterrupt:
		print('\n\n[\033[1m' + 'ᾍδης ἀπέρχεται' + '\033[0m] - GCN listener terminated')
		try:
			sys.exit(130)
		except SystemExit:
			os._exit(130)

	# handle unplanned interruptions
	except Exception as e:
		print(f'{type(e).__name__}: {e}')