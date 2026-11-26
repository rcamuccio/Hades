import numpy as np
import os
import sys
import time
from gcn_kafka import Consumer

from conf import Configuration
from lib.angelos import Angelos
from lib.boethos import Boethos

# terminal preamble
os.system('clear')
script_preamble = Boethos.get_script_preamble('gcn_listener')
print(script_preamble)

st = time.time()

# create a consumer object
consumer = Consumer(client_id=Configuration.CLIENT_ID, client_secret=Configuration.CLIENT_SECRET, **{'log_level': 0})

# get the available notice topics
topic_list, n_topics = Angelos.get_notice_topics()
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
				Angelos.filter_notice(msg_val, msg_top)

	# handle manual interruptions
	except KeyboardInterrupt:
		script_postamble = Boethos.get_script_postamble('gcn_listener')
		print(script_postamble)

		fn = time.time()
		dt = np.around(fn - st, decimals=2)
		print('[' + str(dt) + ' s]')

		try:
			sys.exit(130)
		except SystemExit:
			os._exit(130)

	# handle unplanned interruptions
	except Exception as e:
		print(f'{type(e).__name__}: {e}')