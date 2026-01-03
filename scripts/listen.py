from config import Configuration
from libraries.router import Router
from libraries.utils import Utils
from gcn_kafka import Consumer
import logging
logging.getLogger('gcn').setLevel(logging.WARNING)

consumer = Consumer(client_id=Configuration.CLIENT_ID, client_secret=Configuration.CLIENT_SECRET)

consumer.subscribe(Configuration.AVAILABLE_TOPICS)

while True:
    for message in consumer.consume(timeout=1):
        message_error = message.error()
        message_topic = message.topic()
        message.value = message.value()

        if message_error:
            Utils.log(message_error, 'debug')
            Utils.log(message_topic, 'debug')
            Utils.log(message_value, 'debug')

        else:
            Router.route_alert(message_value, message_topic)