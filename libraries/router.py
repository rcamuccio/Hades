from config import Configuration
from libraries.alerts import Alerts
from libraries.priority import Priority
from libraries.utils import Utils
from astropy import units as u
from astropy.table import QTable, Table
from base64 import b64decode
from io import BytesIO
import astropy_healpix as ah
import json
import numpy as np
import os

class Router:

    @staticmethod
    def decode_classic_notice(value, verbose=False):
        '''This function decodes classic GCN notices.

        :parameter value - 
        :parameter verbose -

        :return params -
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

                if verbose:
                    print(key, item)

        params = dict(zip(keys, item))

        if verbose:
            print('----------')

        return params