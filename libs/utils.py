import os
import re
import logging

def dir_check(dirname):
    if os.path.exists(dirname):
        return
    os.makedirs(dirname)

def file_check(filename):
    return os.path.isfile(filename)

def read_config(configfile):
    config = {}
    with open(configfile) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            if re.match(r'^\s$', line):
                continue
            line = re.sub(r'[\r\n\s]+','',line)
            arr = line.split('=')
            config[arr[0]] = arr[1]
    logging.info('configuration info loaded.')
    return config
