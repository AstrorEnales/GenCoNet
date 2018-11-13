#!/usr/bin/env python3

import json
import os.path
import urllib.request
import io
import csv
from model.network import Network
from model.drug import Drug
import gzip
import shutil

file = '../data/SIDER/meddra_all_indications.tsv'
zip_file = '../data/SIDER/meddra_all_indications.tsv.gz'
url = 'http://sideeffects.embl.de/media/download/meddra_all_indications.tsv.gz'

if not os.path.exists(file):
    print('Database does not exist. Trying to download and extract...')
    if not os.path.exists(zip_file):
        print('Downloading latest archive...')
        with urllib.request.urlopen(url) as response, open(zip_file, 'wb') as f:
            f.write(response.read())
    print('Extracting database file...')
    with gzip.open(zip_file, 'rb') as z:
        with open(file, 'wb') as f:
            shutil.copyfileobj(z, f)

network = Network()

with io.open('../data/SIDER/graph.json', 'w', encoding='utf-8', newline='') as f:
    f.write(json.dumps(network.to_dict(), indent=2))
