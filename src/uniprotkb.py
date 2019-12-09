#!/usr/bin/env python3

import os.path
import urllib.request
import gzip
import io
import csv
from model.network import Network
from model.gene import Gene

id_mapping_file = '../data/UniprotKB/HUMAN_9606_idmapping.dat'
id_mapping_zip_file = '../data/UniprotKB/HUMAN_9606_idmapping.dat.gz'
id_mapping_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz'

if not os.path.exists(id_mapping_file):
    print('Database does not exist. Trying to download and extract...')
    if not os.path.exists(id_mapping_zip_file):
        print('Downloading latest archive...')
        with urllib.request.urlopen(id_mapping_url) as response, open(id_mapping_zip_file, 'wb') as f:
            f.write(response.read())
    print('Extracting database file...')
    with gzip.open(id_mapping_zip_file, 'rb') as f:
        with open(id_mapping_file, 'wb') as out_file:
            out_file.write(f.read())

network = Network()

with io.open(id_mapping_file, 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    for row in reader:
        if row[1] == 'HGNC':
            gene = Gene(['UniProtKB:%s' % row[0], row[2]], [])
            network.add_node(gene)
        elif row[1] == 'Gene_Name':
            gene = Gene(['UniProtKB:%s' % row[0], 'HGNC:%s' % row[2]], [])
            network.add_node(gene)

network.save('../data/UniprotKB/graph.json')
