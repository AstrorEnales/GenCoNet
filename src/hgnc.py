#!/usr/bin/env python3

import os.path
import urllib.request
import io
import csv
from model.network import Network
from model.gene import Gene

file = '../data/HGNC/hgnc_complete_set.txt'
url = 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt'

if not os.path.exists(file):
    print('Database does not exist. Trying to download...')
    with urllib.request.urlopen(url) as response, open(file, 'wb') as f:
        f.write(response.read())

network = Network()

with io.open(file, 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    next(reader, None)
    for row in reader:
        gene_ids = ['HGNC:%s' % row[1]]
        if row[0]:
            gene_ids.append(row[0])
        network.add_node(Gene(gene_ids, [row[2]]))

network.save('../data/HGNC/graph.json')
