#!/usr/bin/env python3

import io
import os
import urllib.request
import urllib.parse

from model.edge import Edge
from model.network import Network
from model.gene import Gene
from model.disease import Disease

file = '../data/HuGE-Navigator/Disease-GeneID.txt'
url = 'https://phgkb.cdc.gov/PHGKB/fileDownload.action'
if not os.path.exists(file):
    print('Database does not exist. Trying to download...')
    data = urllib.parse.urlencode({'downLoadType': 'all_pheno', 'Mysubmit': 'Download'}).encode()
    with urllib.request.urlopen(urllib.request.Request(url, data=data)) as response:
        with open(file, 'wb') as f:
            f.write(response.read())

network = Network()

with io.open(file, 'r', encoding='utf-8', newline='') as f:
    for skip in range(0, 4):
        f.readline()
    for line in f:
        parts = [[y.strip() for y in x[:-1].split('(')] for x in line.strip().split('\t')]
        if len(parts) > 1:
            disease = Disease(['UMLS:%s' % parts[0][1]], [parts[0][0]])
            network.add_node(disease)
            for part in parts[1::]:
                gene = Gene(['HGNC:%s' % part[0], 'Entrez:%s' % part[1]], [part[0]])
                network.add_node(gene)
                rel = {'source': 'HuGE Navigator'}
                network.add_edge(Edge(gene, disease, 'ASSOCIATES_WITH', rel))

network.save('../data/HuGE-Navigator/graph.json')
