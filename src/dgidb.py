#!/usr/bin/env python3

import json
import os.path
import urllib.request
import io
import csv
from model.network import Network
from model.drug import Drug
from model.gene import Gene
from model.edge import Edge

file = '../data/DGIdb/interactions.tsv'
url = 'http://www.dgidb.org/data/interactions.tsv'

if not os.path.exists(file):
    print('Database does not exist. Trying to download...')
    with urllib.request.urlopen(url) as response, open(file, 'wb') as f:
        f.write(response.read())

network = Network()

with io.open(file, 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    next(reader, None)
    for row in reader:
        row = [x.strip() for x in row]
        if row[0] is None or len(row[0]) == 0:
            continue
        if row[7] is None or len(row[7]) == 0:
            continue
        if row[8] is None or len(row[8]) == 0:
            continue
        gene_ids = {'HGNCSymbol:%s' % row[0]}
        if row[2] is not None and len(row[2]) > 0:
            gene_ids.add('Entrez:%s' % row[2])
        gene = Gene(gene_ids, [])
        network.add_node(gene)
        drug = Drug(['ChEMBL:%s' % row[8]], [row[7]])
        network.add_node(drug)
        rel = {
            'source': 'DGIdb,%s' % row[3],
            'actions': [row[4]],
        }
        if row[9] is not None and len(row[9]) > 0:
            pubmed_ids = ','.join(['PMID:%s' % x for x in row[9].strip().split(',')])
            rel['source'] += ',%s' % pubmed_ids
        network.add_edge(Edge(next(iter(drug.ids)), next(iter(gene.ids)), 'TARGETS', rel))

with io.open('../data/DGIdb/graph.json', 'w', encoding='utf-8', newline='') as f:
    f.write(json.dumps(network.to_dict(), indent=2))
