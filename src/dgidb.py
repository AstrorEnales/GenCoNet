#!/usr/bin/env python3

import io
import csv
import os.path
import urllib.request
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
        if not row[0] or not row[7] or not row[8]:
            continue
        gene_ids = {'HGNC:%s' % row[0]}
        if row[2]:
            gene_ids.add('Entrez:%s' % row[2])
        gene = Gene(gene_ids, [])
        network.add_node(gene)
        drug_name = row[7].replace('(%s)' % row[8], '').replace(row[8], '').strip()
        drug = Drug(['ChEMBL:%s' % row[8]], [drug_name] if drug_name else [])
        network.add_node(drug)
        rel = {
            'source': 'DGIdb,%s' % row[3],
            'actions': [row[4]],
        }
        if row[9]:
            pubmed_ids = ','.join(['PMID:%s' % x for x in row[9].strip().split(',')])
            rel['source'] += ',%s' % pubmed_ids
        network.add_edge(Edge(next(iter(drug.ids)), next(iter(gene.ids)), 'TARGETS', rel))

network.save('../data/DGIdb/graph.json')
