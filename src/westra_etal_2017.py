#!/usr/bin/env python3

import json
import os.path
import io
import csv
import urllib.request
import zipfile
from model.network import Network
from model.gene import Gene
from model.variant import Variant
from model.edge import Edge

file_cis = '../data/Westra_etal_2017/2012-12-21-CisAssociationsProbeLevelFDR0.5.txt'
file_trans = '../data/Westra_etal_2017/2012-12-21-TransEQTLsFDR0.5.txt'
zip_file_cis = '../data/Westra_etal_2017/2012-12-21-CisAssociationsProbeLevelFDR0.5.zip'
zip_file_trans = '../data/Westra_etal_2017/2012-12-21-TransEQTLsFDR0.5.zip'
url_cis = 'https://www.genenetwork.nl/bloodeqtlbrowser/2012-12-21-CisAssociationsProbeLevelFDR0.5.zip'
url_trans = 'https://www.genenetwork.nl/bloodeqtlbrowser/2012-12-21-TransEQTLsFDR0.5.zip'
if not os.path.exists(file_cis) or not os.path.exists(file_trans):
    print('Database does not exist. Trying to download and extract...')
    if not os.path.exists(zip_file_cis):
        print('Downloading latest archive...')
        with urllib.request.urlopen(url_cis) as response, open(zip_file_cis, 'wb') as f:
            f.write(response.read())
    if not os.path.exists(zip_file_trans):
        print('Downloading latest archive...')
        with urllib.request.urlopen(url_trans) as response, open(zip_file_trans, 'wb') as f:
            f.write(response.read())
    print('Extracting database file...')
    with zipfile.ZipFile(zip_file_cis) as z:
        with open(file_cis, 'wb') as f:
            f.write(z.read(z.namelist()[0]))
    with zipfile.ZipFile(zip_file_trans) as z:
        with open(file_trans, 'wb') as f:
            f.write(z.read(z.namelist()[0]))


def value_empty(s: str) -> bool:
    return not s or s.strip() == '-'


network = Network()

with io.open(file_cis, 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    next(reader, None)
    for row in reader:
        if value_empty(row[1]) or value_empty(row[13]):
            continue
        variant = Variant(['dbSNP:%s' % row[1]], [])
        network.add_node(variant)
        for gene_id in row[13].split(','):
            gene = Gene(['HGNC:%s' % gene_id], [])
            network.add_node(gene)
            rel = {
                'source': 'PMID:24013639',
                'pvalue': row[0],
                'snp_chr': row[2],
                'cis_trans': row[7]
            }
            network.add_edge(Edge(next(iter(gene.ids)), next(iter(variant.ids)), 'EQTL', rel))

with io.open(file_trans, 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    next(reader, None)
    for row in reader:
        if value_empty(row[1]) or value_empty(row[16]):
            continue
        variant = Variant(['dbSNP:%s' % row[1]], [])
        network.add_node(variant)
        for gene_id in row[16].split(','):
            gene = Gene(['HGNC:%s' % gene_id], [])
            network.add_node(gene)
            rel = {
                'source': 'PMID:24013639',
                'pvalue': row[0],
                'snp_chr': row[2],
                'cis_trans': 'trans'  # row[7] TODO: always "-"?
            }
            network.add_edge(Edge(next(iter(gene.ids)), next(iter(variant.ids)), 'EQTL', rel))

with io.open('../data/Westra_etal_2017/graph.json', 'w', encoding='utf-8', newline='') as f:
    f.write(json.dumps(network.to_dict(), indent=2))
