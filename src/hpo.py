#!/usr/bin/env python3

import os.path
import urllib.request
import io
import csv
from model.network import Network
from model.disease import Disease
from model.gene import Gene

file = '../data/HPO/OMIM_FREQUENT_FEATURES_diseases_to_genes_to_phenotypes.txt'
url = 'http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/' + \
      'OMIM_FREQUENT_FEATURES_diseases_to_genes_to_phenotypes.txt'

if not os.path.exists(file):
    print('Database does not exist. Trying to download...')
    with urllib.request.urlopen(url) as response, open(file, 'wb') as f:
        f.write(response.read())

network = Network()

with io.open(file, 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    next(reader, None)
    for row in reader:
        disease = Disease([row[0]], [])
        network.add_node(disease)
        gene = Gene(['HGNC:%s' % row[1], 'Entrez:%s' % row[2]], [])
        network.add_node(gene)
        hpo_id = row[3]
        hpo_term_name = row[4]
        # TODO

network.save('../data/HPO/graph.json')
