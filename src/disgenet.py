#!/usr/bin/env python3

import os.path
import urllib.request
import gzip
import io
import csv
import json
from model.network import Network
from model.variant import Variant
from model.disease import Disease
from model.gene import Gene
from model.edge import Edge

gene_file = '../data/DisGeNet/curated_gene_disease_associations.tsv'
gene_zip_file = '../data/DisGeNet/curated_gene_disease_associations.tsv.gz'
gene_url = 'http://www.disgenet.org/ds/DisGeNET/results/curated_gene_disease_associations.tsv.gz'
variant_file = '../data/DisGeNet/curated_variant_disease_associations.tsv'
variant_zip_file = '../data/DisGeNet/curated_variant_disease_associations.tsv.gz'
variant_url = 'http://www.disgenet.org/ds/DisGeNET/results/curated_variant_disease_associations.tsv.gz'

for url, file, zip_file in [(gene_url, gene_file, gene_zip_file), (variant_url, variant_file, variant_zip_file)]:
    if not os.path.exists(file):
        print('Database does not exist. Trying to download and extract...')
        if not os.path.exists(zip_file):
            print('Downloading latest archive...')
            with urllib.request.urlopen(url) as response, open(zip_file, 'wb') as f:
                f.write(response.read())
        print('Extracting database file...')
        with gzip.open(zip_file, 'rb') as f:
            with open(file, 'wb') as out_file:
                out_file.write(f.read())

network = Network()

PUBMED_COUNT_THRESHOLD = 2
with io.open('../data/DisGeNet/curated_gene_disease_associations.tsv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    next(reader, None)
    for row in reader:
        # geneId
        # geneSymbol
        # diseaseId
        # diseaseName
        # score
        # NofPmids
        # NofSnps
        # source
        if int(row[5]) >= PUBMED_COUNT_THRESHOLD:
            gene = Gene(['HGNCSymbol:%s' % row[1]], [])
            network.add_node(gene)
            disease = Disease(['UMLS:%s' % row[2]], [row[3]])
            network.add_node(disease)
            rel = {
                'source': 'DisGeNet,%s' % row[7],
                'num_pmids': int(row[5]),
                'num_snps': int(row[6]),
                'score': row[4]
            }
            network.add_edge(Edge(next(iter(gene.ids)), next(iter(disease.ids)), 'ASSOCIATES_WITH', rel))

with io.open('../data/DisGeNet/curated_variant_disease_associations.tsv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    next(reader, None)
    for row in reader:
        # snpId
        # diseaseId
        # diseaseName
        # score
        # NofPmids
        # source
        if int(row[4]) >= PUBMED_COUNT_THRESHOLD:
            variant = Variant(['dbSNP:%s' % row[0]], [])
            network.add_node(variant)
            disease = Disease(['UMLS:%s' % row[1]], [row[2]])
            network.add_node(disease)
            rel = {
                'source': 'DisGeNet,%s' % row[5],
                'num_pmids': int(row[4]),
                'score': row[3]
            }
            network.add_edge(Edge(next(iter(variant.ids)), next(iter(disease.ids)), 'ASSOCIATES_WITH', rel))

with io.open('../data/DisGeNet/graph.json', 'w', encoding='utf-8', newline='') as f:
    f.write(json.dumps(network.to_dict(), indent=2))
