#!/usr/bin/env python3

import os.path
import urllib.request
import gzip
import io
import csv
from model.network import Network
from model.variant import Variant
from model.disease import Disease
from model.gene import Gene
from model.edge import Edge

gene_file = '../data/DisGeNet/curated_gene_disease_associations.tsv'
gene_zip_file = '../data/DisGeNet/curated_gene_disease_associations.tsv.gz'
gene_url = 'https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz'
variant_file = '../data/DisGeNet/curated_variant_disease_associations.tsv'
variant_zip_file = '../data/DisGeNet/curated_variant_disease_associations.tsv.gz'
variant_url = 'https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_variant_disease_associations.tsv.gz'

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
        #  0 - geneId
        #  1 - geneSymbol
        #  2 - DSI
        #  3 - DPI
        #  4 - diseaseId
        #  5 - diseaseName
        #  6 - diseaseType
        #  7 - diseaseClass
        #  8 - diseaseSemanticType
        #  9 - score
        # 10 - EI
        # 11 - YearInitial
        # 12 - YearFinal
        # 13 - NofPmids
        # 14 - NofSnps
        # 15 - source
        if int(row[13]) >= PUBMED_COUNT_THRESHOLD:
            gene = Gene(['HGNC:%s' % row[1]], [])
            network.add_node(gene)
            disease = Disease(['UMLS:%s' % row[4]], [row[5]])
            network.add_node(disease)
            rel = {
                'source': 'DisGeNet,%s' % row[15],
                'num_pmids': int(row[13]),
                'num_snps': int(row[14]),
                'score': row[9]
            }
            network.add_edge(Edge(next(iter(gene.ids)), next(iter(disease.ids)), 'ASSOCIATES_WITH', rel))

with io.open('../data/DisGeNet/curated_variant_disease_associations.tsv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    next(reader, None)
    for row in reader:
        # 0 snpId
        # 1 diseaseId
        # 2 diseaseName
        # 3 score
        # 4 NofPmids
        # 5 source

        #  0 - snpId
        #  1 - chromosome
        #  2 - position
        #  3 - DSI
        #  4 - DPI
        #  5 - diseaseId
        #  6 - diseaseName
        #  7 - diseaseType
        #  8 - diseaseClass
        #  9 - diseaseSemanticType
        # 10 - score
        # 11 - EI
        # 12 - YearInitial
        # 13 - YearFinal
        # 14 - NofPmids
        # 15 - source
        if int(row[14]) >= PUBMED_COUNT_THRESHOLD:
            variant = Variant(['dbSNP:%s' % row[0]], [])
            network.add_node(variant)
            disease = Disease(['UMLS:%s' % row[5]], [row[6]])
            network.add_node(disease)
            rel = {
                'source': 'DisGeNet,%s' % row[15],
                'num_pmids': int(row[14]),
                'score': row[10]
            }
            network.add_edge(Edge(next(iter(variant.ids)), next(iter(disease.ids)), 'ASSOCIATES_WITH', rel))

network.save('../data/DisGeNet/graph.json')
