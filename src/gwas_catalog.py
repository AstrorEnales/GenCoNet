#!/usr/bin/env python3

import re
import io
import csv
import json
import os.path
import urllib.request
from model.network import Network
from model.variant import Variant
from model.gene import Gene
from model.edge import Edge

file = '../data/GWAS-Catalog/gwas_catalog_associations.tsv'
url = 'https://www.ebi.ac.uk/gwas/api/search/downloads/full'

if not os.path.exists(file):
    print('Database does not exist. Trying to download...')
    with urllib.request.urlopen(url) as response, open(file, 'wb') as f:
        f.write(response.read())

network = Network()

# 0 DATE ADDED TO CATALOG
# 1 PUBMEDID
# 2 FIRST AUTHOR
# 3 DATE
# 4 JOURNAL
# 5 LINK
# 6 STUDY
# 7 DISEASE/TRAIT
# 8 INITIAL SAMPLE SIZE
# 9 REPLICATION SAMPLE SIZE
# 10 REGION
# 11 CHR_ID
# 12 CHR_POS
# 13 REPORTED GENE(S)
# 14 MAPPED_GENE
# 15 UPSTREAM_GENE_ID
# 16 DOWNSTREAM_GENE_ID
# 17 SNP_GENE_IDS
# 18 UPSTREAM_GENE_DISTANCE
# 19 DOWNSTREAM_GENE_DISTANCE
# 20 STRONGEST SNP-RISK ALLELE
# 21 SNPS
# 22 MERGED
# 23 SNP_ID_CURRENT
# 24 CONTEXT
# 25 INTERGENIC
# 26 RISK ALLELE FREQUENCY
# 27 P-VALUE
# 28 PVALUE_MLOG
# 29 P-VALUE (TEXT)
# 30 OR or BETA
# 31 95% CI (TEXT)
# 32 PLATFORM [SNPS PASSING QC]
# 33 CNV
loc_pattern = re.compile(r'LOC[0-9]+')
with io.open('../data/GWAS-Catalog/gwas_catalog_associations.tsv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    next(reader, None)
    for row in reader:
        if row[14] is None or len(row[14]) == 0:
            continue
        gene_ids = row[14].replace(' x ', ', ').replace(' - ', ', ').split(', ')
        for gene_id in gene_ids:
            if loc_pattern.fullmatch(gene_id) is not None:
                continue
            gene = Gene(['HGNCSymbol:%s' % gene_id], [])
            network.add_node(gene)
            variant = Variant(['dbSNP:%s' % row[21]], [])
            network.add_node(variant)
            network.add_edge(Edge(gene.id, variant.id, 'CODES', {'source': 'GWASCatalog', 'pmid': row[1]}))

with io.open('../data/GWAS-Catalog/graph.json', 'w', encoding='utf-8', newline='') as f:
    f.write(json.dumps(network.to_dict(), indent=2))