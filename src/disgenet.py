#!/usr/bin/env python3

import os.path
import urllib.request
import gzip

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
