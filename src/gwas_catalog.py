#!/usr/bin/env python3

import os.path
import urllib.request

file = '../data/GWAS-Catalog/gwas_catalog_associations.tsv'
url = 'https://www.ebi.ac.uk/gwas/api/search/downloads/full'

if not os.path.exists(file):
    print('Database does not exist. Trying to download...')
    with urllib.request.urlopen(url) as response, open(file, 'wb') as f:
        f.write(response.read())
