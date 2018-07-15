#!/usr/bin/env python3

import os.path
import urllib.request

file = '../data/HPO/OMIM_FREQUENT_FEATURES_diseases_to_genes_to_phenotypes.txt'
url = 'http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/' + \
      'OMIM_FREQUENT_FEATURES_diseases_to_genes_to_phenotypes.txt'

if not os.path.exists(file):
    print('Database does not exist. Trying to download...')
    with urllib.request.urlopen(url) as response, open(file, 'wb') as f:
        f.write(response.read())
