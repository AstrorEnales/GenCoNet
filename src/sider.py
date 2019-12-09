#!/usr/bin/env python3

import os.path
import urllib.request
import io
import csv

from model.edge import Edge
from model.network import Network
from model.drug import Drug
from model.disease import Disease
import gzip
import shutil

file = '../data/SIDER/meddra_all_indications.tsv'
drug_file = '../data/SIDER/drug_names.tsv'
zip_file = '../data/SIDER/meddra_all_indications.tsv.gz'
url = 'http://sideeffects.embl.de/media/download/meddra_all_indications.tsv.gz'
drug_url = 'http://sideeffects.embl.de/media/download/drug_names.tsv'

if not os.path.exists(file) or not os.path.exists(drug_file):
    print('Database does not exist. Trying to download and extract...')
    if not os.path.exists(zip_file):
        print('Downloading latest archive...')
        with urllib.request.urlopen(url) as response, open(zip_file, 'wb') as f:
            f.write(response.read())
    if not os.path.exists(drug_file):
        with urllib.request.urlopen(drug_url) as response, open(drug_file, 'wb') as f:
            f.write(response.read())
    print('Extracting database file...')
    with gzip.open(zip_file, 'rb') as z:
        with open(file, 'wb') as f:
            shutil.copyfileobj(z, f)

network = Network()

drug_lookup = {}
with io.open(drug_file, 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    for row in reader:
        drug_lookup[row[0].strip()] = row[1].strip()

# 1: STITCH compound id (flat, see above)
# 2: UMLS concept id as it was found on the label
# 3: method of detection: NLP_indication / NLP_precondition / text_mention
# 4: concept name
# 5: MedDRA concept type (LLT = lowest level term, PT = preferred term; in a few cases the term is neither LLT nor PT)
# 6: UMLS concept id for MedDRA term
# 7: MedDRA concept name

# All side effects found on the labels are given as LLT. Additionally, the PT is shown. There is at least one
# PT for every LLT, but sometimes the PT is the same as the LLT.
with io.open(file, 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    for row in reader:
        pubchem_id = row[0][4::].lstrip('0')
        drug = Drug(['PubChem:CID%s' % pubchem_id], [drug_lookup[row[0]]] if row[0] in drug_lookup else [])
        network.add_node(drug)
        disease = Disease(['UMLS:%s' % row[1], 'UMLS:%s' % row[5]], [row[3], row[6]])
        network.add_node(disease)
        network.add_edge(Edge(next(iter(drug.ids)), next(iter(disease.ids)), 'INDICATES', {'source': 'SIDER'}))

network.save('../data/SIDER/graph.json')
