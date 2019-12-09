#!/usr/bin/env python3

import os.path
import urllib.request
import io
import csv
from model.network import Network
from model.drug import Drug
import gzip
import shutil

file = '../data/SuperDrug2/all_drugs_extlinks.csv'
zip_file = '../data/SuperDrug2/all_drugs_extlinks.csv.gz'
url = 'http://cheminfo.charite.de/superdrug2/downloads/zip_files/all_drugs_extlinks.csv.gz'

if not os.path.exists(file):
    print('Database does not exist. Trying to download and extract...')
    if not os.path.exists(zip_file):
        print('Downloading latest archive...')
        with urllib.request.urlopen(url) as response, open(zip_file, 'wb') as f:
            f.write(response.read())
    print('Extracting database file...')
    with gzip.open(zip_file, 'rb') as z:
        with open(file, 'wb') as f:
            shutil.copyfileobj(z, f)

network = Network()

not_available_text = 'Not Available'

with io.open(file, 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter=',', quotechar='"')
    next(reader, None)
    for row in reader:
        # 0 SUPERDRUG_ID
        # 1 PREFERRED_NAME
        # 2 ATC
        # 3 CHEMBL_ID
        # 4 DRUGBANK_ID
        # 5 KEGG_ID
        # 6 PUBCHEM_CID
        # 7 CASRN
        drug_ids = []  # ['SuperDrug:%s' % row[0]]
        # if row[2] != not_available_text:
        #     drug_ids.extend(['AtcCode:%s' % x for x in row[2].split(';')])
        if row[3] != not_available_text:
            for chembl_id in {x.strip() for x in row[3].split(';')}:
                drug_ids.append('ChEMBL:%s' % chembl_id)
        if row[4] != not_available_text:
            drug_ids.append('DrugBank:%s' % row[4])
        else:
            # For now, only use mappings including DrugBank
            continue
        # if row[5] != not_available_text:
        #     drug_ids.append('Kegg:%s' % row[5])
        if row[6] != not_available_text:
            for pubchem_id in {x.strip() for x in row[6].split(';')}:
                drug_ids.append('PubChem:CID%s' % pubchem_id)
        if len(drug_ids) > 1:
            network.add_node(Drug(drug_ids, [row[1]]))

network.save('../data/SuperDrug2/graph.json')
