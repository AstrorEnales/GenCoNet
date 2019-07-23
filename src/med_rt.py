#!/usr/bin/env python3

import json
import os.path
import xml.etree.ElementTree
import io
import urllib.request
import zipfile
import re
from model.network import Network
from model.drug import Drug
from model.disease import Disease
from model.edge import Edge

file = '../data/MED-RT/Core_MEDRT_XML.xml'
zip_file = '../data/MED-RT/Core_MEDRT_XML.zip'
url = 'https://evs.nci.nih.gov/ftp1/MED-RT/Core_MEDRT_XML.zip'
if not os.path.exists(file):
    print('Database does not exist. Trying to download and extract...')
    if not os.path.exists(zip_file):
        print('Downloading latest archive...')
        with urllib.request.urlopen(url) as response, open(zip_file, 'wb') as f:
            f.write(response.read())
    print('Extracting database file...')
    with zipfile.ZipFile(zip_file) as z:
        pattern = re.compile(r'Core_MEDRT_[0-9.]+_XML\.xml')
        date_filename = [x for x in z.namelist() if pattern.fullmatch(x)]
        if len(date_filename) != 1:
            print('Failed to find database file in archive. Extract the file manually as "Core_MEDRT_XML.xml".')
            exit(1)
        else:
            with open(file, 'wb') as f:
                f.write(z.read(date_filename[0]))

network = Network()

root = xml.etree.ElementTree.parse(file).getroot()
added_rxnorm_drugs = set()
added_mesh_diseases = set()
for association in root.findall('association'):
    association_type = association.find('name').text
    from_namespace = association.find('from_namespace').text
    from_id = association.find('from_code').text
    from_name = association.find('from_name').text
    to_namespace = association.find('to_namespace').text
    to_id = association.find('to_code').text
    to_name = association.find('to_name').text
    if from_namespace != 'RxNorm' or to_namespace != 'MeSH':
        continue
    if association_type not in ['induces', 'CI_with', 'may_treat']:
        continue
    drug_id = 'RxNorm:%s' % from_id
    if from_id not in added_rxnorm_drugs:
        drug = Drug([drug_id], [from_name])
        network.add_node(drug)
        added_rxnorm_drugs.add(from_id)
    disease_id = 'MeSH:%s' % to_id
    if to_id not in added_mesh_diseases:
        disease = Disease([disease_id], [to_name])
        network.add_node(disease)
        added_mesh_diseases.add(to_id)
    rel = {'source': 'MEDRT'}
    if association_type == 'induces':
        network.add_edge(Edge(drug_id, disease_id, 'INDUCES', rel))
    elif association_type == 'CI_with':
        network.add_edge(Edge(drug_id, disease_id, 'CONTRAINDICATES', rel))
    elif association_type == 'may_treat':
        network.add_edge(Edge(drug_id, disease_id, 'INDICATES', rel))

with io.open('../data/MED-RT/graph.json', 'w', encoding='utf-8', newline='') as f:
    f.write(json.dumps(network.to_dict(), indent=2))
