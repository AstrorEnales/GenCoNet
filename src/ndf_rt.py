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

file = '../data/NDF-RT/NDFRT_Public_All.xml'
zip_file = '../data/NDF-RT/NDFRT_Public_All 2018-02-05.zip'
url = 'https://evs.nci.nih.gov/ftp1/NDF-RT/Archive/NDFRT_Public_All%202018-02-05.zip'

if not os.path.exists(file):
    print('Database does not exist. Trying to download and extract...')
    if not os.path.exists(zip_file):
        print('Downloading latest archive...')
        with urllib.request.urlopen(url) as response, open(zip_file, 'wb') as f:
            f.write(response.read())
    print('Extracting database file...')
    with zipfile.ZipFile(zip_file) as z:
        pattern = re.compile(r'NDFRT_Public_[0-9.]+/NDFRT_Public_[0-9.]+_TDE\.xml')
        date_filename = [x for x in z.namelist() if pattern.fullmatch(x)]
        if len(date_filename) != 1:
            print('Failed to find database file in archive. Extract the file manually as "NDFRT_Public_All.xml".')
            exit(1)
        else:
            with open(file, 'wb') as f:
                f.write(z.read(date_filename[0]))

root = xml.etree.ElementTree.parse(file).getroot()

kind_defs = {}
kind_defs_rev = {}
for kind in root.findall('kindDef'):
    kind_defs[kind.find('code').text] = kind.find('name').text
    kind_defs_rev[kind.find('name').text] = kind.find('code').text

role_defs = {}
for role in root.findall('roleDef'):
    role_defs[role.find('code').text] = role.find('name').text

property_defs = {}
property_defs_rev = {}
for prop in root.findall('propertyDef'):
    property_defs[prop.find('code').text] = prop.find('name').text
    property_defs_rev[prop.find('name').text] = prop.find('code').text

concept_defs = {}
for concept in root.findall('conceptDef'):
    concept_code = concept.find('code').text
    concept_defs[concept_code] = {
        'name': concept.find('name').text,
        'code': concept.find('code').text,
        'kind': concept.find('kind').text,
        'roles': [],
        'properties': []
    }
    roles = concept.find('definingRoles')
    for role in roles.findall('role'):
        concept_defs[concept_code]['roles'].append([role.find('name').text, role.find('value').text])
    properties = concept.find('properties')
    for prop in properties.findall('property'):
        concept_defs[concept_code]['properties'].append([prop.find('name').text, prop.find('value').text])

network = Network()
for concept_code in sorted(concept_defs.keys()):
    concept = concept_defs[concept_code]
    if concept['kind'] == kind_defs_rev['DRUG_KIND']:
        drug_ids = ['NDF-RT:%s' % concept['code']]
        drug_names = [concept['name']]
        for prop in concept['properties']:
            if property_defs[prop[0]] == 'RxNorm_CUI':
                drug_ids.append('RxNorm:%s' % prop[1])
            elif property_defs[prop[0]] == 'RxNorm_Name':
                drug_names.append(prop[1])
            elif property_defs[prop[0]] == 'UMLS_CUI':
                drug_ids.append('UMLS:%s' % prop[1])
            elif property_defs[prop[0]] == 'Synonym':
                drug_names.append(prop[1])
        drug_names = [x.replace('[VA Product]', '').strip() for x in drug_names]
        drug = Drug(drug_ids, drug_names)
        network.add_node(drug)
        for role in concept['roles']:
            role_name = role_defs[role[0]]
            rel = {'source': 'NDF-RT'}
            if role_name == 'induces {NDFRT}':
                network.add_edge(Edge(next(iter(drug.ids)), 'NDF-RT:%s' % role[1], 'INDUCES', rel))
            elif role_name == 'CI_with {NDFRT}':
                network.add_edge(Edge(next(iter(drug.ids)), 'NDF-RT:%s' % role[1], 'CONTRAINDICATES', rel))
            elif role_name == 'may_treat {NDFRT}':
                network.add_edge(Edge(next(iter(drug.ids)), 'NDF-RT:%s' % role[1], 'INDICATES', rel))
    elif concept['kind'] == kind_defs_rev['DISEASE_KIND']:
        disease_ids = ['NDF-RT:%s' % concept['code']]
        disease_names = [concept['name']]
        for prop in concept['properties']:
            if property_defs[prop[0]] == 'SNOMED_CID':
                disease_ids.append('SnomedCT:%s' % prop[1])
            elif property_defs[prop[0]] == 'UMLS_CUI':
                disease_ids.append('UMLS:%s' % prop[1])
            elif property_defs[prop[0]] == 'MeSH_CUI':
                disease_ids.append('MESH:%s' % prop[1])
            elif property_defs[prop[0]] == 'Synonym':
                disease_names.append(prop[1])
        disease_names = [x.replace('[Disease/Finding]', '').strip() for x in disease_names]
        disease = Disease(disease_ids, disease_names)
        network.add_node(disease)

with io.open('../data/NDF-RT/graph.json', 'w', encoding='utf-8', newline='') as f:
    f.write(json.dumps(network.to_dict(), indent=2))
