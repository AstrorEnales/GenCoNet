#!/usr/bin/env python3

import os.path
import xml.etree.ElementTree
import io
import csv
import urllib.request
import zipfile
import re

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
for kind in root.findall('kindDef'):
    kind_defs[kind.find('code').text] = kind.find('name').text

role_defs = {}
for role in root.findall('roleDef'):
    role_defs[role.find('code').text] = role.find('name').text

concept_defs = {}
for concept in root.findall('conceptDef'):
    concept_code = concept.find('code').text
    concept_defs[concept_code] = {
        'name': concept.find('name').text,
        'code': concept.find('code').text,
        'kind': concept.find('kind').text,
        'roles': []
    }
    roles = concept.find('definingRoles')
    for role in roles.findall('role'):
        concept_defs[concept_code]['roles'].append([role.find('name').text, role.find('value').text])

with io.open('../data/NDF-RT/concepts.csv', 'w', encoding='utf-8', newline='') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"')
    writer.writerow(['concept_code', 'concept_kind', 'concept_name', 'role_code', 'role_name', 'role_value_code',
                     'role_value_name'])
    for concept_code in sorted(concept_defs.keys()):
        concept = concept_defs[concept_code]
        for role in concept['roles']:
            role_name = role_defs[role[0]]
            role_value_name = concept_defs[role[1]]['name']
            writer.writerow([concept_code, concept['kind'], concept['name'], role[0], role_name, role[1],
                             role_value_name])
