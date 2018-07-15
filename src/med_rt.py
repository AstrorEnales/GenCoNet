#!/usr/bin/env python3

import os.path
import xml.etree.ElementTree
import io
import csv
import urllib.request
import zipfile
import re

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

root = xml.etree.ElementTree.parse(file).getroot()
for t in [('MED_RT_contraindications', 'CI_with'), ('MED_RT_indications', 'may_treat'),
          ('MED_RT_inductions', 'induces')]:
    with io.open('../data/MED-RT/%s.csv' % t[0], 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow(['rxnorm_id', 'mesh_id', 'mesh_term'])
        for association in root.findall('association'):
            association_type = association.find('name').text
            from_namespace = association.find('from_namespace').text
            from_id = association.find('from_code').text
            to_namespace = association.find('to_namespace').text
            to_id = association.find('to_code').text
            to_name = association.find('to_name').text
            if from_namespace != 'RxNorm' or to_namespace != 'MeSH':
                continue
            if association_type == t[1]:
                writer.writerow([from_id, to_id, to_name])
