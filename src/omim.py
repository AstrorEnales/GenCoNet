#!/usr/bin/env python3

import io
import csv
import json
import requests
from model.network import Network
from model.gene import Gene
from model.disease import Disease
from model.edge import Edge

filtered_results = []

with io.open('../data/OMIM/genemap2.txt', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    for row in reader:
        if row[0].startswith('#'):
            continue
        if row[8] is not None and len(row[8]) > 0 and row[12] is not None and len(row[12]) > 0:
            phenotypes = [x.strip() for x in row[12].split(';')]
            for phenotype in phenotypes:
                last_left_brace = phenotype.rfind('(')
                last_right_brace = phenotype.rfind(')')
                if last_right_brace != last_left_brace + 2:
                    continue
                mapping_key = int(phenotype[last_left_brace + 1])
                inheritance = phenotype[last_right_brace + 1::].strip(', ')
                phenotype = phenotype[0:last_left_brace].strip()
                if ',' not in phenotype:
                    phenotype_mim = None
                else:
                    last_comma = phenotype.rfind(',')
                    mim_part = phenotype[last_comma + 1::].strip()
                    if mim_part.isdecimal():
                        phenotype_mim = int(mim_part)
                        phenotype = phenotype[0:last_comma]
                    else:
                        phenotype_mim = None
                # 3 - Location
                # Phenotype
                # Phenotype MIM number
                # Inheritance
                # Phenotype mapping key
                # 8 - Gene/Locus
                # 5 - Gene/Locus MIM number
                if phenotype_mim is not None:
                    filtered_results.append([row[3], phenotype, phenotype_mim, inheritance, mapping_key,
                                             row[8], row[5]])

mappings = {}
for row in filtered_results:
    mappings['OMIM:%s' % row[2]] = [None, None]
sorted_keys = sorted(mappings.keys())
chunk_size = 20
i = 0
while i < len(sorted_keys):
    chunk = sorted_keys[i:i + chunk_size]
    print('[', i, ',', i + chunk_size, ']', '/', len(sorted_keys))
    r = requests.post('https://www.ebi.ac.uk/spot/oxo/api/search',
                      json={"ids": chunk, "mappingTarget": ["UMLS"], "distance": "1"})
    if r.status_code == 500:
        if chunk_size == 1:
            i += 1
        else:
            chunk_size = 1
    else:
        query_results = r.json()['_embedded']['searchResults']
        for row in query_results:
            if len(row['mappingResponseList']) > 0:
                target = row['mappingResponseList'][0]
                mappings[row['curie']] = [target['curie'], target['label']]
        i += chunk_size
        chunk_size = 20

with io.open('../data/OMIM/omim_to_umls.csv', 'w', encoding='utf-8', newline='') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"')
    writer.writerow(['OMIM', 'UMLS', 'UMLS name'])
    for key in sorted(mappings.keys()):
        writer.writerow([key] + mappings[key])

with io.open('../data/OMIM/filtered_associations.csv', 'w', encoding='utf-8', newline='') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"')
    writer.writerow(['Location', 'Phenotype', 'Phenotype MIM number', 'Inheritance', 'Phenotype mapping key',
                     'Gene/Locus', 'Gene/Locus MIM number', 'Phenotype UMLS', 'Phenotype UMLS name'])
    for row in filtered_results:
        writer.writerow(row + mappings['OMIM:%s' % row[2]])

network = Network()

# 0 Location
# 1 Phenotype
# 2 Phenotype MIM number
# 3 Inheritance
# 4 Phenotype mapping key
# 5 Gene/Locus
# 6 Gene/Locus MIM number
# 7 Phenotype UMLS
# 8 Phenotype UMLS name
with io.open('../data/OMIM/filtered_associations.csv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter=',', quotechar='"')
    next(reader, None)
    for row in reader:
        disease_ids = ['OMIM:%s' % row[2]]
        if row[7] is not None and len(row[7]) > 0:
            disease_ids.append(row[7])
        else:
            # TODO: currently can't use a disease without UMLS id
            continue
        disease_names = [] if row[8] is not None and len(row[8]) > 0 else []
        disease = Disease(disease_ids, disease_names)
        network.add_node(disease)
        gene = Gene(['HGNCSymbol:%s' % row[5], 'OMIM:%s' % row[6]], [])
        network.add_node(gene)
        rel = {
            'source': 'OMIM',
            'location': row[0],
            'phenotype': row[1],
            'inheritance': row[2],
            'phenotype_mapping_key': row[4]
        }
        network.add_edge(Edge(next(iter(gene.ids)), next(iter(disease.ids)), 'ASSOCIATES_WITH', rel))

with io.open('../data/OMIM/graph.json', 'w', encoding='utf-8', newline='') as f:
    f.write(json.dumps(network.to_dict(), indent=2))
