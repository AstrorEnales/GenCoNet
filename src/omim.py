#!/usr/bin/env python3

import io
import csv
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
        if row[8] and row[12]:
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

with io.open('../data/OMIM/filtered_associations.csv', 'w', encoding='utf-8', newline='') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"')
    writer.writerow(['Location', 'Phenotype', 'Phenotype MIM number', 'Inheritance', 'Phenotype mapping key',
                     'Gene/Locus', 'Gene/Locus MIM number'])
    for row in filtered_results:
        writer.writerow(row)

network = Network()

# 0 Location
# 1 Phenotype
# 2 Phenotype MIM number
# 3 Inheritance
# 4 Phenotype mapping key
# 5 Gene/Locus
# 6 Gene/Locus MIM number
with io.open('../data/OMIM/filtered_associations.csv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter=',', quotechar='"')
    next(reader, None)
    for row in reader:
        disease = Disease(['OMIM:%s' % row[2]], [])
        network.add_node(disease)
        gene = Gene(['HGNC:%s' % row[5]], [])  # , 'OMIM:%s' % row[6]
        network.add_node(gene)
        rel = {
            'source': 'OMIM',
            'location': row[0],
            'phenotype': row[1],
            'inheritance': row[2],
            'phenotype_mapping_key': row[4]
        }
        network.add_edge(Edge(gene, disease, 'ASSOCIATES_WITH', rel))

network.save('../data/OMIM/graph.json')
