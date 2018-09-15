#!/usr/bin/env python3

import json
import io
import csv
from model.network import Network
from model.drug import Drug
from model.disease import Disease
from model.edge import Edge

network = Network()
with io.open('../data/DrugCentral/drugcentral_mappings.csv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter=',', quotechar='"')
    next(reader, None)
    for row in reader:
        ids = ['DrugCentral:%s' % row[0], 'DrugBank:%s' % row[1]]
        if row[2] is not None and len(row[2]) > 0:
            ids.append('RxNorm:%s' % row[2])
        network.add_node(Drug(ids, [row[3]]))

with io.open('../data/DrugCentral/drugcentral_indications.csv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter=',', quotechar='"')
    next(reader, None)
    for row in reader:
        disease = Disease(['SnomedCT:%s' % row[2], 'UMLS:%s' % row[3]], [row[1]])
        network.add_node(disease)
        e = Edge('DrugBank:%s' % row[0], 'SnomedCT:%s' % row[2], 'INDICATES', {'source': 'DrugCentral'})
        network.add_edge(e)

with io.open('../data/DrugCentral/drugcentral_contraindications.csv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter=',', quotechar='"')
    next(reader, None)
    for row in reader:
        network.add_node(Disease(['SnomedCT:%s' % row[2], 'UMLS:%s' % row[3]], [row[1]]))
        e = Edge('DrugBank:%s' % row[0], 'SnomedCT:%s' % row[2], 'CONTRAINDICATES', {'source': 'DrugCentral'})
        network.add_edge(e)

with io.open('../data/DrugCentral/graph.json', 'w', encoding='utf-8', newline='') as f:
    f.write(json.dumps(network.to_dict(), indent=2))