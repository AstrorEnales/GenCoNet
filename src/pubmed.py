#!/usr/bin/env python3

import io
import csv
from model.network import Network
from model.drug import Drug
from model.disease import Disease
from model.edge import Edge

network = Network()

with io.open('../data/PubMed/drug_disease.csv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter=',', quotechar='"')
    next(reader, None)
    for row in reader:
        drug = Drug(['DrugBank:%s' % row[1]], [row[0]])
        disease = Disease([row[4]], [row[3]])
        network.add_node(drug)
        network.add_node(disease)
        network.add_edge(Edge(drug, disease, row[2], {'source': 'PubMed', 'pmid': row[5]}))

network.save('../data/PubMed/graph.json')
