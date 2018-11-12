#!/usr/bin/env python3

import json
import io

from model.edge import Edge
from model.network import Network
from model.gene import Gene
from model.disease import Disease

file = '../data/HuGE-Navigator/Disease-GeneID.txt'

network = Network()

with io.open(file, 'r', encoding='utf-8', newline='') as f:
    for skip in range(0, 4):
        f.readline()
    for line in f:
        parts = [[y.strip() for y in x[:-1].split('(')] for x in line.strip().split('\t')]
        if len(parts) > 1:
            disease = Disease(['UMLS:%s' % parts[0][1]], [parts[0][0]])
            network.add_node(disease)
            for part in parts[1::]:
                gene = Gene(['HGNCSymbol:%s' % part[0], 'Entrez:%s' % part[1]], [part[0]])
                network.add_node(gene)
                rel = {'source': 'HuGE Navigator'}
                network.add_edge(Edge(next(iter(gene.ids)), next(iter(disease.ids)), 'ASSOCIATES_WITH', rel))

with io.open('../data/HuGE-Navigator/graph.json', 'w', encoding='utf-8', newline='') as f:
    f.write(json.dumps(network.to_dict(), indent=2))
