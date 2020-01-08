#!/usr/bin/env python3

import os
import urllib.request
import gzip
import io
import csv
import xml.etree.ElementTree

from model.edge import Edge
from model.network import Network
from model.gene import Gene
from model.go_class import GOClass

annotations_file = '../data/GO/goa_human.gaf'
annotations_zip_file = '../data/GO/goa_human.gaf.gz'
annotations_url = 'http://geneontology.org/gene-associations/goa_human.gaf.gz'
ontology_file = '../data/GO/go.owl'
ontology_url = 'http://purl.obolibrary.org/obo/go.owl'

if not os.path.exists(annotations_file):
    print('Database does not exist. Trying to download and extract...')
    if not os.path.exists(annotations_zip_file):
        print('Downloading latest archive...')
        with urllib.request.urlopen(annotations_url) as response, open(annotations_zip_file, 'wb') as f:
            f.write(response.read())
    print('Extracting database file...')
    with gzip.open(annotations_zip_file, 'rb') as f:
        with open(annotations_file, 'wb') as out_file:
            out_file.write(f.read())

if not os.path.exists(ontology_file):
    print('Ontology does not exist. Trying to download...')
    with urllib.request.urlopen(ontology_url) as response, open(ontology_file, 'wb') as f:
        f.write(response.read())

network = Network()

go_class_redirects = {}
go_class_ns_lookup = {}
owl_ns = '{http://www.w3.org/2002/07/owl#}'
obo_in_owl_ns = '{http://www.geneontology.org/formats/oboInOwl#}'
rdfs_ns = '{http://www.w3.org/2000/01/rdf-schema#}'
root = xml.etree.ElementTree.parse(ontology_file).getroot()
for owl_class in root.findall(owl_ns + 'Class'):
    id_node = owl_class.find(obo_in_owl_ns + 'id')
    obo_ns_node = owl_class.find(obo_in_owl_ns + 'hasOBONamespace')
    label_node = owl_class.find(rdfs_ns + 'label')
    if id_node is not None and obo_ns_node is not None:
        go_class = GOClass([id_node.text], [label_node.text])
        network.add_node(go_class)
        go_class_ns_lookup[id_node.text] = obo_ns_node.text
        for alternative_id_node in owl_class.findall(obo_in_owl_ns + 'hasAlternativeId'):
            go_class_redirects[alternative_id_node.text] = id_node.text

with io.open(annotations_file, 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    for row in reader:
        if not row[0][0].startswith('!') and row[12] == 'taxon:9606':
            gene = Gene(['UniProtKB:%s' % row[1], 'HGNC:%s' % row[2]], [])
            network.add_node(gene)
            if row[4] not in go_class_ns_lookup:
                # print('[WARN] GO id %s is obsolete, redirecting to %s' % (row[4], go_class_redirects[row[4]]))
                row[4] = go_class_redirects[row[4]]
            label = go_class_ns_lookup[row[4]].upper()
            if label == 'MOLECULAR_FUNCTION':
                label = 'HAS_' + label
            elif label == 'BIOLOGICAL_PROCESS':
                label = 'BELONGS_TO_' + label
            elif label == 'CELLULAR_COMPONENT':
                label = 'IN_' + label
            e = Edge(gene, network.get_node_by_id(row[4], 'GOClass'), label, {'source': 'GO,%s' % row[5]})
            network.add_edge(e)

network.save('../data/GO/graph.json')
