#!/usr/bin/env python3

import os.path
import urllib.request
import io
import json


def map_to_simple_id(s: str) -> str or None:
    if s.startswith('http://linkedlifedata.com/resource/umls/id/'):
        return 'UMLS:%s' % s.split('umls/id/')[1]
    if s.startswith('http://identifiers.org/omim/'):
        return 'OMIM:%s' % s.split('omim/')[1]
    if s.startswith('http://identifiers.org/snomedct/'):
        return 'SnomedCT:%s' % s.split('snomedct/')[1]
    if s.startswith('http://identifiers.org/mesh/'):
        return 'MESH:%s' % s.split('mesh/')[1]
    if s.startswith('http://purl.obolibrary.org/obo/DOID_'):
        return 'DO:%s' % s.split('DOID_')[1]
    if s.startswith('http://identifiers.org/meddra/'):
        return 'MEDDRA:%s' % s.split('meddra/')[1]
    # http://identifiers.org/medgen/
    # http://www.orpha.net/ORDO/Orphanet_
    # http://purl.obolibrary.org/obo/NCIT_
    return None


file = '../data/MONDO/mondo.json'
url = 'http://purl.obolibrary.org/obo/mondo.json'

if not os.path.exists(file):
    print('Database does not exist. Trying to download...')
    with urllib.request.urlopen(url) as response, open(file, 'wb') as f:
        f.write(response.read())

lookup = {}
reverse_lookup = {}

with io.open(file, 'r', encoding='utf-8', newline='') as f:
    root = json.load(f)['graphs']
    for graph in root:
        '''
        nodes
        edges
        id
        meta
        equivalentNodesSets
        logicalDefinitionAxioms
        domainRangeAxioms
        propertyChainAxioms
        '''
        if graph['id'] != 'http://purl.obolibrary.org/obo/mondo.owl':
            continue
        for node in graph['nodes']:
            if 'MONDO_' not in node['id']:
                continue
            node_id = 'MONDO:%s' % node['id'].split('MONDO_')[1]
            node_label = node['lbl'] if 'lbl' in node else None
            n = {'label': node_label, 'refs': []}
            lookup[node_id] = n
            if 'meta' not in node or 'basicPropertyValues' not in node['meta']:
                continue
            for reference in node['meta']['basicPropertyValues']:
                if reference['pred'] == 'http://www.w3.org/2004/02/skos/core#exactMatch':
                    simple_id = map_to_simple_id(reference['val'])
                    if simple_id is not None:
                        n['refs'].append(simple_id)
                        if simple_id not in reverse_lookup:
                            reverse_lookup[simple_id] = []
                        reverse_lookup[simple_id].append(node_id)

with io.open('../data/MONDO/lookup.json', 'w', encoding='utf-8', newline='') as f:
    f.write(json.dumps(lookup, indent=2))

with io.open('../data/MONDO/reverse_lookup.json', 'w', encoding='utf-8', newline='') as f:
    f.write(json.dumps(reverse_lookup, indent=2))
