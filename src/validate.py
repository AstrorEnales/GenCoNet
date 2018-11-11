#!/usr/bin/env python3

import os
import io
import json

from utils import name_utils

from model.network import Network

if __name__ == '__main__':
    with io.open('../data/config.json', 'r', encoding='utf-8', newline='') as f:
        config = json.load(f)

    network = Network()
    with io.open(os.path.join(config['output-path'], 'graph.json'), 'r', encoding='utf-8', newline='') as f:
        g = json.loads(f.read())
        network.load_from_dict(g)

    drug_count = 0
    drug_check_failed_count = 0
    for drug in network.get_nodes_by_label('Drug'):
        drug_count += 1
        indications = {x.target for x in network.get_node_edges_by_label(drug, 'INDICATES')}
        contraindications = {x.target for x in network.get_node_edges_by_label(drug, 'CONTRAINDICATES')}
        if not indications.isdisjoint(contraindications):
            drug_check_failed_count += 1
            print('[WARN] Plausibility check failed: Drug has same disease as indication and contraindication.')
            print('\tDrug:', drug)
            for intersection in indications.intersection(contraindications):
                disease = network.get_node_by_id(intersection)
                print('\tIntersection:', disease)
                for indication in network.get_edges_from_to(drug, disease, 'INDICATES'):
                    print('\t\t', indication.attributes['source'], indication)
                for contraindication in network.get_edges_from_to(drug, disease, 'CONTRAINDICATES'):
                    print('\t\t', contraindication.attributes['source'], contraindication)

    drug_multiple_names_count = 0
    for drug in network.get_nodes_by_label('Drug'):
        if not name_utils.drug_names_synonym(drug.names):
            drug_multiple_names_count += 1
            print('Drug:', drug.ids)
            for name in drug.names:
                print('\t', name)

    print('%s/%s drugs failed check' % (drug_check_failed_count, drug_count))
    print('%s/%s drugs multiple names' % (drug_multiple_names_count, drug_count))

    for node in network.get_nodes():
        id_prefix = {}
        for node_id in node.ids:
            prefix, suffix = node_id.split(':')
            if prefix not in id_prefix:
                id_prefix[prefix] = 0
            id_prefix[prefix] += 1
        if any([id_prefix[prefix] > 1 for prefix in id_prefix]):
            print('[WARN] Possible loss of distinction for node:', node)
