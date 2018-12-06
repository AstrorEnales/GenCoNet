#!/usr/bin/env python3

import os
import io
import json

from utils import name_utils

from model.network import Network


def node_ids_to_links(node_ids):
    result = ''
    for node_id in node_ids:
        prefix, suffix = node_id.split(':')
        result += (', ' if len(result) > 0 else '')
        if prefix == 'DrugBank':
            result += '<a href="https://www.drugbank.ca/%s" target="_blank">%s</a>' % (suffix, node_id)
        elif prefix == 'DrugCentral':
            result += '<a href="http://drugcentral.org/drugcard/%s" target="_blank">%s</a>' % (suffix, node_id)
        elif prefix == 'MESH':
            result += '<a href="http://identifiers.org/mesh/%s" target="_blank">%s</a>' % (suffix, node_id)
        elif prefix == 'UMLS':
            result += '<a href="http://linkedlifedata.com/resource/umls/id/%s" target="_blank">%s</a>' % (suffix, node_id)
        else:
            result += node_id
    return result


if __name__ == '__main__':
    with io.open('../data/config.json', 'r', encoding='utf-8', newline='') as f:
        config = json.load(f)

    network = Network()
    with io.open(os.path.join(config['output-path'], 'graph.json'), 'r', encoding='utf-8', newline='') as f:
        g = json.loads(f.read())
        network.load_from_dict(g)

    with io.open(os.path.join(config['output-path'], 'validation_report.html'), 'w', encoding='utf-8', newline='') as f:
        f.write('<!doctype html>\n<html>\n<head>\n</head>\n<body>\n')
        f.write('<ul>\n')
        f.write('<li><a href="#indi-contra-section">Drugs indicating and contraindicating same disease</a></li>\n')
        f.write('<li><a href="#synonym-section">Drugs with synonyms</a></li>\n')
        f.write('<li><a href="#loss-distinction-section">Possible loss of distinction</a></li>\n')
        f.write('<li><a href="#single-id-nodes">Nodes with single IDs</a></li>\n')
        f.write('</ul>\n')
        f.write('<h1><a name="indi-contra-section">Drugs indicating and contraindicating same disease:</a></h1>\n')
        f.write('<table border="1">\n<thead>\n<tr><th>Drug</th><th>Indications</th><th>Contraindications</th></tr>\n</thead>\n')
        f.write('<tbody>\n')
        drug_count = 0
        drug_check_failed_count = 0
        for drug in network.get_nodes_by_label('Drug'):
            drug_count += 1
            indications = {x.target for x in network.get_node_edges_by_label(drug, 'INDICATES')}
            contraindications = {x.target for x in network.get_node_edges_by_label(drug, 'CONTRAINDICATES')}
            if not indications.isdisjoint(contraindications):
                drug_check_failed_count += 1

                for intersection in indications.intersection(contraindications):
                    disease = network.get_node_by_id(intersection)
                    drug_text = '%s<br/>%s' % (node_ids_to_links(drug.ids), '<br/>'.join(drug.names))
                    indications_text = '<br/>'.join(['%s: %s -> %s' % (x.attributes['source'], node_ids_to_links([x.source]), node_ids_to_links([x.target])) for x in network.get_edges_from_to(drug, disease, 'INDICATES')])
                    contraindications_text = '<br/>'.join(['%s: %s -> %s' % (x.attributes['source'], node_ids_to_links([x.source]), node_ids_to_links([x.target])) for x in network.get_edges_from_to(drug, disease, 'CONTRAINDICATES')])
                    f.write('<tr><td>%s</td><td style="white-space: nowrap;">%s</td><td style="white-space: nowrap;">%s</td></tr>\n' % (drug_text, indications_text, contraindications_text))
        f.write('</tbody>\n</table>\n')

        f.write('<h1><a name="synonym-section">Drugs with synonyms</a></h1>\n')
        f.write('<table border="1">\n<thead>\n<tr><th>Drug</th><th>Names</th></tr>\n</thead>\n')
        f.write('<tbody>\n')
        drug_multiple_names_count = 0
        for drug in network.get_nodes_by_label('Drug'):
            if not name_utils.node_names_synonym(drug.names):
                drug_multiple_names_count += 1
                f.write('<tr><td>%s</td><td>%s</td></tr>\n' % (node_ids_to_links(drug.ids), '<br/>'.join(drug.names)))
        f.write('</tbody>\n</table>\n')

        f.write('<h1><a name="loss-distinction-section">Possible loss of distinction:</a></h1>\n')
        f.write('<table border="1">\n<thead>\n<tr><th>Node type</th><th>Node IDs</th><th>Names</th></tr>\n</thead>\n')
        f.write('<tbody>\n')
        for node in network.get_nodes():
            id_prefix = {}
            for node_id in node.ids:
                prefix, suffix = node_id.split(':')
                if prefix not in id_prefix:
                    id_prefix[prefix] = 0
                id_prefix[prefix] += 1
            if any([id_prefix[prefix] > 1 for prefix in id_prefix]):
                f.write('<tr><td>%s</td><td>%s</td><td>%s</td></tr>\n' % (node.label, node_ids_to_links(sorted(node.ids)), '<br/>'.join(node.names)))
        f.write('</tbody>\n</table>\n')

        f.write('<h1><a name="single-id-nodes">Nodes with single IDs:</a></h1>\n')
        f.write('<table border="1">\n<thead>\n<tr><th>Node type</th><th>Node IDs</th><th>Names</th></tr>\n</thead>\n')
        f.write('<tbody>\n')
        for node in network.get_nodes():
            if len(node.ids) <= 1:
                f.write('<tr><td>%s</td><td>%s</td><td>%s</td></tr>\n' % (node.label, node_ids_to_links(sorted(node.ids)), '<br/>'.join(node.names)))
        f.write('</tbody>\n</table>\n')

        f.write('</body>\n</html>\n')
