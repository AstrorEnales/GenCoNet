#!/usr/bin/env python3

import csv
import io
import os
import json
import shutil

import mondo_mapper

from utils import name_utils

from model.disease import Disease
from model.network import Network


def cleanup_output(output_path: str):
    # Cleanup previous export
    if os.path.exists(output_path) and os.path.isdir(output_path):
        shutil.rmtree(output_path)
    os.mkdir(output_path)


def merge_duplicate_node_names(network: Network):
    for node in network.nodes.values():
        node.names = name_utils.normalize_drug_names(node.names)


def save_network(network: Network, config):
    output_path = config['output-path']
    # Save nodes
    with io.open(os.path.join(output_path, 'nodes.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow(['_id:ID(Node-ID)', 'ids:string[]', 'names:string[]', ':LABEL'])
        for n in set(network.nodes.values()):
            writer.writerow([n.id, ';'.join(n.ids), ';'.join(n.names), n.label])

    # Save INDICATES relationships
    with io.open(os.path.join(output_path, 'rel_INDICATES.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Node-ID)', 'source:string', ':END_ID(Node-ID)', ':TYPE'])
        for e in network.get_edges_by_label('INDICATES'):
            writer.writerow([network.get_node_by_id(e.source).id, e.attributes['source'],
                             network.get_node_by_id(e.target).id, e.label])

    # Save CONTRAINDICATES relationships
    with io.open(os.path.join(output_path, 'rel_CONTRAINDICATES.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Node-ID)', 'source:string', ':END_ID(Node-ID)', ':TYPE'])
        for e in network.get_edges_by_label('CONTRAINDICATES'):
            writer.writerow([network.get_node_by_id(e.source).id, e.attributes['source'],
                             network.get_node_by_id(e.target).id, e.label])

    # Save TARGETS relationships
    with io.open(os.path.join(output_path, 'rel_TARGETS.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Node-ID)', 'source:string', 'known_action:boolean', 'actions:string[]',
                         'simplified_action:string', ':END_ID(Node-ID)', ':TYPE'])
        for e in network.get_edges_by_label('TARGETS'):
            writer.writerow([network.get_node_by_id(e.source).id, e.attributes['source'],
                             'true' if e.attributes['known_action'] else 'false',
                             ';'.join(e.attributes['actions']), e.attributes['simplified_action'],
                             network.get_node_by_id(e.target).id, e.label])

    # Save ASSOCIATES_WITH relationships
    with io.open(os.path.join(output_path, 'rel_ASSOCIATES_WITH.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Node-ID)', 'source:string', 'num_pmids:int', 'num_snps:int', 'score:string',
                         ':END_ID(Node-ID)', ':TYPE'])
        for e in network.get_edges_by_label('ASSOCIATES_WITH'):
            num_pmids = e.attributes['num_pmids'] if 'num_pmids' in e.attributes else None
            num_snps = e.attributes['num_snps'] if 'num_snps' in e.attributes else None
            score = e.attributes['score'] if 'score' in e.attributes else None
            writer.writerow([network.get_node_by_id(e.source).id, e.attributes['source'], num_pmids, num_snps, score,
                             network.get_node_by_id(e.target).id, e.label])

    # Save INDUCES relationships
    with io.open(os.path.join(output_path, 'rel_INDUCES.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Node-ID)', 'source:string', ':END_ID(Node-ID)', ':TYPE'])
        for e in network.get_edges_by_label('INDUCES'):
            writer.writerow([network.get_node_by_id(e.source).id, e.attributes['source'],
                             network.get_node_by_id(e.target).id, e.label])

    # Save CODES relationships
    with io.open(os.path.join(output_path, 'rel_CODES.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Node-ID)', 'source:string', 'pmid:int', ':END_ID(Node-ID)', ':TYPE'])
        for e in network.get_edges_by_label('CODES'):
            writer.writerow([network.get_node_by_id(e.source).id, e.attributes['source'], e.attributes['pmid'],
                             network.get_node_by_id(e.target).id, e.label])

    # Save INTERACTS relationships
    with io.open(os.path.join(output_path, 'rel_INTERACTS.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Node-ID)', 'source:string', 'description:string', ':END_ID(Node-ID)', ':TYPE'])
        for e in network.get_edges_by_label('INTERACTS'):
            writer.writerow([network.get_node_by_id(e.source).id, e.attributes['source'], e.attributes['description'],
                             network.get_node_by_id(e.target).id, e.label])

    with io.open(os.path.join(output_path, 'create_indices.cypher'), 'w', encoding='utf-8', newline='') as f:
        f.write('create constraint on (p:Drug) assert p._id is unique;\n')
        f.write('create constraint on (p:Gene) assert p._id is unique;\n')
        f.write('create constraint on (p:Variant) assert p._id is unique;\n')
        f.write('create constraint on (p:Disease) assert p._id is unique;\n')
    with io.open(os.path.join(output_path, 'import_admin.bat'), 'w', encoding='utf-8', newline='') as f:
        f.write('@echo off\n')
        f.write('net stop neo4j\n')
        f.write('rmdir /s "%s"\n' % os.path.join(config['Neo4j']['database-path'], config['Neo4j']['database-name']))
        f.write(os.path.join(config['Neo4j']['bin-path'], 'neo4j-admin'))
        f.write(' import ' +
                '--database %s ' % config['Neo4j']['database-name'] +
                '--nodes nodes.csv ' +
                '--relationships rel_INDICATES.csv ' +
                '--relationships rel_CONTRAINDICATES.csv ' +
                '--relationships rel_TARGETS.csv ' +
                '--relationships rel_ASSOCIATES_WITH.csv ' +
                '--relationships rel_INDUCES.csv ' +
                '--relationships rel_CODES.csv ' +
                '--relationships rel_INTERACTS.csv > import.log\n')
        f.write('net start neo4j\n')
        f.write(os.path.join(config['Neo4j']['bin-path'], 'cypher-shell'))
        f.write(' -u neo4j -p root --non-interactive < create_indices.cypher 1>> import.log 2>&1\n')
    with io.open(os.path.join(output_path, 'import_admin.sh'), 'w', encoding='utf-8', newline='') as f:
        f.write(os.path.join(config['Neo4j']['bin-path'], 'neo4j-admin'))
        f.write(' import ' +
                '--database %s ' % config['Neo4j']['database-name'] +
                '--nodes nodes.csv ' +
                '--relationships rel_INDICATES.csv ' +
                '--relationships rel_CONTRAINDICATES.csv ' +
                '--relationships rel_TARGETS.csv ' +
                '--relationships rel_ASSOCIATES_WITH.csv ' +
                '--relationships rel_CODES.csv ' +
                '--relationships rel_INTERACTS.csv > import.log\n')


if __name__ == '__main__':
    with io.open('../data/config.json', 'r', encoding='utf-8', newline='') as f:
        config = json.load(f)

    network = Network()
    # Import
    graphs = [
        '../data/DisGeNet/graph.json',
        '../data/DrugBank/graph.json',
        '../data/DrugCentral/graph.json',
        '../data/GWAS-Catalog/graph.json',
        '../data/HGNC/graph.json',
        '../data/HPO/graph.json',
        '../data/MED-RT/graph.json',
        '../data/NDF-RT/graph.json',
        '../data/OMIM/graph.json',
        '../data/HuGE-Navigator/graph.json',
        '../data/SIDER/graph.json',
        # '../data/SuperDrug2/graph.json',
        # '../data/PubMed/graph.json',
    ]
    # Fusion
    for graph in graphs:
        with io.open(graph, 'r', encoding='utf-8', newline='') as f:
            g = json.loads(f.read())
            network.load_from_dict(g)
    # Mapping
    all_disease_ids = set()
    for node in network.get_nodes_by_label('Disease'):
        all_disease_ids.update(node.ids)
    for disease_id in all_disease_ids:
        mapped_ids, mapped_names = mondo_mapper.map_from(disease_id)
        if len(mapped_ids) > 0:
            network.add_node(Disease(mapped_ids, mapped_names))
    # Cleanup
    network.prune()
    merge_duplicate_node_names(network)
    # Export
    cleanup_output(config['output-path'])
    with io.open(os.path.join(config['output-path'], 'graph.json'), 'w', encoding='utf-8', newline='') as f:
        f.write(json.dumps(network.to_dict(), separators=(',', ':')))
    save_network(network, config)

'''
MATCH p=(disease2:Disease)-[*1..2]-(:Gene)-[:TARGETS]-(drug:Drug)-[:INDICATES]-(disease1:Disease)
    WHERE disease1<>disease2
    AND ANY(name IN disease1.names WHERE toLower(name) =~ ".*hypertension.*")
    AND ANY(name IN disease2.names WHERE toLower(name) =~ ".*asthma.*")
    RETURN drug._id, drug.names[0], disease1.names[0], collect(distinct disease2.names[0])
    LIMIT 100
'''
