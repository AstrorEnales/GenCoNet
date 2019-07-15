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
        for f in os.listdir(output_path):
            file_path = os.path.join(output_path, f)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(e)
    else:
        os.mkdir(output_path)


def merge_duplicate_node_names(network: Network):
    for node in network.nodes.values():
        node.names = name_utils.normalize_node_names(node.names)


def save_network(network: Network, config):
    output_path = config['output-path']
    # Save nodes
    with io.open(os.path.join(output_path, 'nodes.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow(['_id:ID(Node-ID)', 'ids:string[]', 'names:string[]', ':LABEL'])
        for n in set(network.nodes.values()):
            writer.writerow([n.id, ';'.join(n.ids), ';'.join(n.names), n.label])

    # Save HAS_MOLECULAR_FUNCTION relationships
    with io.open(os.path.join(output_path, 'rel_HAS_MOLECULAR_FUNCTION.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Node-ID)', 'source:string', ':END_ID(Node-ID)', ':TYPE'])
        for e in network.get_edges_by_label('HAS_MOLECULAR_FUNCTION'):
            writer.writerow([network.get_node_by_id(e.source).id, e.attributes['source'],
                             network.get_node_by_id(e.target).id, e.label])

    # Save BELONGS_TO_BIOLOGICAL_PROCESS relationships
    with io.open(os.path.join(output_path, 'rel_BELONGS_TO_BIOLOGICAL_PROCESS.csv'), 'w', encoding='utf-8',
                 newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Node-ID)', 'source:string', ':END_ID(Node-ID)', ':TYPE'])
        for e in network.get_edges_by_label('BELONGS_TO_BIOLOGICAL_PROCESS'):
            writer.writerow([network.get_node_by_id(e.source).id, e.attributes['source'],
                             network.get_node_by_id(e.target).id, e.label])

    # Save IN_CELLULAR_COMPONENT relationships
    with io.open(os.path.join(output_path, 'rel_IN_CELLULAR_COMPONENT.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Node-ID)', 'source:string', ':END_ID(Node-ID)', ':TYPE'])
        for e in network.get_edges_by_label('IN_CELLULAR_COMPONENT'):
            writer.writerow([network.get_node_by_id(e.source).id, e.attributes['source'],
                             network.get_node_by_id(e.target).id, e.label])

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
            known_action = ('true' if e.attributes['known_action'] else 'false') \
                if 'known_action' in e.attributes else None
            simplified_action = e.attributes['simplified_action'] if 'simplified_action' in e.attributes else None
            writer.writerow([network.get_node_by_id(e.source).id, e.attributes['source'],
                             known_action, ';'.join(e.attributes['actions']), simplified_action,
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

    # Save EQTL relationships
    with io.open(os.path.join(output_path, 'rel_EQTL.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Node-ID)', 'source:string', 'pvalue:string', 'snp_chr:string', 'cis_trans:string',
                         ':END_ID(Node-ID)', ':TYPE'])
        for e in network.get_edges_by_label('EQTL'):
            writer.writerow([network.get_node_by_id(e.source).id, e.attributes['source'], e.attributes['pvalue'],
                             e.attributes['snp_chr'], e.attributes['cis_trans'], network.get_node_by_id(e.target).id,
                             e.label])

    # Save INTERACTS relationships
    with io.open(os.path.join(output_path, 'rel_INTERACTS.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Node-ID)', 'source:string', 'description:string', ':END_ID(Node-ID)', ':TYPE'])
        for e in network.get_edges_by_label('INTERACTS'):
            writer.writerow([network.get_node_by_id(e.source).id, e.attributes['source'], e.attributes['description'],
                             network.get_node_by_id(e.target).id, e.label])

    with io.open(os.path.join(output_path, 'create_indices.cypher'), 'w', encoding='utf-8', newline='') as f:
        for node_label in network.node_labels():
            f.write('create constraint on (p:%s) assert p._id is unique;\n' % node_label)
    with io.open(os.path.join(output_path, 'import_admin.bat'), 'w', encoding='utf-8', newline='') as f:
        f.write('@echo off\n')
        f.write('net stop neo4j\n')
        f.write('rmdir /s "%s"\n' % os.path.join(config['Neo4j']['database-path'], config['Neo4j']['database-name']))
        f.write('CALL ' + os.path.join(config['Neo4j']['bin-path'], 'neo4j-admin'))
        f.write(' import ' +
                '--database %s ' % config['Neo4j']['database-name'] +
                '--nodes nodes.csv ' +
                ' '.join(['--relationships rel_%s.csv' % x for x in network.edge_labels()]) +
                ' > import.log\n')
        f.write('net start neo4j\n')
        f.write(os.path.join(config['Neo4j']['bin-path'], 'cypher-shell'))
        f.write(' -u %s -p %s --non-interactive < create_indices.cypher 1>> import.log 2>&1\n'
                % (config['Neo4j']['user'], config['Neo4j']['password']))
    with io.open(os.path.join(output_path, 'import_admin.sh'), 'w', encoding='utf-8', newline='') as f:
        f.write(os.path.join(config['Neo4j']['bin-path'], 'neo4j-admin'))
        f.write(' import ' +
                '--database %s ' % config['Neo4j']['database-name'] +
                '--nodes nodes.csv ' +
                ' '.join(['--relationships rel_%s.csv' % x for x in network.edge_labels()]) +
                ' > import.log\n')


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
        '../data/DGIdb/graph.json',
        '../data/Westra_etal_2017/graph.json',
        '../data/SuperDrug2/graph.json',
        '../data/UniprotKB/graph.json',
        '../data/GO/graph.json',
        # '../data/PubMed/graph.json',
    ]
    # Fusion
    print('[INFO] Network fusion')
    for graph in graphs:
        print('[INFO] Add network', graph)
        with io.open(graph, 'r', encoding='utf-8', newline='') as f:
            g = json.loads(f.read())
            network.load_from_dict(g)
    # Mapping
    print('[INFO] Add disease mappings')
    all_disease_ids = set()
    for node in network.get_nodes_by_label('Disease'):
        all_disease_ids.update(node.ids)
    for disease_id in all_disease_ids:
        mapped_ids, mapped_names = mondo_mapper.map_from(disease_id)
        if len(mapped_ids) > 0:
            network.add_node(Disease(mapped_ids, mapped_names))
    # Cleanup
    print('[INFO] Prune network')
    network.prune()
    print('[INFO] Merge duplicate node names')
    merge_duplicate_node_names(network)
    print('[INFO] Merge duplicate edges')
    network.merge_duplicate_edges()
    # Export
    print('[INFO] Export network')
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
