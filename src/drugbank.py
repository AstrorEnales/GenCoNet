#!/usr/bin/env python3

import json
import os.path
import xml.etree.ElementTree
import io
import zipfile
import csv
from model.network import Network
from model.drug import Drug
from model.gene import Gene
from model.edge import Edge


def get_drugbank_id(drug_node, ns):
    ids = drug_node.findall(ns + 'drugbank-id')
    if len(ids) > 0:
        primary_ids = [x for x in ids if 'primary' in x.attrib and x.attrib['primary'] == 'true']
        return primary_ids[0].text if len(primary_ids) > 0 else ids[0].text
    return None


file = '../data/DrugBank/full database.xml'
zip_file = '../data/DrugBank/drugbank_all_full_database.xml.zip'
if not os.path.exists(file):
    print('Extracting database file...')
    with zipfile.ZipFile(zip_file) as z:
        date_filename = [x for x in z.namelist() if x == 'full database.xml']
        if len(date_filename) != 1:
            print('Failed to find database file in archive. Extract the file manually as "full database.xml".')
            exit(1)
        else:
            with open(file, 'wb') as f:
                f.write(z.read(date_filename[0]))

output_file = '../data/DrugBank/drugs_target_human_genes.csv'
if not os.path.exists(output_file):
    positive = ['inducer', 'agonist', 'activator', 'partial agonist', 'stimulator', 'positive modulator',
                'positive allosteric modulator']
    negative = ['blocker', 'antagonist', 'antibody', 'weak inhibitor', 'suppressor', 'partial antagonist',
                'negative modulator', 'inverse agonist', 'inhibitor', 'inactivator', 'antisense oligonucleotide',
                'inhibitory allosteric modulator', 'inhibitory immune response']

    drugs = {}
    ns = '{http://www.drugbank.ca}'
    root = xml.etree.ElementTree.parse(file).getroot()
    for drug_node in root.findall(ns + 'drug'):
        drugbank_id = get_drugbank_id(drug_node, ns)
        if drugbank_id is None:
            continue
        drug_name = drug_node.find(ns + 'name').text
        # Collect all gene targets for the drug
        targets = []
        targets_node = drug_node.find(ns + 'targets')
        if targets_node is not None:
            for target_node in targets_node.findall(ns + 'target'):
                actions = set()
                actions_node = target_node.find(ns + 'actions')
                if actions_node is not None:
                    for action_node in actions_node.findall(ns + 'action'):
                        actions.add(action_node.text)
                known_action_node = target_node.find(ns + 'known-action')
                known_action = None
                if known_action_node is not None:
                    known_action = 1 if known_action_node.text == 'yes' else 0
                for polypeptide_node in target_node.findall(ns + 'polypeptide'):
                    hgnc_id = None
                    external_ids_node = polypeptide_node.find(ns + 'external-identifiers')
                    for external_id_node in external_ids_node.findall(ns + 'external-identifier'):
                        if external_id_node.find(ns + 'resource').text == 'HUGO Gene Nomenclature Committee (HGNC)':
                            hgnc_id = external_id_node.find(ns + 'identifier').text
                            break
                    simplified_action = sum({1 if x in positive else (-1 if x in negative else 0) for x in actions})
                    simplified_action = 'positive' if simplified_action == 1 else (
                        'negative' if simplified_action == -1 else 'other')
                    target = {
                        'gene': polypeptide_node.find(ns + 'gene-name').text,
                        'gene_name': polypeptide_node.find(ns + 'name').text,
                        'hgnc_id': hgnc_id,
                        'taxon_id': polypeptide_node.find(ns + 'organism').attrib['ncbi-taxonomy-id'],
                        'known_action': known_action,
                        'simplified_action': simplified_action,
                        'actions': sorted(actions)
                    }
                    targets.append(target)
        if len(targets) > 0:
            drugs[drugbank_id] = [drug_name, targets]

    results = []
    for drugbank_id in sorted(drugs.keys()):
        drug = drugs[drugbank_id]
        for target in drug[1]:
            if target['taxon_id'] == '9606':
                results.append([drugbank_id, drug[0], target['gene'], target['gene_name'], target['hgnc_id'],
                                target['known_action'], ','.join(target['actions']), target['simplified_action']])

    with io.open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow(['drugbank_id', 'drug_name', 'gene', 'gene_name', 'hgnc_id', 'known_action', 'actions',
                         'simplified_action'])
        for row in results:
            writer.writerow(row)
else:
    results = []
    with io.open(output_file, 'r', encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter=',', quotechar='"')
        next(reader, None)
        for row in reader:
            results.append(row)

network = Network()
for row in results:
    drug = Drug(['DrugBank:%s' % row[0]], [row[1]])
    network.add_node(drug)
    gene_ids = ['HGNCSymbol:%s' % row[2]]
    if row[4] is not None and len(row[4]) > 0:
        gene_ids.append(row[4])
    gene = Gene(gene_ids, [row[3]])
    network.add_node(gene)
    rel = {
        'source': 'DrugBank',
        'known_action': row[5] == 1,
        'actions': row[6].split(',') if row[6] is not None and len(row[6]) > 0 else [],
        'simplified_action': row[7]
    }
    network.add_edge(Edge(next(iter(drug.ids)), next(iter(gene.ids)), 'TARGETS', rel))

with io.open('../data/DrugBank/graph.json', 'w', encoding='utf-8', newline='') as f:
    f.write(json.dumps(network.to_dict(), indent=2))
