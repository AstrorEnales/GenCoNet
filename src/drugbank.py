#!/usr/bin/env python3

import io
import csv
import os.path
import zipfile
import lxml.etree as etree

from model.adr import AdverseDrugReaction
from model.network import Network
from model.drug import Drug
from model.gene import Gene
from model.edge import Edge
from model.variant import Variant


def get_drugbank_id(drug_node, ns):
    ids = drug_node.findall(ns + 'drugbank-id')
    if ids:
        primary_ids = [x for x in ids if 'primary' in x.attrib and x.attrib['primary'] == 'true']
        return primary_ids[0].text if primary_ids else ids[0].text
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

identifiers_output_file = '../data/DrugBank/drug_identifiers.csv'
targets_output_file = '../data/DrugBank/drugs_target_human_genes.csv'
interactions_output_file = '../data/DrugBank/drug_interactions.csv'
snp_adr_output_file = '../data/DrugBank/drug_snp_adrs.csv'
if any([not os.path.exists(x) for x in
        [targets_output_file, interactions_output_file, identifiers_output_file, snp_adr_output_file]]):
    positive = {'inducer', 'agonist', 'activator', 'partial agonist', 'stimulator', 'positive modulator',
                'positive allosteric modulator'}
    negative = {'blocker', 'antagonist', 'antibody', 'weak inhibitor', 'suppressor', 'partial antagonist',
                'negative modulator', 'inverse agonist', 'inhibitor', 'inactivator', 'antisense oligonucleotide',
                'inhibitory allosteric modulator', 'inhibitory immune response'}

    drugs = {}
    ns = '{http://www.drugbank.ca}'
    for event, drug_node in etree.iterparse(file, events=('end',), tag=ns + 'drug'):
        # Skip non high level drug nodes
        if 'created' not in drug_node.attrib:
            continue
        drugbank_id = get_drugbank_id(drug_node, ns)
        if drugbank_id is None:
            continue
        drug_name = drug_node.find(ns + 'name').text
        # Collect all external identifiers
        external_ids = set()
        external_identifiers_node = drug_node.find(ns + 'external-identifiers')
        if external_identifiers_node is not None:
            for external_identifier_node in external_identifiers_node:
                # ChEMBL, Wikipedia, UniProtKB, PharmGKB, KEGG Drug, PubChem Compound
                source = external_identifier_node.find(ns + 'resource').text.strip()
                identifier = external_identifier_node.find(ns + 'identifier').text.strip()
                if source == 'ChEMBL':
                    external_ids.add('ChEMBL:%s' % identifier)
                elif source == 'KEGG Drug':
                    external_ids.add('Kegg:%s' % identifier)
                elif source == 'PubChem Compound':
                    external_ids.add('PubChem:CID%s' % identifier)
        # Collect all gene targets for the drug
        targets = []
        targets_node = drug_node.find(ns + 'targets')
        if targets_node is not None:
            for target_node in targets_node:
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
        # Collect all interactions for the drug
        interactions = []
        interactions_node = drug_node.find(ns + 'drug-interactions')
        if interactions_node is not None:
            for interaction_node in interactions_node:
                id2 = interaction_node.find(ns + 'drugbank-id').text
                description = interaction_node.find(ns + 'description').text
                interactions.append([id2, description])
        # Collect all SNP ADRs for the drug
        snp_adrs = []
        snp_adrs_node = drug_node.find(ns + 'snp-adverse-drug-reactions')
        if snp_adrs_node is not None:
            for snp_adr_node in snp_adrs_node:
                gene_symbol = snp_adr_node.find(ns + 'gene-symbol').text
                rs_id = snp_adr_node.find(ns + 'rs-id').text
                adverse_reaction = snp_adr_node.find(ns + 'adverse-reaction').text
                description = snp_adr_node.find(ns + 'description').text
                pubmed_id = snp_adr_node.find(ns + 'pubmed-id').text
                snp_adrs.append([gene_symbol, rs_id, adverse_reaction, description, pubmed_id])
        drugs[drugbank_id] = [drug_name, targets, interactions, external_ids, snp_adrs]
        drug_node.clear()

    targets_results = []
    interactions_results = []
    external_id_results = []
    snp_adrs_results = []
    for drugbank_id in sorted(drugs.keys()):
        drug = drugs[drugbank_id]
        for target in drug[1]:
            if target['taxon_id'] == '9606':
                targets_results.append([drugbank_id, drug[0], target['gene'], target['gene_name'], target['hgnc_id'],
                                        target['known_action'], ','.join(target['actions']),
                                        target['simplified_action']])
        for interaction in drug[2]:
            if interaction[0] in drugs:
                interactions_results.append(
                    [drugbank_id, drug[0], interaction[0], drugs[interaction[0]][0], interaction[1]])

        chembl_id = None
        kegg_id = None
        pubchem_id = None
        for external_id in drug[3]:
            if external_id.startswith('ChEMBL:'):
                chembl_id = external_id
            elif external_id.startswith('Kegg:'):
                kegg_id = external_id
            elif external_id.startswith('PubChem:'):
                pubchem_id = external_id
        external_id_results.append([drugbank_id, chembl_id, kegg_id, pubchem_id])

        for snp_adr in drug[4]:
            snp_adrs_results.append([drugbank_id] + snp_adr)

    with io.open(targets_output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow(['drugbank_id', 'drug_name', 'gene', 'gene_name', 'hgnc_id', 'known_action', 'actions',
                         'simplified_action'])
        for row in targets_results:
            writer.writerow(row)
    with io.open(interactions_output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow(['drugbank_id1', 'drug_name1', 'drugbank_id2', 'drug_name2', 'description'])
        for row in interactions_results:
            writer.writerow(row)
    with io.open(identifiers_output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow(['drugbank_id', 'ChEMBL', 'KEGG Drug', 'PubChem Compound'])
        for row in external_id_results:
            writer.writerow(row)
    with io.open(snp_adr_output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow(['drugbank_id', 'gene_symbol', 'rs_id', 'adverse_reaction', 'description', 'pubmed_id'])
        for row in snp_adrs_results:
            writer.writerow(row)
else:
    targets_results = []
    with io.open(targets_output_file, 'r', encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter=',', quotechar='"')
        next(reader, None)
        for row in reader:
            targets_results.append(row)
    interactions_results = []
    with io.open(interactions_output_file, 'r', encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter=',', quotechar='"')
        next(reader, None)
        for row in reader:
            interactions_results.append(row)
    external_id_results = []
    with io.open(identifiers_output_file, 'r', encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter=',', quotechar='"')
        next(reader, None)
        for row in reader:
            external_id_results.append(row)
    snp_adrs_results = []
    with io.open(snp_adr_output_file, 'r', encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter=',', quotechar='"')
        next(reader, None)
        for row in reader:
            snp_adrs_results.append(row)

external_id_lookup = {}
for row in external_id_results:
    external_id_lookup[row[0]] = [x for x in row[1::] if x]

network = Network()
for row in targets_results:
    drug_ids = ['DrugBank:%s' % row[0]]
    if row[0] in external_id_lookup:
        drug_ids.extend(external_id_lookup[row[0]])
    drug = Drug(drug_ids, [row[1]])
    network.add_node(drug)
    gene_ids = ['HGNC:%s' % row[2]]
    if row[4]:
        gene_ids.append(row[4])
    gene = Gene(gene_ids, [row[3]])
    network.add_node(gene)
    rel = {
        'source': 'DrugBank',
        'known_action': row[5] == 1,
        'actions': row[6].split(',') if row[6] else [],
        'simplified_action': row[7]
    }
    network.add_edge(Edge(drug, gene, 'TARGETS', rel))
for row in interactions_results:
    drug1 = Drug(['DrugBank:%s' % row[0]], [row[1]])
    network.add_node(drug1)
    drug2 = Drug(['DrugBank:%s' % row[2]], [row[3]])
    network.add_node(drug2)
    rel = {
        'source': 'DrugBank',
        'description': row[4]
    }
    network.add_edge(Edge(drug1, drug2, 'INTERACTS', rel))
adr_id_counter = 1
for row in snp_adrs_results:
    # drugbank_id, gene_symbol, rs_id, adverse_reaction, description, pubmed_id
    adr = AdverseDrugReaction(['GenCoNet:DrugBank_ADR_%s' % adr_id_counter], [])
    adr_id_counter += 1
    adr.attributes = {
        'source': 'DrugBank',
        'adverse_reaction': row[3],
        'description': row[4],
        'pmid': row[5]
    }
    network.add_node(adr)
    drug = Drug(['DrugBank:%s' % row[0]], [])
    network.add_node(drug)
    network.add_edge(Edge(drug, adr, 'HAS_ADR', {'source': 'DrugBank'}))
    if row[2] and len(row[2]) > 0:
        variant = Variant(['dbSNP:%s' % row[2]], [])
        network.add_node(variant)
        network.add_edge(Edge(variant, adr, 'ASSOCIATED_WITH_ADR', {'source': 'DrugBank'}))
    gene = Gene(['HGNC:%s' % row[1]], [])
    network.add_node(gene)
    network.add_edge(Edge(gene, adr, 'ASSOCIATED_WITH_ADR', {'source': 'DrugBank'}))

network.save('../data/DrugBank/graph.json')
