#!/usr/bin/env python3

import csv
import io
import os
import shutil

from model.disease import Disease
from model.drug import Drug
from model.gene import Gene
from model.network import Network
from model.variant import Variant


def import_drugcentral(network: Network):
    with io.open('../data/DrugCentral/drugcentral_mappings.csv', 'r', encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter=',', quotechar='"')
        next(reader, None)
        for row in reader:
            drug = Drug(['DrugCentral:%s' % row[0], 'DrugBank:%s' % row[1], 'RxNorm:%s' % row[2]], [row[3]])
            network.add_drug(drug)

    with io.open('../data/DrugCentral/drugcentral_indications.csv', 'r', encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter=',', quotechar='"')
        next(reader, None)
        for row in reader:
            disease = Disease(['SnomedCT:%s' % row[2], 'UMLS:%s' % row[3]], [row[1]])
            network.add_disease(disease)
            network.drug_indicates_disease.append(
                ['DrugBank:%s' % row[0], 'SnomedCT:%s' % row[2], {'source': 'DrugCentral'}])

    with io.open('../data/DrugCentral/drugcentral_contraindications.csv', 'r', encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter=',', quotechar='"')
        next(reader, None)
        for row in reader:
            disease = Disease(['SnomedCT:%s' % row[2], 'UMLS:%s' % row[3]], [row[1]])
            network.add_disease(disease)
            network.drug_contraindicates_disease.append(
                ['DrugBank:%s' % row[0], 'SnomedCT:%s' % row[2], {'source': 'DrugCentral'}])


def import_drugbank(network: Network):
    with io.open('../data/DrugBank/drugs_target_human_genes.csv', 'r', encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter=',', quotechar='"')
        next(reader, None)
        for row in reader:
            drug = Drug(['DrugBank:%s' % row[0]], [row[1]])
            network.add_drug(drug)
            gene = Gene([row[4], 'HGNCSymbol:%s' % row[2]], [row[3]])
            network.add_gene(gene)
            rel = {'source': 'DrugBank', 'known_action': row[5] == 1, 'actions': row[6].split(','),
                   'simplified_action': row[7]}
            network.drug_targets_gene.append([next(iter(drug.ids)), next(iter(gene.ids)), rel])


def import_disgenet(network: Network):
    PUBMED_COUNT_THRESHOLD = 2

    with io.open('../data/DisGeNet/curated_gene_disease_associations.tsv', 'r', encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter='\t', quotechar='"')
        next(reader, None)
        for row in reader:
            # geneId
            # geneSymbol
            # diseaseId
            # diseaseName
            # score
            # NofPmids
            # NofSnps
            # source
            if int(row[5]) >= PUBMED_COUNT_THRESHOLD:
                gene = Gene(['HGNCSymbol:%s' % row[1]], [])
                network.add_gene(gene)
                disease = Disease(['UMLS:%s' % row[2]], [row[3]])
                network.add_disease(disease)
                rel = {
                    'source': 'DisGeNet,%s' % row[7],
                    'num_pubmed_ids': int(row[5]),
                    'num_snps': int(row[6]),
                    'score': row[4]
                }
                network.gene_associates_with_disease.append([next(iter(gene.ids)), next(iter(disease.ids)), rel])

    with io.open('../data/DisGeNet/curated_variant_disease_associations.tsv', 'r', encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter='\t', quotechar='"')
        next(reader, None)
        for row in reader:
            # snpId
            # diseaseId
            # diseaseName
            # score
            # NofPmids
            # source
            if int(row[4]) >= PUBMED_COUNT_THRESHOLD:
                variant = Variant(['dbSNP:%s' % row[0]], [])
                network.add_variant(variant)
                disease = Disease(['UMLS:%s' % row[1]], [row[2]])
                network.add_disease(disease)
                rel = {
                    'source': 'DisGeNet,%s' % row[5],
                    'num_pubmed_ids': int(row[4]),
                    'score': row[5]
                }
                network.variant_associates_with_disease.append([next(iter(variant.ids)), next(iter(disease.ids)), rel])


def save_network(network: Network, output_path: str):
    # Cleanup previous export
    if os.path.exists(output_path) and os.path.isdir(output_path):
        shutil.rmtree(output_path)
    os.mkdir(output_path)
    # Save drug nodes
    with io.open(os.path.join(output_path, 'drugs.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow(['_id:ID(Drug-ID)', 'ids:string[]', 'names:string[]', ':LABEL'])
        for drug in set(network.drugs.values()):
            writer.writerow([drug.get_id(), ';'.join(drug.ids), ';'.join(drug.names), 'Drug'])
    # Save disease nodes
    with io.open(os.path.join(output_path, 'diseases.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow(['_id:ID(Disease-ID)', 'ids:string[]', 'names:string[]', ':LABEL'])
        for disease in set(network.diseases.values()):
            writer.writerow([disease.get_id(), ';'.join(disease.ids), ';'.join(disease.names), 'Disease'])
    # Save gene nodes
    with io.open(os.path.join(output_path, 'genes.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow(['_id:ID(Gene-ID)', 'ids:string[]', 'names:string[]', ':LABEL'])
        for gene in set(network.genes.values()):
            writer.writerow([gene.get_id(), ';'.join(gene.ids), ';'.join(gene.names), 'Gene'])
    # Save variant nodes
    with io.open(os.path.join(output_path, 'variants.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow(['_id:ID(Variant-ID)', 'ids:string[]', 'names:string[]', ':LABEL'])
        for variant in set(network.variants.values()):
            writer.writerow([variant.get_id(), ';'.join(variant.ids), ';'.join(variant.names), 'Variant'])
    # Save drug -[INDICATES]-> disease relationships
    with io.open(os.path.join(output_path, 'drug_indicates_disease.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Drug-ID)', 'source:string', ':END_ID(Disease-ID)', ':TYPE'])
        for row in network.drug_indicates_disease:
            writer.writerow([network.get_drug_by_id(row[0]).get_id(), row[2]['source'],
                             network.get_disease_by_id(row[1]).get_id(), 'INDICATES'])
    # Save drug -[CONTRAINDICATES]-> disease relationships
    with io.open(os.path.join(output_path, 'drug_contraindicates_disease.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Drug-ID)', 'source:string', ':END_ID(Disease-ID)', ':TYPE'])
        for row in network.drug_contraindicates_disease:
            writer.writerow([network.get_drug_by_id(row[0]).get_id(), row[2]['source'],
                             network.get_disease_by_id(row[1]).get_id(), 'CONTRAINDICATES'])
    # Save drug -[TARGETS]-> gene relationships
    with io.open(os.path.join(output_path, 'drug_targets_gene.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Drug-ID)', 'source:string', 'known_action:boolean', 'actions:string[]',
                         'simplified_action:string', ':END_ID(Gene-ID)', ':TYPE'])
        for row in network.drug_targets_gene:
            writer.writerow([network.get_drug_by_id(row[0]).get_id(), row[2]['source'],
                             'true' if row[2]['known_action'] else 'false', ';'.join(row[2]['actions']),
                             row[2]['simplified_action'], network.get_gene_by_id(row[1]).get_id(), 'TARGETS'])
    # Save gene -[ASSOCIATES_WITH]-> disease relationships
    with io.open(os.path.join(output_path, 'gene_associates_with_disease.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Gene-ID)', 'source:string', 'num_pubmed_ids:int', 'num_snps:int', 'score:string',
                         ':END_ID(Disease-ID)', ':TYPE'])
        for row in network.gene_associates_with_disease:
            writer.writerow([network.get_gene_by_id(row[0]).get_id(), row[2]['source'], row[2]['num_pubmed_ids'],
                             row[2]['num_snps'], row[2]['score'], network.get_disease_by_id(row[1]).get_id(),
                             'ASSOCIATES_WITH'])
    # Save variant -[ASSOCIATES_WITH]-> disease relationships
    with io.open(os.path.join(output_path, 'variant_associates_with_disease.csv'), 'w', encoding='utf-8', newline='') \
            as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Variant-ID)', 'source:string', 'num_pubmed_ids:int', 'score:string',
                         ':END_ID(Disease-ID)', ':TYPE'])
        for row in network.variant_associates_with_disease:
            writer.writerow([network.get_variant_by_id(row[0]).get_id(), row[2]['source'], row[2]['num_pubmed_ids'],
                             row[2]['score'], network.get_disease_by_id(row[1]).get_id(), 'ASSOCIATES_WITH'])


if __name__ == '__main__':
    network = Network()
    # Import
    import_drugcentral(network)
    import_drugbank(network)
    import_disgenet(network)
    # Cleanup
    # TODO: network.prune()
    # Export
    save_network(network, '../output/')
