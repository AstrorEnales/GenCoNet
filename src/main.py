#!/usr/bin/env python3

import csv
import io
import os
import re
import shutil

from exporter.graphml_exporter import GraphMLExporter
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
            ids = ['DrugCentral:%s' % row[0], 'DrugBank:%s' % row[1]]
            if row[2] is not None and len(row[2]) > 0:
                ids.append('RxNorm:%s' % row[2])
            drug = Drug(ids, [row[3]])
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
                    'num_pmids': int(row[5]),
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
                    'num_pmids': int(row[4]),
                    'score': row[5]
                }
                network.variant_associates_with_disease.append([next(iter(variant.ids)), next(iter(disease.ids)), rel])


def import_gwas_catalog(network: Network):
    # 0 DATE ADDED TO CATALOG
    # 1 PUBMEDID
    # 2 FIRST AUTHOR
    # 3 DATE
    # 4 JOURNAL
    # 5 LINK
    # 6 STUDY
    # 7 DISEASE/TRAIT
    # 8 INITIAL SAMPLE SIZE
    # 9 REPLICATION SAMPLE SIZE
    # 10 REGION
    # 11 CHR_ID
    # 12 CHR_POS
    # 13 REPORTED GENE(S)
    # 14 MAPPED_GENE
    # 15 UPSTREAM_GENE_ID
    # 16 DOWNSTREAM_GENE_ID
    # 17 SNP_GENE_IDS
    # 18 UPSTREAM_GENE_DISTANCE
    # 19 DOWNSTREAM_GENE_DISTANCE
    # 20 STRONGEST SNP-RISK ALLELE
    # 21 SNPS
    # 22 MERGED
    # 23 SNP_ID_CURRENT
    # 24 CONTEXT
    # 25 INTERGENIC
    # 26 RISK ALLELE FREQUENCY
    # 27 P-VALUE
    # 28 PVALUE_MLOG
    # 29 P-VALUE (TEXT)
    # 30 OR or BETA
    # 31 95% CI (TEXT)
    # 32 PLATFORM [SNPS PASSING QC]
    # 33 CNV
    loc_pattern = re.compile(r'LOC[0-9]+')
    with io.open('../data/GWAS-Catalog/gwas_catalog_associations.tsv', 'r', encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter='\t', quotechar='"')
        next(reader, None)
        for row in reader:
            if row[14] is None or len(row[14]) == 0:
                continue
            gene_ids = row[14].replace(' x ', ', ').replace(' - ', ', ').split(', ')
            for gene_id in gene_ids:
                if loc_pattern.fullmatch(gene_id) is not None:
                    continue
                gene = Gene(['HGNCSymbol:%s' % gene_id], [])
                network.add_gene(gene)
                variant = Variant(['dbSNP:%s' % row[21]], [])
                network.add_variant(variant)
                network.gene_codes_variant.append(
                    [gene.get_id(), variant.get_id(), {'source': 'GWASCatalog', 'pmid': row[1]}])


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
        writer.writerow([':START_ID(Gene-ID)', 'source:string', 'num_pmids:int', 'num_snps:int', 'score:string',
                         ':END_ID(Disease-ID)', ':TYPE'])
        for row in network.gene_associates_with_disease:
            writer.writerow([network.get_gene_by_id(row[0]).get_id(), row[2]['source'], row[2]['num_pmids'],
                             row[2]['num_snps'], row[2]['score'], network.get_disease_by_id(row[1]).get_id(),
                             'ASSOCIATES_WITH'])
    # Save variant -[ASSOCIATES_WITH]-> disease relationships
    with io.open(os.path.join(output_path, 'variant_associates_with_disease.csv'), 'w', encoding='utf-8', newline='') \
            as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Variant-ID)', 'source:string', 'num_pmids:int', 'score:string',
                         ':END_ID(Disease-ID)', ':TYPE'])
        for row in network.variant_associates_with_disease:
            writer.writerow([network.get_variant_by_id(row[0]).get_id(), row[2]['source'], row[2]['num_pmids'],
                             row[2]['score'], network.get_disease_by_id(row[1]).get_id(), 'ASSOCIATES_WITH'])
    # Save gene -[CODES]-> variant relationships
    with io.open(os.path.join(output_path, 'gene_codes_variant.csv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow([':START_ID(Gene-ID)', 'source:string', 'pmid:int', ':END_ID(Variant-ID)', ':TYPE'])
        for row in network.gene_codes_variant:
            writer.writerow([network.get_gene_by_id(row[0]).get_id(), row[2]['source'], row[2]['pmid'],
                             network.get_variant_by_id(row[1]).get_id(), 'CODES'])
    with io.open(os.path.join(output_path, 'create_indices.cypher'), 'w', encoding='utf-8', newline='') as f:
        f.write('create constraint on (p:Drug) assert p._id is unique;\n')
        f.write('create constraint on (p:Gene) assert p._id is unique;\n')
        f.write('create constraint on (p:Variant) assert p._id is unique;\n')
        f.write('create constraint on (p:Disease) assert p._id is unique;\n')
    with io.open(os.path.join(output_path, 'import_admin.bat'), 'w', encoding='utf-8', newline='') as f:
        f.write('@echo off\n')
        f.write('E:/runtime/neo4j-community-3.3.1/bin/neo4j-admin import --nodes diseases.csv --nodes drugs.csv ' +
                '--nodes genes.csv --nodes variants.csv --relationships drug_contraindicates_disease.csv ' +
                '--relationships drug_indicates_disease.csv --relationships drug_targets_gene.csv ' +
                '--relationships gene_associates_with_disease.csv --relationships gene_codes_variant.csv ' +
                '--relationships variant_associates_with_disease.csv\n')


if __name__ == '__main__':
    network = Network()
    # Import
    import_drugcentral(network)
    import_drugbank(network)
    import_disgenet(network)
    import_gwas_catalog(network)
    # Cleanup
    network.prune()
    # Export
    save_network(network, '../output/')

    # GraphMLExporter(network).save('../output/graph.graphml')
