#!/usr/bin/env python3

import os
import io
import sys
import csv
import json
import zipfile
import urllib.request
from typing import List, Set, Tuple

from model.gene import Gene
from model.network import Network
from model.drug import Drug
from model.disease import Disease
from model.edge import Edge
from model.variant import Variant

maxInt = sys.maxsize
while True:
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt / 10)


def split_list(value: str) -> List[str]:
    if value is None or len(value) <= 0:
        return []
    return value.replace('","', ';').replace(',"', ';').strip('"').split(';')


def process_drug_cross_references(ids: List[str]) -> Set[str]:
    # TODO consider more ids:
    # BindingDB:10015
    # ChEBI:CHEBI:100147
    # Chemical Abstracts Service:100427-26-7
    # ChemSpider:10080
    # ClinicalTrials.gov:NCT00005520
    # Drugs Product Database (DPD):00000027
    # FDA Drug Label at DailyMed:0013824b-6aee-4da4-affd-35bc6bf19d91
    # GenBank:AF280603
    # HET:097
    # HMDB:HMDB00042
    # IUPHAR Ligand:100
    # KEGG Compound:C00002
    # KEGG Drug:D00006
    # National Drug Code Directory:0002-0604-40
    # PDB:097
    # Therapeutic Targets Database:DAP000001
    # UniProtKB:P00451
    # URL:http://en.wikipedia.org/wiki/Abacavir

    # Not used:
    # PubChem Substance:10299466

    # Used ids:
    # DrugBank:DB00001
    # PubChem Compound:10090750
    filtered_ids = set()
    for drug_id in ids:
        if drug_id.startswith('DrugBank:'):
            filtered_ids.add(drug_id)
        elif drug_id.startswith('PubChem Compound:'):
            filtered_ids.add('PubChem:CID%s' % drug_id.split(':')[1])
    return filtered_ids


def process_gene_cross_references(ids: List[str]) -> Set[str]:
    # TODO consider more ids:
    # ALFRED:LO000002C
    # Comparative Toxicogenomics Database:10017
    # GenAtlas:A1BG
    # GenBank:AY545216.1
    # GO:GO:0000082
    # HumanCyc Gene:HS00051
    # IUPHAR Receptor:114
    # ModBase:A0A183
    # MutDB:A1BG
    # NCBI Gene:100009601
    # RefSeq DNA:NG_000007
    # RefSeq Protein:BAE09150
    # RefSeq RNA:NM_000016
    # UCSC Genome Browser:NM_000014

    # Not used:
    # URL:http://www.imm.ki.se/CYPalleles/cyp3a4.htm

    # Used ids:
    # Ensembl:ENSG00000000460
    # GeneCard:A1BG
    # HGNC:10013
    # OMIM:100640
    # UniProtKB:A0AV47
    filtered_ids = set()
    for gene_id in ids:
        if any([gene_id.startswith(x) for x in ['HGNC:', 'Ensembl:', 'OMIM:', 'GeneCard:', 'UniProtKB:']]):
            filtered_ids.add(gene_id)
    return filtered_ids


def process_disease_external_vocabulary(ids: List[str]) -> Set[Tuple[str, str]]:
    # Used ids:
    # MedDRA:10000247(Abrasion of teeth)
    # MeSH:C536109(N-acetyl glutamate synthetase deficiency)
    # NDFRT:N0000000263(Abetalipoproteinemia [Disease/Finding])
    # SnoMedCT:109769000(Necrotizing sialometaplasia)
    # UMLS:C0001349(Acute-Phase Reaction [Disease/Finding]) or UMLS:C0001261(C0001261)
    filtered_ids = set()
    for disease_id in ids:
        if any([disease_id.startswith(x) for x in ['MedDRA:', 'MeSH:', 'NDFRT:', 'SnoMedCT:', 'UMLS:']]):
            index = disease_id.index('(')
            _id = disease_id[0:index]
            name = disease_id[index + 1:-1]
            filtered_ids.add((_id, None if name == _id.split(':')[1] else name))
    return filtered_ids


network = Network()
request_headers = {
    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_9_3) AppleWebKit/537.36 (KHTML, like Gecko) ' +
                  'Chrome/35.0.1916.47 Safari/537.36'
}
for name in ['genes', 'variants', 'drugs', 'phenotypes']:
    url = 'https://s3.pgkb.org/data/%s.zip' % name
    zip_file = '../data/PharmGKB/%s.zip' % name
    if not os.path.exists(zip_file):
        print('Downloading latest archive for %s...' % name)
        request = urllib.request.Request(url, headers=request_headers)
        with urllib.request.urlopen(request) as response, open(zip_file, 'wb') as f:
            f.write(response.read())
    with zipfile.ZipFile(zip_file) as z:
        with open('../data/PharmGKB/%s.tsv' % name, 'wb') as f:
            f.write(z.read('%s.tsv' % name))

with io.open('../data/PharmGKB/drugs.tsv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    next(reader, None)
    for row in reader:
        # Only parse drugs and not drug classes for now
        if row[5] == 'Drug':
            drug_ids = {'PharmGKB:%s' % row[0]}
            drug_ids.update(process_drug_cross_references(split_list(row[6])))
            for rx_norm_id in split_list(row[21]):
                drug_ids.add('RxNorm:%s' % rx_norm_id)
            for atc_code in split_list(row[22]):
                drug_ids.add('AtcCode:%s' % atc_code)
            for compound_id in split_list(row[23]):
                drug_ids.add('PubChem:CID%s' % compound_id)
            drug = Drug(drug_ids, [row[1]])
            network.add_node(drug)

with io.open('../data/PharmGKB/genes.tsv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    next(reader, None)
    for row in reader:
        gene_ids = {'PharmGKB:%s' % row[0]}
        if row[2]:
            gene_ids.add('HGNC:%s' % row[2])
        for ensembl_id in split_list(row[3]):
            gene_ids.add('Ensembl:%s' % ensembl_id)
        if row[5]:
            gene_ids.add('HGNC:%s' % row[5])
        gene_ids.update(process_gene_cross_references(split_list(row[10])))
        gene = Gene(gene_ids, [row[4]])
        network.add_node(gene)

with io.open('../data/PharmGKB/variants.tsv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    next(reader, None)
    for row in reader:
        variant_ids = {'PharmGKB:%s' % row[0]}
        if row[1]:
            variant_ids.add('dbSNP:%s' % row[1])
        variant = Variant(variant_ids, [])
        network.add_node(variant)

with io.open('../data/PharmGKB/phenotypes.tsv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    next(reader, None)
    for row in reader:
        disease_ids = {'PharmGKB:%s' % row[0]}
        disease_names = {row[1]}
        for id_name_pair in process_disease_external_vocabulary(split_list(row[4])):
            disease_ids.add(id_name_pair[0])
            if id_name_pair[1] is not None:
                disease_names.add(id_name_pair[1])
        disease = Disease(disease_ids, disease_names)
        network.add_node(disease)

with io.open('../data/PharmGKB/graph.json', 'w', encoding='utf-8', newline='') as f:
    f.write(json.dumps(network.to_dict(), indent=2))
