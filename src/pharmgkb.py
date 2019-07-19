#!/usr/bin/env python3

import json
import io
import sys
import csv
import zipfile
from typing import List, Set

from model.gene import Gene
from model.network import Network
from model.drug import Drug
from model.disease import Disease
from model.edge import Edge

maxInt = sys.maxsize
while True:
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)


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


network = Network()
for name in ['genes', 'variants', 'drugs', 'phenotypes']:
    with zipfile.ZipFile('../data/PharmGKB/%s.zip' % name) as z:
        with open('../data/PharmGKB/%s.tsv' % name, 'wb') as f:
            f.write(z.read('%s.tsv' % name))

with io.open('../data/PharmGKB/drugs.tsv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    next(reader, None)
    for row in reader:
        # Only parse drugs and not drug classes for now
        if row[5] == 'Drug':
            drug_ids = {'PharmGKB:%s' % row[0]}
            if row[6] is not None and len(row[6]) > 0:
                drug_ids.update(process_drug_cross_references(row[6].split('","')))
            if row[21] is not None and len(row[21]) > 0:
                for rx_norm_id in row[21].split('","'):
                    drug_ids.add('RxNorm:%s' % rx_norm_id)
            if row[22] is not None and len(row[22]) > 0:
                for atc_code in row[22].split('","'):
                    drug_ids.add('AtcCode:%s' % atc_code)
            if row[23] is not None and len(row[23]) > 0:
                for compound_id in row[23].split('","'):
                    drug_ids.add('PubChem:CID%s' % compound_id)
            drug = Drug(drug_ids, [row[1]])
            network.add_node(drug)

with io.open('../data/PharmGKB/genes.tsv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    next(reader, None)
    for row in reader:
        gene_ids = {'PharmGKB:%s' % row[0]}
        gene = Gene(gene_ids, [row[4]])
        # TODO
        network.add_node(gene)
        pass

with io.open('../data/PharmGKB/graph.json', 'w', encoding='utf-8', newline='') as f:
    f.write(json.dumps(network.to_dict(), indent=2))
