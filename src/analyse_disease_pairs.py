#!/usr/bin/env python3
# pip install neo4j-driver
from neo4j import GraphDatabase, basic_auth
import json
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import io
from model.gene import Gene
from model.edge import Edge
from model.variant import Variant
from model.disease import Disease
from model.drug import Drug
from model.network import Network
from model.circrna import CircRNA
from model.erna import ERNA
from model.lncrna import LncRNA
from model.mirna import MiRNA
from model.ncrna import NcRNA
from model.pirna import PiRNA
from model.pseudogene import Pseudogene
from model.ribozyme import Ribozyme
from model.rna import RNA
from model.rrna import RRNA
from model.scarna import ScaRNA
from model.scrna import ScRNA
from model.snorna import SnoRNA
from model.snrna import SnRNA

# ID pairs of all diseases, returns array in which each pair is an array
# writes pairs with IDs and names to file if the function is called with the parameter True
def get_disease_pairs(writing_files):
    query = """ MATCH (d:Disease) RETURN d.`_id`, d.names """
    results = session.run(query)
    diseases = []
    for disease in results:
        diseases.append([disease['d.`_id`'], disease['d.names']])
    print('Number of diseases: ', len(diseases))
    disease_pairs = []
    writing_pairs = []
    for d1 in range(len(diseases)):
        for d2 in range(d1 + 1, len(diseases)):
            disease_pairs.append([diseases[d1][0], diseases[d2][0]])
            if writing_files:
                writing_pairs.append([diseases[d1][0], diseases[d1][1],  diseases[d2][0], diseases[d2][1]])

    if writing_files:
        with io.open('../analysis/disease_pairs/disease_pairs.tsv', 'w', encoding='utf-8', newline='') as disease_pairs_file:
            disease_pairs_file.write('#Disease A ID\tDisease A names\tDisease B ID\tDisease B names\n')
            for pair in writing_pairs:
                d1_names = str(pair[1])
                d1_names = d1_names.replace('[', '')
                d1_names = d1_names.replace(']', '')
                d1_names = d1_names.replace('\'', '')
                d2_names = str(pair[3])
                d2_names = d2_names.replace('[', '')
                d2_names = d2_names.replace(']', '')
                d2_names = d2_names.replace('\'', '')
                disease_pairs_file.write(pair[0] + '\t' + d1_names + '\t' + pair[2] + '\t' + d2_names + '\n')

    print('Number of disease pairs: ', len(disease_pairs))
    return disease_pairs


# returns two sets of gene IDs, for each disease of the disease pair one set with the disease associated gene IDs is returned
def get_genes(disease_pair):
    d1_genes = set()
    query = """ MATCH p=(d:Disease)-[a]-(g:Gene) WHERE {d1_id} IN d.ids RETURN g.`_id` """
    d1_results = session.run(query, parameters={'d1_id': disease_pair[0]})
    for gene in d1_results:
        d1_genes.add(gene['g.`_id`'])

    d2_genes = set()
    query = """ MATCH p=(d:Disease)-[a]-(g:Gene) WHERE {d2_id} IN d.ids RETURN g.`_id` """
    d2_results = session.run(query, parameters={'d2_id': disease_pair[1]})
    for gene in d2_results:
        d2_genes.add(gene['g.`_id`'])

    return d1_genes, d2_genes


# an array is returned, in which for each disease pair a network is stored
# in the network common genes are associated to both diseases
# if the parameter True is handed over, the common genes are written to file
def get_common_genes(disease_pairs, networks, writing_files):
    new_networks = []
    for index, disease_pair in enumerate(disease_pairs):
        network = networks[index]
        d1_genes, d2_genes = get_genes(disease_pair)
        common_genes = d1_genes.intersection(d2_genes)
        d1 = Disease([disease_pair[0]], [])
        network.add_node(d1)
        d2 = Disease([disease_pair[1]], [])
        network.add_node(d2)
        for g_id in common_genes:
            gene = Gene([g_id], [])
            network.add_node(gene)
            network.add_edge(Edge(gene, d1, 'ASSOCIATES_WITH', {}))
            network.add_edge(Edge(gene, d2, 'ASSOCIATES_WITH', {}))
        if len(common_genes) > 0 and writing_files:
            temp_id1 = disease_pair[0].replace(':', '-')
            temp_id2 = disease_pair[1].replace(':', '-')
            path = '../analysis/disease_pairs/' + temp_id1 + '_' + temp_id2
            try:
                os.mkdir(path)
            except FileExistsError:
                pass
            with io.open(path + '/' + temp_id1 + '_' + temp_id2 + '_common_genes.tsv', 'w', encoding='utf-8', newline='') as common_genes_file:
                common_genes_file.write('#Common genes of ' + disease_pair[0] + ' and ' + disease_pair[1] + '\n')
                for gene in common_genes:
                    common_genes_file.write(gene + '\n')
        new_networks.append(network)
    print('Done getting genes')
    return new_networks


# an array is returned, in which for each disease pair a network is stored
# in the network common drugs are associated to both diseases or the drug is associated to the genes of the disease
# if the parameter True is handed over, the common drugs are written to file
def get_common_drugs(disease_pairs, networks, writing_files):
    new_networks = []
    for index, disease_pair in enumerate(disease_pairs):
        network = networks[index]
        d1 = Disease([disease_pair[0]], [])
        network.add_node(d1)
        d2 = Disease([disease_pair[1]], [])
        network.add_node(d2)

        # the drug INDICATES, CONTRAINDICATES or INDUCES both diseases
        common_drugs = set()
        query = """ MATCH (d1:Disease)-[a]-(n:Drug)--(d2:Disease) WHERE {d1_id} IN d1.ids AND {d2_id} IN d2.ids RETURN distinct(type(a)), n.`_id` """
        results = session.run(query, parameters={'d1_id': disease_pair[0], 'd2_id': disease_pair[1]})
        for result in results:
            drug_id = result['n.`_id`']
            type = result['(type(a))']
            common_drugs.add(drug_id)
            drug = Drug([drug_id], [])
            network.add_node(drug)
            network.add_edge(Edge(drug, d1, type, {}))
        query = """ MATCH (d1:Disease)--(n:Drug)-[a]-(d2:Disease) WHERE {d1_id} IN d1.ids AND {d2_id} IN d2.ids RETURN distinct(type(a)), n.`_id` """
        results = session.run(query, parameters={'d1_id': disease_pair[0], 'd2_id': disease_pair[1]})
        for result in results:
            drug_id = result['n.`_id`']
            type = result['(type(a))']
            common_drugs.add(drug_id)
            drug = Drug([drug_id], [])
            network.add_node(drug)
            network.add_edge(Edge(drug, d2, type, {}))

        # the drug targets a gene of one disease and is associated to the other disease
        query = """ MATCH (d1:Disease)-[a]-(n:Drug)-[:TARGETS]-(g:Gene)-[:ASSOCIATES_WITH]-(d2:Disease) WHERE {d1_id} IN d1.ids AND {d2_id} IN d2.ids RETURN distinct(type(a)), n.`_id`, g.`_id` """
        results = session.run(query, parameters={'d1_id': disease_pair[0], 'd2_id': disease_pair[1]})
        for result in results:
            drug_id = result['n.`_id`']
            type = result['(type(a))']
            common_drugs.add(drug_id)
            drug = Drug([drug_id], [])
            network.add_node(drug)
            network.add_edge(Edge(drug, d1, type, {}))
            gene_id = result['g.`_id`']
            gene = Gene([gene_id], [])
            network.add_node(gene)
            network.add_edge(Edge(drug, gene, 'TARGETS', {'actions': []}))
            network.add_edge(Edge(gene, d2, 'ASSOCIATES_WITH', {}))
        query = """ MATCH (d2:Disease)-[a]-(n:Drug)-[:TARGETS]-(g:Gene)-[:ASSOCIATES_WITH]-(d1:Disease) WHERE {d1_id} IN d1.ids AND {d2_id} IN d2.ids RETURN distinct(type(a)), n.`_id`, g.`_id` """
        results = session.run(query, parameters={'d1_id': disease_pair[0], 'd2_id': disease_pair[1]})
        for result in results:
            drug_id = result['n.`_id`']
            type = result['(type(a))']
            common_drugs.add(drug_id)
            drug = Drug([drug_id], [])
            network.add_node(drug)
            network.add_edge(Edge(drug, d2, type, {}))
            gene_id = result['g.`_id`']
            gene = Gene([gene_id], [])
            network.add_node(gene)
            network.add_edge(Edge(drug, gene, 'TARGETS', {'actions': []}))
            network.add_edge(Edge(gene, d1, 'ASSOCIATES_WITH', {}))

        # the drug targets one gene which is associated to both diseases or the drug targets two different genes
        # where each gene is associated to one of the diseases
        query = """ MATCH (d1:Disease)-[:ASSOCIATES_WITH]-(g1:Gene)-[:TARGETS]-(n:Drug)-[:TARGETS]-(g2:Gene)-
        [:ASSOCIATES_WITH]-(d2:Disease) WHERE {d1_id} IN d1.ids AND {d2_id} IN d2.ids RETURN n.`_id`, g1.`_id`, g2.`_id` """
        results = session.run(query, parameters={'d1_id': disease_pair[0], 'd2_id': disease_pair[1]})
        for result in results:
            drug_id = result['n.`_id`']
            common_drugs.add(drug_id)
            g1_id = result['g1.`_id`']
            g2_id = result['g2.`_id`']
            g1 = Gene([g1_id], [])
            network.add_node(g1)
            network.add_edge(Edge(g1, d1, 'ASSOCIATES_WITH', {}))
            drug = Drug([drug_id], [])
            network.add_node(drug)
            network.add_edge(Edge(drug, g1, 'TARGETS', {'actions': []}))
            g2 = Gene([g2_id], [])
            network.add_node(g2)
            network.add_edge(Edge(drug, g2, 'TARGETS', {'actions': []}))
            network.add_edge(Edge(g2, d2, 'ASSOCIATES_WITH', {}))

        new_networks.append(network)

        if len(common_drugs) > 0 and writing_files:
            temp_id1 = disease_pair[0].replace(':', '-')
            temp_id2 = disease_pair[1].replace(':', '-')
            path = '../analysis/disease_pairs/' + temp_id1 + '_' + temp_id2
            try:
                os.mkdir(path)
            except FileExistsError:
                pass
            with io.open(path + '/' + temp_id1 + '_' + temp_id2 + '_common_drugs.tsv', 'w', encoding='utf-8', newline='') as common_drugs_file:
                common_drugs_file.write('#Common drugs of ' + disease_pair[0] + ' and ' + disease_pair[1] + '\n')
                for drug in common_drugs:
                    common_drugs_file.write(drug + '\n')
    print('Done getting drugs')
    return new_networks


# an array is returned, in which for each disease pair a network is stored
# in the network common RNAs regulate at least a common gene, two disease specific genes or a common and a disease specific gene
# if the parameter True is handed over, the common RNAs are written to file as well as the number and names of the genes
# they regulate and the RNAs which regulate the RNAs
def get_common_rnas(disease_pairs, networks, writing_files):
    new_networks = []
    for index, disease_pair in enumerate(disease_pairs):
        network = networks[index]
        d1 = Disease([disease_pair[0]], [])
        network.add_node(d1)
        d2 = Disease([disease_pair[1]], [])
        network.add_node(d2)

        d1_genes_ids, d2_genes_ids = get_genes(disease_pair)

        # this differentiation is done to get the correct number of regulated, in this subgraph present genes
        d1_only_genes_ids = d1_genes_ids.difference(d2_genes_ids)
        d2_only_genes_ids = d2_genes_ids.difference(d1_genes_ids)
        common_genes_ids = d1_genes_ids.intersection(d2_genes_ids)

        common_rnas = {}    #dict with the RNA name as key and the regulated genes as an array as value
        for gene_id in common_genes_ids:
            query = """ MATCH (g:Gene)-[:REGULATES]-(r:RNA) WHERE {gene_id} IN g.ids RETURN distinct(r.`_id`) """
            results = session.run(query, parameters={'gene_id': gene_id})
            for result in results:
                rna_id = result['(r.`_id`)']
                if rna_id in common_rnas:
                    gene_ids = common_rnas[rna_id]
                    gene_ids.append(gene_id)
                    common_rnas[rna_id] = gene_ids
                else:
                    common_rnas[rna_id] = [gene_id]
                    gene = Gene([gene_id], [])
                    network.add_node(gene)
                    rna = RNA([rna_id], [])
                    network.add_node(rna)
                    network.add_edge(Edge(rna, gene, 'REGULATES', {}))
                    network.add_edge(Edge(gene, d1, 'ASSOCIATES_WITH', {}))
                    network.add_edge(Edge(gene, d2, 'ASSOCIATES_WITH', {}))

        rnas_d1_only_genes = {}
        for gene_id in d1_only_genes_ids:
            query = """ MATCH (g:Gene)-[:REGULATES]-(r:RNA) WHERE {gene_id} IN g.ids RETURN distinct(r.`_id`) """
            results = session.run(query, parameters={'gene_id': gene_id})
            for result in results:
                rna_id = result['(r.`_id`)']
                if rna_id in rnas_d1_only_genes:
                    gene_ids = rnas_d1_only_genes[rna_id]
                    gene_ids.append(gene_id)
                    rnas_d1_only_genes[rna_id] = gene_ids
                else:
                    rnas_d1_only_genes[rna_id] = [gene_id]

        rnas_d2_only_genes = {}
        for gene_id in d2_only_genes_ids:
            query = """ MATCH (g:Gene)-[:REGULATES]-(r:RNA) WHERE {gene_id} IN g.ids RETURN distinct(r.`_id`) """
            results = session.run(query, parameters={'gene_id': gene_id})
            for result in results:
                rna_id = result['(r.`_id`)']
                if rna_id in rnas_d2_only_genes:
                    gene_ids = rnas_d2_only_genes[rna_id]
                    gene_ids.append(gene_id)
                    rnas_d2_only_genes[rna_id] = gene_ids
                else:
                    rnas_d2_only_genes[rna_id] = [gene_id]

        #common_rnas = {'A':1, 'B':1, 'D':1}
        #rnas_d1_only_genes = {'A':2, 'B':1, 'E':1}
        #rnas_d2_only_genes = {'A':2, 'C':1, 'E':1}

        for rna_id in rnas_d1_only_genes:
            if rna_id in common_rnas:
                # common_rnas have already been added to the network, here the number of regulated genes is updated
                common_rnas[rna_id] = common_rnas[rna_id] + rnas_d1_only_genes[rna_id]
            elif rna_id in rnas_d2_only_genes:
                # RNA regulates genes associated to d1 and genes associated to d2, RNA does not regulate a common gene
                common_rnas[rna_id] = rnas_d1_only_genes[rna_id] + rnas_d2_only_genes[rna_id]
                g1_ids = rnas_d1_only_genes[rna_id]
                g2_ids = rnas_d2_only_genes[rna_id]
                rna = RNA([rna_id], [])
                network.add_node(rna)
                for g_id in g1_ids:
                    gene = Gene([g_id], [])
                    network.add_node(gene)
                    network.add_edge(Edge(gene, d1, 'ASSOCIATES_WITH', {}))
                    network.add_edge(Edge(rna, gene, 'REGULATES', {}))
                for g_id in g2_ids:
                    gene = Gene([g_id], [])
                    network.add_node(gene)
                    network.add_edge(Edge(gene, d2, 'ASSOCIATES_WITH', {}))
                    network.add_edge(Edge(rna, gene, 'REGULATES', {}))
                del rnas_d2_only_genes[rna_id]
        for rna_id in rnas_d2_only_genes:
            if rna_id in common_rnas:
                # common_rnas have already been added to the network, here the number of regulated genes is updated
                common_rnas[rna_id] = common_rnas[rna_id] + rnas_d2_only_genes[rna_id]

        # for each RNA add an array of RNAs, which regulate this RNA. MRNAs are not included
        for rna_id in common_rnas:
            second_rnas = []
            query = """MATCH (r:RNA)-[:REGULATES]-(n:RNA) WHERE {r_id} IN r.ids AND NOT n.label_id CONTAINS "MRNA" RETURN distinct(n.`_id`) """
            results = session.run(query, parameters={'r_id': rna_id})
            rna = RNA([rna_id], [])
            network.add_node(rna)
            for result in results:
                second_rna_id = result['(n.`_id`)']
                second_rnas.append(second_rna_id)
                second_rna = RNA([second_rna_id], [])
                network.add_node(second_rna)
                network.add_edge(Edge(second_rna, rna, 'REGULATES', {}))
            # the value of common_rnas is now changed to an array where at the first position the array with the regulated
            # genes from this subgraph is stored and at the second position the array with RNAs regulating the RNA is stored
            common_rnas[rna_id] = [common_rnas[rna_id], second_rnas]

        new_networks.append(network)

        if len(common_rnas) > 0 and writing_files:
            temp_id1 = disease_pair[0].replace(':', '-')
            temp_id2 = disease_pair[1].replace(':', '-')
            path = '../analysis/disease_pairs/' + temp_id1 + '_' + temp_id2
            try:
                os.mkdir(path)
            except FileExistsError:
                pass
            with io.open(path + '/' + temp_id1 + '_' + temp_id2 + '_common_rnas.tsv', 'w', encoding='utf-8', newline='') as common_rnas_file:
                common_rnas_file.write('#Common rnas of ' + disease_pair[0] + ' and ' + disease_pair[1] + '\tsorted by number of regulated genes\tRegulated genes\tRNAs regulating the RNA\n')
                for key, value in sorted(common_rnas.items(), key=lambda item: len(item[1][0]), reverse=True):
                    # sort by the number of genes in this subgraph which are regulated by the RNA
                    regulated_genes = str(value[0])
                    regulated_genes = regulated_genes.replace('[', '')
                    regulated_genes = regulated_genes.replace(']', '')
                    regulated_genes = regulated_genes.replace('\'', '')
                    second_rnas = str(value[1])
                    second_rnas = second_rnas.replace('[', '')
                    second_rnas = second_rnas.replace(']', '')
                    second_rnas = second_rnas.replace('\'', '')
                    common_rnas_file.write(key + '\t' + str(len(value[0])) + '\t' + regulated_genes + '\t' + second_rnas + '\n')
    print('Done getting RNAs')
    return new_networks


# an array is returned, in which for each disease pair a network is stored
# in the network common variants can be disease or gene associated
# if the parameter True is handed over, the common variants are written to file
def get_common_variants(disease_pairs, networks, writing_files):
    new_networks = []
    for index, disease_pair in enumerate(disease_pairs):
        network = networks[index]
        d1 = Disease([disease_pair[0]], [])
        network.add_node(d1)
        d2 = Disease([disease_pair[1]], [])
        network.add_node(d2)

        common_variants = []    # each variant is an array
        query = """ MATCH (d1:Disease)--(v:Variant)--(d2:Disease) WHERE {d1_id} in d1.ids AND {d2_id} in d2.ids RETURN v.`_id` """
        results = session.run(query, parameters={'d1_id': disease_pair[0], 'd2_id': disease_pair[1]})
        for result in results:
            v_id = result['v.`_id`']
            common_variants.append([v_id, 'disease associated'])
            variant = Variant([v_id], [])
            network.add_node(variant)
            network.add_edge(Edge(d1, variant, 'ASSOCIATES_WITH', {}))
            network.add_edge(Edge(d2, variant, 'ASSOCIATES_WITH', {}))

        # variants associated to common genes
        d1_genes, d2_genes = get_genes(disease_pair)
        common_genes_ids = d1_genes.intersection(d2_genes)
        for gene_id in common_genes_ids:
            query = """ MATCH (g:Gene)-[a]-(v:Variant) WHERE {g_id} in g.ids RETURN v.`_id`, type(a) """
            results = session.run(query, parameters={'g_id': gene_id})
            for result in results:
                v_id = result['v.`_id`']
                type = result['type(a)']    # can be CODES or EQTL
                variant_pair = v_id + '-' + gene_id
                common_variants.append([variant_pair, 'gene associated'])
                variant = Variant([v_id], [])
                network.add_node(variant)
                gene = Gene([gene_id], [])
                network.add_node(gene)
                network.add_edge(Edge(gene, variant, type, {}))
                network.add_edge(Edge(gene, d1, 'ASSOCIATES_WITH', {}))
                network.add_edge(Edge(gene, d2, 'ASSOCIATES_WITH', {}))

        new_networks.append(network)

        if len(common_variants) > 0 and writing_files:
            temp_id1 = disease_pair[0].replace(':', '-')
            temp_id2 = disease_pair[1].replace(':', '-')
            path = '../analysis/disease_pairs/' + temp_id1 + '_' + temp_id2
            try:
                os.mkdir(path)
            except FileExistsError:
                pass
            with io.open(path + '/' + temp_id1 + '_' + temp_id2 + '_common_variants.tsv', 'w', encoding='utf-8',
                         newline='') as common_variants_file:
                common_variants_file.write(
                    '#Common variants associated with ' + disease_pair[0] + ' and ' + disease_pair[1] + '\n')
                for variant in common_variants:
                    common_variants_file.write(variant[0] + '\t' + variant[1] + '\n')
    print('Done getting variants')
    return new_networks


# returns all information in an array
# for each disease pair a network with common genes, common drugs, common RNAs and common variant is returned
# and the network is saved as a json file and thus can be used as database
# if the parameter True is handed over, all information is written to files
def get_disease_pairs_info(disease_pairs, writing_files):
    networks = []
    for disease_pair in disease_pairs:
        networks.append(Network())
    networks = get_common_genes(disease_pairs, networks, writing_files)
    networks = get_common_drugs(disease_pairs, networks, writing_files)
    networks = get_common_rnas(disease_pairs, networks, writing_files)
    networks = get_common_variants(disease_pairs, networks, writing_files)
    if writing_files:
        for index, disease_pair in enumerate(disease_pairs):
            temp_id1 = disease_pair[0].replace(':', '-')
            temp_id2 = disease_pair[1].replace(':', '-')
            path = '../analysis/disease_pairs/' + temp_id1 + '_' + temp_id2
            try:
                os.mkdir(path)
            except FileExistsError:
                pass
            network = networks[index]
            network.save(path + '/' + temp_id1 + '_' + temp_id2 + '_full_graph.json')
    return networks


# an array is returned, in which for each disease pair an array with a network for each drug subgraph is stored
# the network is a drug related subgraph with the most important genes, RNAs and variants
# thus, the number of visualized variants, RNAs and genes can be limited
# all relevant information is written to a results file, a json file and visualised as a PNG file
def get_given_drugs_related_info(disease_pairs, drugs):   # first disease pair with first drug array
    all_networks = []   # contains an array for each disease pair
    for index, disease_pair in enumerate(disease_pairs):
        networks_per_drug = []  # contains a network for each drug
        pair_drugs_ids = drugs[index]
        temp_id1 = disease_pair[0].replace(':', '-')
        temp_id2 = disease_pair[1].replace(':', '-')
        path = '../analysis/disease_pairs/' + temp_id1 + '_' + temp_id2
        for drug_id in pair_drugs_ids:
            try:
                os.mkdir(path)
            except FileExistsError:
                pass

            network = Network()
            d1 = Disease([disease_pair[0]], [])
            network.add_node(d1)
            d2 = Disease([disease_pair[1]], [])
            network.add_node(d2)
            drug = Drug([drug_id], [])
            network.add_node(drug)
            temp_drug_id = drug_id.replace(':', '-')
            with io.open(path + '/' + temp_id1 + '_' + temp_id2 + '_' + temp_drug_id + '_results.txt', 'w', encoding='utf-8', newline='') as results_file:
                results_file.write('In this file all information about the connection between ' + disease_pair[0] +
                                   ' and ' + disease_pair[1] + ' and the drug ' + drug_id + ' is summarized:\n')

                # the drug INDICATES, CONTRAINDICATES or INDUCES the disease
                query = """ MATCH (d:Disease)-[a]-(n:Drug) WHERE {d1_id} IN d.ids AND {n_id} in n.ids RETURN distinct(type(a)) """
                d1_results = session.run(query, parameters={'d1_id': disease_pair[0], 'n_id': drug_id})
                for result in d1_results:
                    results_file.write(drug_id + ' ' + result['(type(a))'] + ' ' + disease_pair[0] + '\n')
                    network.add_edge(Edge(drug, d1, result['(type(a))'], {}))
                query = """ MATCH (d:Disease)-[a]-(n:Drug) WHERE {d2_id} IN d.ids AND {n_id} in n.ids RETURN distinct(type(a)) """
                d2_results = session.run(query, parameters={'d2_id': disease_pair[1], 'n_id': drug_id})
                for result in d2_results:
                    results_file.write(drug_id + ' ' + result['(type(a))'] + ' ' + disease_pair[1] + '\n')
                    network.add_edge(Edge(drug, d2, result['(type(a))'], {}))

                # the drug targets a gene which is associated to the disease
                d1_genes = set()
                query = """ MATCH (n:Drug)-[:TARGETS]-(g:Gene)-[:ASSOCIATES_WITH]-(d:Disease) WHERE {d1_id} IN d.ids AND {n_id} in n.ids RETURN g.`_id` """
                d1_results = session.run(query, parameters={'d1_id': disease_pair[0], 'n_id': drug_id})
                for gene in d1_results:
                    d1_genes.add(gene['g.`_id`'])
                    g = Gene([gene['g.`_id`']], [])
                    network.add_node(g)
                    network.add_edge(Edge(drug, g, 'TARGETS', {'actions': []})) #TODO
                    network.add_edge(Edge(g, d1, 'ASSOCIATES_WITH', {}))
                d2_genes = set()
                query = """ MATCH (n:Drug)-[:TARGETS]-(g:Gene)-[:ASSOCIATES_WITH]-(d:Disease) WHERE {d2_id} IN d.ids AND {n_id} in n.ids RETURN g.`_id` """
                d2_results = session.run(query, parameters={'d2_id': disease_pair[1], 'n_id': drug_id})
                for gene in d2_results:
                    d2_genes.add(gene['g.`_id`'])
                    g = Gene([gene['g.`_id`']], [])
                    network.add_node(g)
                    network.add_edge(Edge(drug, g, 'TARGETS', {'actions': []})) #TODO
                    network.add_edge(Edge(g, d2, 'ASSOCIATES_WITH', {}))

                common_drug_genes = d1_genes.intersection(d2_genes) # genes associated to the drug and both diseases
                # relevant_genes are all genes associated to at least one disease and the drug, below the common genes
                # with the most disease associated references are added
                relevant_genes = d1_genes.union(d2_genes)
                if len(d1_genes) > 0:
                    nbr = str(len(d1_genes))
                    d1_genes = str(d1_genes)
                    d1_genes = d1_genes.replace('{', '')
                    d1_genes = d1_genes.replace('}', '')
                    d1_genes = d1_genes.replace('\'', '')
                    results_file.write(drug_id + ' targets following ' + nbr + ' genes which are associated to ' + disease_pair[0] + ': ' + d1_genes + '\n')
                if len(d2_genes) > 0:
                    nbr = str(len(d2_genes))
                    d2_genes = str(d2_genes)
                    d2_genes = d2_genes.replace('{', '')
                    d2_genes = d2_genes.replace('}', '')
                    d2_genes = d2_genes.replace('\'', '')
                    results_file.write(drug_id + ' targets following ' + nbr + ' genes which are associated to ' + disease_pair[1] + ': ' + d2_genes + '\n')
                if len(common_drug_genes) > 0:
                    nbr = str(len(common_drug_genes))
                    cdgs = str(common_drug_genes)
                    cdgs = cdgs.replace('{', '')
                    cdgs = cdgs.replace('}', '')
                    cdgs = cdgs.replace('\'', '')
                    results_file.write('The disease pair has ' + nbr + ' common genes which are targeted by the drug: ' + cdgs + '\n')

                # add the common genes with the most disease associated references
                # no given num_pmids is similar to num_pmids = 0
                all_d1_genes, all_d2_genes = get_genes(disease_pair)
                all_common_genes = all_d1_genes.intersection(all_d2_genes)
                relevant_common_genes = []  # the genes with the most cited gene-disease association, threshold 10
                if len(all_common_genes) > 0:
                    results_file.write('The disease pair has ' + str(len(all_common_genes)) + ' common genes, not considering the connection to the drug.'
                                        ' Following genes have the most references regarding their connection to both diseases:\n')
                    for gene in all_common_genes:
                        query = """ MATCH (d1:Disease)-[a]-(g:Gene) WHERE {g_id} IN g.ids AND {d1_id} IN d1.ids RETURN a.num_pmids """
                        results = session.run(query, parameters={'g_id': gene, 'd1_id': disease_pair[0]})
                        num_pmids = 0
                        for result in results:  # multiple edges to the same gene
                            temp = result['a.num_pmids']
                            if temp is not None:
                                num_pmids = num_pmids + temp
                        query = """ MATCH (d2:Disease)-[a]-(g:Gene) WHERE {g_id} IN g.ids AND {d2_id} IN d2.ids RETURN a.num_pmids """
                        results = session.run(query, parameters={'g_id': gene, 'd2_id': disease_pair[1]})
                        for result in results:  # multiple edges to the same gene
                            temp = result['a.num_pmids']
                            if temp is not None:
                                num_pmids = num_pmids + temp
                        relevant_common_genes.append([gene, num_pmids])
                    # sort by number of pmids
                    relevant_common_genes = sorted(relevant_common_genes, key=lambda item: item[1], reverse=True)
                    relevant_common_genes = relevant_common_genes[:10]  # threshold
                    rcgs = str(relevant_common_genes)
                    rcgs = rcgs[1:-1]
                    rcgs = rcgs.replace('\'', '')
                    results_file.write(rcgs + '\n')
                    for g in relevant_common_genes:
                        gene = Gene([g[0]], [])
                        network.add_node(gene)
                        network.add_edge(Edge(gene, d1, 'ASSOCIATES_WITH', {}))
                        network.add_edge(Edge(gene, d2, 'ASSOCIATES_WITH', {}))
                        relevant_genes.add(g[0])

                # add the common disease associated variants with most references
                # no given num_pmids is similar to num_pmids = 0
                disease_variants = {}
                query = """ MATCH (d1:Disease)-[a]-(v:Variant)--(d2:Disease) WHERE {d1_id} in d1.ids AND {d2_id} in d2.ids RETURN distinct(a.num_pmids), v.`_id` """
                results = session.run(query, parameters={'d1_id': disease_pair[0], 'd2_id': disease_pair[1]})
                for variant in results:
                    num_pmids = variant['(a.num_pmids)']
                    if num_pmids is None:
                        num_pmids = 0
                    var_id = variant['v.`_id`']
                    if var_id in disease_variants:
                        temp = disease_variants[var_id]
                        disease_variants[var_id] = temp + num_pmids
                    else:
                        disease_variants[var_id] = num_pmids
                query = """ MATCH (d2:Disease)-[a]-(v:Variant)--(d1:Disease) WHERE {d1_id} in d1.ids AND {d2_id} in d2.ids RETURN distinct(a.num_pmids), v.`_id` """
                results = session.run(query, parameters={'d1_id': disease_pair[0], 'd2_id': disease_pair[1]})
                for variant in results:
                    num_pmids = variant['(a.num_pmids)']
                    if num_pmids is None:
                        num_pmids = 0
                    var_id = variant['v.`_id`']
                    if var_id in disease_variants:
                        temp = disease_variants[var_id]
                        disease_variants[var_id] = temp + num_pmids
                    else:
                        disease_variants[var_id] = num_pmids
                dvs = ''
                i = 0
                for key, value in sorted(disease_variants.items(), key=lambda item: item[1], reverse=True):
                    if i < 9:   # threshold
                        num_pmids = disease_variants[key]
                        variant = Variant([key], [])
                        network.add_node(variant)
                        network.add_edge(Edge(variant, d1, 'ASSOCIATES_WITH', {}))
                        network.add_edge(Edge(variant, d2, 'ASSOCIATES_WITH', {}))
                        dvs = dvs + key + ':' + str(num_pmids) + ' PMIDs, '
                        i += 1
                dvs = dvs[:-2]

                # add the gene associated variants with smallest pvalues
                # if no pvalue is given, pvalue is set to 1
                gene_variants = []
                for gene in relevant_genes:
                    query = """ MATCH (g:Gene)-[a]-(v:Variant) WHERE {g_id} in g.ids RETURN v.`_id`, a.pvalue, type(a) """
                    results = session.run(query, parameters={'g_id': gene})
                    for variant in results:
                        pvalue = variant['a.pvalue']
                        if pvalue is None:
                            pvalue = 1
                        else:
                            pvalue = float(pvalue)
                        gene_variants.append([variant['v.`_id`'] + '-' + gene, pvalue, variant['type(a)']])
                gene_variants = sorted(gene_variants, key=lambda item: item[1])
                gene_variants = gene_variants[:10]  # threshold
                for v in gene_variants:
                    temp = v[0].split('-')
                    v_id = temp[0]
                    g_id = temp[1]
                    variant = Variant([v_id], [])
                    network.add_node(variant)
                    gene = Gene([g_id], [])
                    network.add_node(gene)
                    network.add_edge(Edge(gene, variant, v[2], {'pvalue': v[1]}))
                if len(gene_variants) > 0:
                    gvs = str(gene_variants)
                    gvs = gvs[1:-1]
                    gvs = gvs.replace('\'', '')
                else:
                    gvs = ''

                if len(disease_variants) > 0 or len(gene_variants) > 0:
                    results_file.write('The disease pair has at least ' + str(i) + ' variants associated to both diseases: ' +
                                           dvs + ' and at least ' + str(len(gene_variants)) + ' gene associated variants: ' + gvs + '\n')

                # dict with RNA name as key and an array as value
                # first array position is the number of regulated genes, second position is an array with the gene names
                relevant_rnas = {}
                for gene in relevant_genes:
                    query = """ MATCH (r:RNA)--(g:Gene) WHERE {g_id} in g.ids AND NOT r.label_id CONTAINS "MRNA" return r.`_id` """
                    results = session.run(query, parameters={'g_id': gene})
                    for result in results:
                        key = result['r.`_id`']
                        if key in relevant_rnas:
                            value = relevant_rnas[key]
                            genes = value[1]
                            if gene not in genes:
                                genes.add(gene)
                                relevant_rnas[key] = [value[0] + 1, genes]
                        else:
                            genes = set()
                            genes.add(gene)
                            relevant_rnas[key] = [1, genes]

                if len(relevant_rnas) > 0:
                    i = 0
                    for key, value in sorted(relevant_rnas.items(), key=lambda item: item[1], reverse=True):
                    # sort by the number of regulated genes
                        if i > 9:   # threshold
                            break
                        elif value[0] > 1:  # only add and print RNAs which regulate more than one gene
                            if i == 0:
                                results_file.write('RNAs with the number and names of the genes they regulate: \n')
                            rna_id = key
                            for gene_id in value[1]:
                                rna = RNA([rna_id], [])
                                network.add_node(rna)
                                gene = Gene([gene_id], [])
                                network.add_node(gene)
                                network.add_edge(Edge(rna, gene, 'REGULATES', {}))
                            regulated_genes = str(value[1])
                            regulated_genes = regulated_genes[1:-1]
                            regulated_genes = regulated_genes.replace('\'', '')
                            results_file.write(rna_id + '\t' + str(value[0]) + '\t' + regulated_genes + '\n')
                            i += 1

                    # append regulating RNAs to one RNA which regulates the most genes, MRNAs are not added
                    for key, value in sorted(relevant_rnas.items(), key=lambda item: item[1], reverse=True):
                        if value[0] > 1:
                            most_relevant_rna = RNA([key], [])
                            network.add_node(most_relevant_rna)
                            query = """ MATCH (r:RNA)--(n:RNA) WHERE {r_id} in r.ids AND NOT n.label_id CONTAINS "MRNA" RETURN n.`_id`, labels(n) """
                            results = session.run(query, parameters={'r_id': key})
                            reg_rnas = ''
                            for result in results:
                                rna_id = result['n.`_id`']
                                types = result['labels(n)']
                                for type in types:
                                    if type != 'RNA':
                                        if type == 'CircRNA':
                                            rna = CircRNA([rna_id], [])
                                        if type == 'ERNA':
                                            rna = ERNA([rna_id], [])
                                        if type == 'LncRNA':
                                            rna = LncRNA([rna_id], [])
                                        if type == 'MiRNA':
                                            rna = MiRNA([rna_id], [])
                                        if type == 'NcRNA':
                                            rna = NcRNA([rna_id], [])
                                        if type == 'PiRNA':
                                            rna = PiRNA([rna_id], [])
                                        if type == 'Pseudogene':
                                            rna = Pseudogene([rna_id], [])
                                        if type == 'Ribozyme':
                                            rna = Ribozyme([rna_id], [])
                                        if type == 'RRNA':
                                            rna = RRNA([rna_id], [])
                                        if type == 'ScaRNA':
                                            rna = ScaRNA([rna_id], [])
                                        if type == 'ScRNA':
                                            rna = ScRNA([rna_id], [])
                                        if type == 'SnoRNA':
                                            rna = SnoRNA([rna_id], [])
                                        if type == 'SnRNA':
                                            rna = SnRNA([rna_id], [])
                                        network.add_node(rna)
                                        network.add_edge(Edge(rna, most_relevant_rna, 'REGULATES', {}))
                                        reg_rnas = reg_rnas + rna_id + ', '
                            reg_rnas = reg_rnas[:-2]
                            results_file.write(key + ' is the RNA which regulates the most genes in this subgraph. It is regulated by ' + reg_rnas + '.\n')
                        break
            json_file = path + '/' + temp_id1 + '_' + temp_id2 + '_' + temp_drug_id + '_graph.json'
            network.save(json_file)
            draw_drug_subgraph(json_file)
            networks_per_drug.append(network)
        all_networks.append(networks_per_drug)
    return all_networks


# if no drugs are specified, common drugs for each disease pair are determined and the function
# get_given_drugs_related_info is called with the parameters disease pairs and an array of drug IDs for each pair
def get_drugs_related_info(disease_pairs):
    networks = []
    for i in disease_pairs:
        networks.append(Network())
    networks = get_common_drugs(disease_pairs, networks, True)
    drugs = []
    for network in networks:
        pair_drugs = network.get_nodes_by_label('Drug')
        pair_drugs_ids = []
        for d in pair_drugs:
            pair_drugs_ids.append(d.id)
        drugs.append(pair_drugs_ids)
    networks = get_given_drugs_related_info(disease_pairs, drugs)
    return networks


# read a json file and visualise the network as PNG file
def draw_drug_subgraph(json_file):
    figure_file = json_file.replace('json', 'png')
    with open(json_file) as f:
        plt.figure(figsize=(18, 8))
        js_graph = json.load(f)
        G = nx.Graph()
        node_color = []
        for n in js_graph['nodes']:
            node_id = n['_id']
            label = n['_label']
            if label == 'Gene':
                node_color.append('red')
            elif label == 'RNA':
                node_color.append('gold')
            elif label == 'Variant':
                node_color.append('cornflowerblue')
            elif label == 'Disease':
                node_color.append('darkorange')
            elif label == 'Drug':
                node_color.append('grey')
            else:
                node_color.append('mediumseagreen')
            G.add_node(node_id)
        for n in js_graph['edges']:
            source_id = n['_source_id']
            target_id = n['_target_id']
            label = n['_label']
            G.add_edge(source_id, target_id, label = label)
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, font_weight='bold', node_color=node_color)
    for edge in G.edges:
        label = G.get_edge_data(edge[0], edge[1])
        label = label['label']
        nx.draw_networkx_edge_labels(G, pos, edge_labels={(edge[0], edge[1]): label}, font_color='black', label_pos=0.5)
    disease_patch = mpatches.Patch(color='darkorange', label='Disease')
    drug_patch = mpatches.Patch(color='grey', label='Drug')
    gene_patch = mpatches.Patch(color='red', label='Gene')
    rna_patch = mpatches.Patch(color='gold', label='RNA')
    variant_patch = mpatches.Patch(color='cornflowerblue', label='Variant')
    other_patch = mpatches.Patch(color='mediumseagreen', label='Other')
    plt.legend(handles=[disease_patch, drug_patch, gene_patch, rna_patch, variant_patch, other_patch])
    plt.savefig(figure_file)
   # plt.show()


# networks should be an array of networks
# this function should be used to analyse all drug subgraphs of a disease pair
# the most prevalent genes, RNAs and variants are printed to the console
# results can vary from run to run due to random selection of nodes in the function get_given_drugs_related_info if
# there are too many nodes to visualize all
def analysis_across_networks(networks):
    gene_count = {}
    rna_count = {}
    variant_count = {}
    for network in networks:
        genes = network.get_nodes_by_label('Gene')
        for gene in genes:
            gene_id = gene.id
            if gene_id in gene_count:
                gene_count[gene_id] += 1
            else:
                gene_count[gene_id] = 1

        rnas = network.get_nodes_by_label('RNA')
        edges = network.get_edges_by_label('REGULATES')
        for rna in rnas:
            rna_id = rna.id
           # edges = network.get_node_edges_by_label(rna, 'REGULATES')
            target_node_ids = set()
            for edge in edges:
                if edge.source_node_id == rna_id:
                    target_node_ids.add(edge.target_node_id)
            if rna_id in rna_count:
                all_genes = rna_count[rna_id]
                all_genes = all_genes.union(target_node_ids)
                rna_count[rna_id] = all_genes
            else:
                rna_count[rna_id] = target_node_ids

        variants = network.get_nodes_by_label('Variant')
        for variant in variants:
            variant_id = variant.id
            if variant_id in variant_count:
                variant_count[variant_id] += 1
            else:
                variant_count[variant_id] = 1

    i = 0
    for key, value in sorted(gene_count.items(), key=lambda item: item[1], reverse=True):
        if i < 10:
            print(key, value)
        else:
            break
        i += 1

    i = 0
    for key, value in sorted(rna_count.items(), key=lambda item: len(item[1]), reverse=True):
        if i < 10:
            print(key, len(value), value)
        else:
            break
        i += 1

    i = 0
    for key, value in sorted(variant_count.items(), key=lambda item: item[1], reverse=True):
        if i < 10:
            print(key, value)
        else:
            break
        i += 1


# creates the folders analysis and analysis/disease_pairs,
# IDs of all disease pairs can be returned
# all disease pairs or specified disease pairs can be analysed: a network per disease pair is returned and
# all regarding information can be written to files and json graphs
# subgraphs for each disease pair with specified drugs or all common drugs can be returned as network and json graph
# and a summary file is written
if __name__ == '__main__':
    driver = GraphDatabase.driver('bolt://localhost', auth=basic_auth('neo4j', '123'))
    session = driver.session()
    try:
        os.mkdir('../analysis')
    except FileExistsError:
        pass
    try:
        os.mkdir('../analysis/disease_pairs')
    except FileExistsError:
        pass

   # disease_pairs = get_disease_pairs(False)  # 12.916 Diseases = 83.405.070 pairs
    disease_pairs = [
        ['UMLS:C0011847', 'UMLS:C0085129'],  #  diabetes, asthma DrugBank:DB01223
     #  ['UMLS:C0085580', 'UMLS:C0085129'], # bluthochdruck, asthma

    ]
  #  networks = get_disease_pairs_info(disease_pairs, True)  #all infos
    all_drug_networks = get_drugs_related_info(disease_pairs)

    # or
    drugs = [
       # ['DrugBank:DB01223', 'DrugBank:DB01223'],
      #  ['DrugBank:DB00999']
    ['DrugBank:DB01136', 'DrugBank:DB01174']
    ]
   # all_drug_networks = get_given_drugs_related_info(disease_pairs, drugs)    #first pair with first drug list,
    pair_networks = all_drug_networks[0]
    analysis_across_networks(pair_networks)
