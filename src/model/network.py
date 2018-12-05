from model.disease import Disease
from model.drug import Drug
from model.gene import Gene
from model.node import Node
from model.edge import Edge
from model.variant import Variant
from typing import List, Dict, Iterator


class Network:
    def __init__(self):
        self.nodes: Dict[str, Node] = {}
        self.edges: Dict[int, Edge] = {}
        self.edge_lookup: Dict[str, Dict[int, Edge]] = {}
        self.edge_source_lookup: Dict[str, Dict[int, Edge]] = {}
        self.edge_target_lookup: Dict[str, Dict[int, Edge]] = {}

    def add_node(self, node: Node):
        matches = [self.nodes[x] for x in node.ids if x in self.nodes]
        if any([node.label != x.label for x in matches]):
            print('[WARN] Label mismatch for id overlap:', node, [str(x) for x in matches])
        for match in matches:
            node.merge(match)
        for x in node.ids:
            self.nodes[x] = node

    def get_node_by_id(self, _id: str) -> Node:
        return self.nodes[_id] if _id in self.nodes else None

    def get_nodes(self) -> Iterator[Node]:
        for node in set(self.nodes.values()):
            yield node

    def get_nodes_by_label(self, label: str) -> List[Node]:
        result = set()
        for node in self.nodes.values():
            if node.label == label:
                result.add(node)
        return list(result)

    def add_edge(self, edge: Edge):
        self.edges[edge.id] = edge
        if edge.label not in self.edge_lookup:
            self.edge_lookup[edge.label] = {}
        self.edge_lookup[edge.label][edge.id] = edge
        if edge.source not in self.edge_source_lookup:
            self.edge_source_lookup[edge.source] = {}
        self.edge_source_lookup[edge.source][edge.id] = edge
        if edge.target not in self.edge_target_lookup:
            self.edge_target_lookup[edge.target] = {}
        self.edge_target_lookup[edge.target][edge.id] = edge

    def get_edges_by_label(self, label: str) -> List[Edge]:
        return list(self.edge_lookup[label].values()) if label in self.edge_lookup else []

    def get_node_edges_by_label(self, node: Node, label: str) -> List[Edge]:
        result = []
        if label in self.edge_lookup:
            for _id in node.ids:
                if _id in self.edge_source_lookup:
                    result.extend([e for e in self.edge_source_lookup[_id].values() if e.label == label])
                if _id in self.edge_target_lookup:
                    result.extend([e for e in self.edge_target_lookup[_id].values() if e.label == label])
        return result

    def get_edges_from_to(self, node_from: Node, node_to: Node, label: str) -> List[Edge]:
        result = []
        if label in self.edge_lookup:
            for _id in node_from.ids:
                if _id in self.edge_source_lookup:
                    result.extend([e for e in self.edge_source_lookup[_id].values()
                                   if e.label == label and e.target in node_to.ids])
        return result

    def delete_node(self, node: Node):
        edges = []
        for _id in node.ids:
            del self.nodes[_id]
            if _id in self.edge_source_lookup:
                edges.extend(self.edge_source_lookup[_id].values())
                del self.edge_source_lookup[_id]
            if _id in self.edge_target_lookup:
                edges.extend(self.edge_target_lookup[_id].values())
                del self.edge_target_lookup[_id]
        for edge in edges:
            if edge.id in self.edges:
                del self.edges[edge.id]
            if edge.id in self.edge_lookup[edge.label]:
                del self.edge_lookup[edge.label][edge.id]

    def prune(self):
        '''
        # Remove genes of no interest
        targeted_genes_id = {x.target for x in self.edge_lookup['TARGETS'].values()}
        for gene in set(self.get_nodes_by_label('Gene')):
            if targeted_genes_id.isdisjoint(gene.ids):
                self.delete_node(gene)
        # Remove variants of no interest
        coded_variants_id = {x.target for x in self.edge_lookup['CODES'].values()}
        for variant in set(self.get_nodes_by_label('Variant')):
            if coded_variants_id.isdisjoint(variant.ids):
                self.delete_node(variant)
        '''
        # Remove singletons
        for node in set(self.nodes.values()):
            if not any([_id in self.edge_source_lookup or _id in self.edge_target_lookup for _id in node.ids]):
                self.delete_node(node)

    def to_dict(self) -> {}:
        result = {
            'nodes': [],
            'edges': []
        }
        for node in set(self.nodes.values()):
            n = {'ids': sorted(node.ids), 'names': sorted(node.names), '_id': node.id, '_label': node.label}
            result['nodes'].append(n)
        for edge in self.edges.values():
            e = {'_label': edge.label, '_source': edge.source, '_target': edge.target}
            for key in edge.attributes:
                e[key] = edge.attributes[key]
            result['edges'].append(e)
        return result

    def load_from_dict(self, source: {}):
        for node in source['nodes']:
            if node['_label'] == 'Drug':
                self.add_node(Drug(node['ids'], node['names']))
            elif node['_label'] == 'Gene':
                self.add_node(Gene(node['ids'], node['names']))
            elif node['_label'] == 'Variant':
                self.add_node(Variant(node['ids'], node['names']))
            elif node['_label'] == 'Disease':
                self.add_node(Disease(node['ids'], node['names']))
        for edge in source['edges']:
            params = dict(edge)
            del params['_source']
            del params['_target']
            del params['_label']
            self.add_edge(Edge(edge['_source'], edge['_target'], edge['_label'], params))
