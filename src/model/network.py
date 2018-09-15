from model.node import Node
from model.edge import Edge


class Network:
    def __init__(self):
        self.nodes: {str: Node} = {}
        self.edges: [Edge] = []

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

    def get_nodes_by_label(self, label: str) -> [Node]:
        result = []
        for node in self.nodes.values():
            if node.label == label:
                result.append(node)
        return result

    def add_edge(self, edge: Edge):
        self.edges.append(edge)

    def get_edges_by_label(self, label: str) -> [Edge]:
        result = []
        for edge in self.edges:
            if edge.label == label:
                result.append(edge)
        return result

    def prune(self):
        """
        Remove nodes that are not connected to drugs and therefore of no interest.
        """
        pass
        ''' TODO
        # Remove genes of no interest
        removed_gene_ids = set()
        targeted_genes_id = {x[1] for x in self.get_edges_by_label('TARGETS')}
        for gene in set(self.get_nodes_by_label('Gene')):
            if targeted_genes_id.isdisjoint(gene.ids):
                for gene_id in gene.ids:
                    del self.nodes[gene_id]
                    removed_gene_ids.add(gene_id)
        for i in range(len(self.gene_associates_with_disease) - 1, -1, -1):
            if self.gene_associates_with_disease[i][0] in removed_gene_ids:
                del self.gene_associates_with_disease[i]
        for i in range(len(self.gene_codes_variant) - 1, -1, -1):
            if self.gene_codes_variant[i][0] in removed_gene_ids:
                del self.gene_codes_variant[i]
        # Remove variants of no interest
        removed_variant_ids = set()
        coded_variants_id = {x[1] for x in self.get_edges_by_label('CODES')}
        for variant in set(self.get_nodes_by_label('Variant')):
            if coded_variants_id.isdisjoint(variant.ids):
                for variant_id in variant.ids:
                    del self.nodes[variant_id]
                    removed_variant_ids.add(variant_id)
        for i in range(len(self.variant_associates_with_disease) - 1, -1, -1):
            if self.variant_associates_with_disease[i][0] in removed_variant_ids:
                del self.variant_associates_with_disease[i]
        '''

    def to_dict(self) -> {}:
        result = {
            'nodes': [],
            'edges': []
        }
        for node in set(self.nodes.values()):
            n = {'ids': sorted(node.ids), 'names': sorted(node.names), '_id': node.id, '_label': node.label}
            result['nodes'].append(n)
        for edge in self.edges:
            e = {'_label': edge.label, '_source': edge.source, '_target': edge.target}
            for key in edge.attributes:
                e[key] = edge.attributes[key]
            result['edges'].append(e)
        return result
