import io
from xml.dom import minidom
from xml.dom.minidom import Element, Document

from model.network import Network


class GraphMLExporter:
    def __init__(self, network: Network):
        self.network = network
        self.attribute_counter = 0

    def create_attribute_node(self, doc: Document, name: str, type: str, for_type: str) -> Element:
        n = doc.createElement('key')
        n.setAttribute('id', 'a%s' % self.attribute_counter)
        self.attribute_counter += 1
        n.setAttribute('for', for_type)
        n.setAttribute('attr.name', name)
        n.setAttribute('attr.type', type)
        return n

    @staticmethod
    def create_vertex_node(doc: Document, properties: dict) -> Element:
        n = doc.createElement('node')
        n.setAttribute('id', properties['id'])
        for key in properties:
            if key != 'id':
                data = doc.createElement('data')
                data.setAttribute('key', key)
                data.appendChild(doc.createTextNode(str(properties[key])))
                if key == 'label':
                    n.appendChild(data)
        return n

    @staticmethod
    def create_edge_node(doc: Document, source: str, target: str, label: str, properties: dict) -> Element:
        e = doc.createElement('edge')
        e.setAttribute('source', source)
        e.setAttribute('target', target)
        for key in properties:
            data = doc.createElement('data')
            data.setAttribute('key', key)
            data.appendChild(doc.createTextNode(str(properties[key])))
            #e.appendChild(data)
        data = doc.createElement('data')
        data.setAttribute('key', 'label')
        data.appendChild(doc.createTextNode(label))
        e.appendChild(data)
        return e

    def save(self, output_filepath):
        self.attribute_counter = 0
        doc = minidom.Document()
        root = doc.createElement('graphml')
        doc.appendChild(root)
        # Add attributes
        root.appendChild(self.create_attribute_node(doc, 'label', 'string', 'node'))
        root.appendChild(self.create_attribute_node(doc, 'ids', 'string', 'node'))
        root.appendChild(self.create_attribute_node(doc, 'names', 'string', 'node'))
        root.appendChild(self.create_attribute_node(doc, 'source', 'string', 'edge'))
        root.appendChild(self.create_attribute_node(doc, 'pmid', 'integer', 'edge'))
        root.appendChild(self.create_attribute_node(doc, 'score', 'string', 'edge'))
        root.appendChild(self.create_attribute_node(doc, 'num_pmids', 'integer', 'edge'))
        root.appendChild(self.create_attribute_node(doc, 'num_snps', 'integer', 'edge'))
        root.appendChild(self.create_attribute_node(doc, 'known_action', 'string', 'edge'))
        root.appendChild(self.create_attribute_node(doc, 'actions', 'string', 'edge'))
        root.appendChild(self.create_attribute_node(doc, 'simplified_action', 'string', 'edge'))
        # Add graph node
        graph_node = doc.createElement('graph')
        graph_node.setAttribute('id', 'GenCoNet')
        graph_node.setAttribute('edgedefault', 'directed')
        root.appendChild(graph_node)
        # Add nodes
        for gene in set(self.network.genes.values()):
            properties = {'label': 'Gene', 'id': gene.get_id(), 'ids': gene.ids, 'names': gene.names}
            graph_node.appendChild(self.create_vertex_node(doc, properties))
        for variant in set(self.network.variants.values()):
            properties = {'label': 'Variant', 'id': variant.get_id(), 'ids': variant.ids, 'names': variant.names}
            graph_node.appendChild(self.create_vertex_node(doc, properties))
        for drug in set(self.network.drugs.values()):
            properties = {'label': 'Drug', 'id': drug.get_id(), 'ids': drug.ids, 'names': drug.names}
            graph_node.appendChild(self.create_vertex_node(doc, properties))
        for disease in set(self.network.diseases.values()):
            properties = {'label': 'Disease', 'id': disease.get_id(), 'ids': disease.ids, 'names': disease.names}
            graph_node.appendChild(self.create_vertex_node(doc, properties))
        # Add edges
        for row in self.network.drug_indicates_disease:
            source = self.network.get_drug_by_id(row[0]).get_id()
            target = self.network.get_disease_by_id(row[1]).get_id()
            graph_node.appendChild(self.create_edge_node(doc, source, target, 'INDICATES', row[2]))
        for row in self.network.drug_contraindicates_disease:
            source = self.network.get_drug_by_id(row[0]).get_id()
            target = self.network.get_disease_by_id(row[1]).get_id()
            graph_node.appendChild(self.create_edge_node(doc, source, target, 'CONTRAINDICATES', row[2]))
        for row in self.network.drug_targets_gene:
            source = self.network.get_drug_by_id(row[0]).get_id()
            target = self.network.get_gene_by_id(row[1]).get_id()
            graph_node.appendChild(self.create_edge_node(doc, source, target, 'TARGETS', row[2]))
        for row in self.network.gene_associates_with_disease:
            source = self.network.get_gene_by_id(row[0]).get_id()
            target = self.network.get_disease_by_id(row[1]).get_id()
            graph_node.appendChild(self.create_edge_node(doc, source, target, 'ASSOCIATES_WITH', row[2]))
        for row in self.network.variant_associates_with_disease:
            source = self.network.get_variant_by_id(row[0]).get_id()
            target = self.network.get_disease_by_id(row[1]).get_id()
            graph_node.appendChild(self.create_edge_node(doc, source, target, 'ASSOCIATES_WITH', row[2]))
        for row in self.network.gene_codes_variant:
            source = self.network.get_gene_by_id(row[0]).get_id()
            target = self.network.get_variant_by_id(row[1]).get_id()
            graph_node.appendChild(self.create_edge_node(doc, source, target, 'CODES', row[2]))
        # Save graph to file
        with io.open(output_filepath, 'w', encoding='utf-8') as f:
            f.write(doc.toprettyxml(indent='  '))
