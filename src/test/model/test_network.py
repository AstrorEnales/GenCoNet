import io
import os
import json
import unittest
from model import node, edge, network


class DummyNode1(node.Node):
    def __init__(self, ids: [str], names: [str]):
        super().__init__(ids, names)
        self.primary_id_prefix = 'TEST1'


class DummyNode2(node.Node):
    def __init__(self, ids: [str], names: [str]):
        super().__init__(ids, names)
        self.primary_id_prefix = 'OTHER'


class DummyNode3(node.Node):
    def __init__(self, ids: [str], names: [str]):
        super().__init__(ids, names)
        self.primary_id_prefix = 'TEST2'


class TestMethods(unittest.TestCase):
    def setUp(self):
        self.temp_file_path = 'genconet_testgraph.json'

    def test_label(self):
        n1 = DummyNode1(['TEST1:123', 'OTHER:abc'], [])
        n2 = DummyNode2(['OTHER:abc'], [])
        n3 = DummyNode3(['TEST2:456', 'SOME:rtf'], [])
        n4 = DummyNode3(['TEST2:456', 'SOME:rtf'], [])
        graph = network.Network()
        graph.add_node(n1)
        graph.add_node(n2)
        graph.add_node(n3)
        graph.add_node(n4)
        graph.add_edge(edge.Edge(n3, n1, 'LINKS', {}))
        graph.add_edge(edge.Edge(n3, n2, 'LINKS', {}))
        graph.add_edge(edge.Edge(n1, n2, 'LINKS', {}))
        self.assertEqual(len(list(graph.get_nodes())), 3)
        self.assertEqual(len(set(graph.nodes.values())), 3)
        self.assertEqual(len(list(graph.get_edges_by_label('LINKS'))), 3)
        graph.save(self.temp_file_path)
        graph = network.Network()
        with io.open(self.temp_file_path, 'r', encoding='utf-8', newline='') as f:
            g = json.loads(f.read())
            graph.load_from_dict(g)
        self.assertEqual(len(list(graph.get_nodes())), 3)
        self.assertEqual(len(list(graph.get_edges_by_label('LINKS'))), 3)

    def tearDown(self):
        os.remove(self.temp_file_path)
