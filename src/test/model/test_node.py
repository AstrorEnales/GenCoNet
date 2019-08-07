import unittest
from model import node


class DummyNode(node.Node):
    def __init__(self, ids: [str], names: [str]):
        super().__init__(ids, names)
        self.primary_id_prefix = 'TEST'


class TestMethods(unittest.TestCase):
    def test_label(self):
        n = DummyNode([], [])
        self.assertEqual(n.label, 'DummyNode')

    def test_str(self):
        n = DummyNode(['TEST:1'], ['test name'])
        self.assertEqual(str(n), 'DummyNode={ids: [TEST:1], names: ["test name"]}')
