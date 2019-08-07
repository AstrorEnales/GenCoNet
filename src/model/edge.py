from typing import Dict, Any


class Edge:
    internal_id_counter = 1

    def __init__(self, source_node_id: str, target_node_id: str, label: str, attributes: Dict[str, Any]):
        self.source_node_id = source_node_id
        self.target_node_id = target_node_id
        self.label = label
        self.attributes = attributes
        self._id = Edge.internal_id_counter
        Edge.internal_id_counter += 1

    @property
    def id(self) -> int:
        return self._id

    def __str__(self) -> str:
        return 'Edge={label: %s, source: %s, target: %s}' % (self.label, self.source_node_id, self.target_node_id)
