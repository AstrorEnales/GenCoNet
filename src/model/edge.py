class Edge:
    internal_id = 1

    def __init__(self, source: str, target: str, label: str, attributes: {}):
        self.source = source
        self.target = target
        self.label = label
        self.attributes = attributes
        self._id = Edge.internal_id
        Edge.internal_id += 1

    @property
    def id(self) -> int:
        return self._id

    def __str__(self) -> str:
        return 'Edge={label: %s, source: %s, target: %s}' % (self.label, self.source, self.target)
