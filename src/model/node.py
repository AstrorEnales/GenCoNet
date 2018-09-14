class Node:
    def __init__(self, ids: [str], names: [str]):
        self.ids = set(ids)
        self.names = set(names)
        # Cache the hash value
        self.hash = ','.join(self.ids).__hash__()
        self.primary_id_prefix = ''

    def __str__(self) -> str:
        return '%s={ids: [%s], names: [%s]}' % (self.__class__.__name__, ','.join(sorted(self.ids)),
                                                ','.join(sorted(self.names)))

    def __eq__(self, o: object) -> bool:
        if o is not None and isinstance(o, type(self)):
            return len(self.ids.intersection(o.ids)) > 0
        return False

    def __hash__(self) -> int:
        return self.hash

    @property
    def label(self) -> str:
        return self.__class__.__name__

    @property
    def id(self) -> str or None:
        for x in self.ids:
            if x.startswith('%s:' % self.primary_id_prefix):
                return x
        return None

    def merge(self, o):
        self.ids.update(o.ids)
        self.names.update(o.names)
        # Update the hash value with new ids
        self.hash = ','.join(self.ids).__hash__()
