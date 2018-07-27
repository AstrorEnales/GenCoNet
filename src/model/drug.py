class Drug:
    def __init__(self, ids: [str], names: [str]):
        self.ids = set(ids)
        self.names = set(names)
        # Cache the hash value
        self.hash = ','.join(self.ids).__hash__()

    def __str__(self) -> str:
        return 'Drug={ids: [%s], names: [%s]}' % (','.join(sorted(self.ids)), ','.join(sorted(self.names)))

    def __eq__(self, o: object) -> bool:
        if o is not None and isinstance(o, type(Drug)):
            return len(self.ids.intersection(o.ids)) > 0
        return False

    def __hash__(self) -> int:
        return self.hash

    def get_id(self) -> str:
        return [x for x in self.ids if x.startswith('DrugBank:')][0]

    def merge(self, drug):
        self.ids.update(drug.ids)
        self.names.update(drug.names)
        # Update the hash value with new ids
        self.hash = ','.join(self.ids).__hash__()
