import io
import json

with io.open('../data/MONDO/lookup.json', 'r', encoding='utf-8', newline='') as f:
    lookup = json.load(f)
with io.open('../data/MONDO/reverse_lookup.json', 'r', encoding='utf-8', newline='') as f:
    reverse_lookup = json.load(f)


def map_from_to(source_id: str, target_prefix: str) -> str or None:
    if source_id not in reverse_lookup:
        return None
    for mondo_id in reverse_lookup[source_id]:
        for reference in lookup[mondo_id]['refs']:
            if reference.startswith(target_prefix):
                return reference
    return None


def map_from(source_id: str) -> [str]:
    result = set()
    if source_id in reverse_lookup:
        for mondo_id in reverse_lookup[source_id]:
            result.update(lookup[mondo_id]['refs'])
    return sorted(result)
