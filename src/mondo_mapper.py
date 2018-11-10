import io
import json
from typing import List

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


def map_from(source_id: str) -> (List[str], List[str]):
    result_ids = set()
    result_names = set()
    if source_id in reverse_lookup:
        for mondo_id in reverse_lookup[source_id]:
            result_ids.update(lookup[mondo_id]['refs'])
            result_names.add(lookup[mondo_id]['label'])
    return sorted(result_ids), result_names
