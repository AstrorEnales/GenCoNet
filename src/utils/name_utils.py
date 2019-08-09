import re
import difflib
from typing import Set


def normalize_node_names(names: Set[str]) -> Set[str]:
    if len(names) <= 1:
        return names
    similar_names = {}
    for name in names:
        name_lower = name.lower()
        if name_lower not in similar_names:
            similar_names[name_lower] = []
        similar_names[name_lower].append(name)
    result = set()
    for key in similar_names:
        max_upper_count = -1
        max_index = 0
        candidates = similar_names[key]
        for i in range(0, len(candidates)):
            upper_count = sum([1 if x[0].isupper() and not x.isupper() else 0 for x in candidates[i].split(' ')])
            if upper_count > max_upper_count:
                max_upper_count = upper_count
                max_index = i
        result.add(candidates[max_index])
    lower_results = {x.lower() for x in result}
    for name in list(result):
        if ', ' in name and ' '.join(name.split(', ')[::-1]).lower() in lower_results:
            result.remove(name)
    return result


def node_names_synonym(names: Set[str]) -> bool:
    names = sorted({x.lower() for x in names})
    if len(names) == 1:
        return True
    valid = True
    for i in range(0, len(names) - 1):
        a = names[i]
        for j in range(i + 1, len(names)):
            b = names[j]
            # Check for specific diff cases
            diff = list(difflib.ndiff(a, b))
            changes = [(i, x) for i, x in enumerate(diff) if x[0] != ' ']
            if len(changes) == 1:
                if changes[0][0] > 0 and changes[0][1] == '- h' and a[changes[0][0] - 1] == 't':
                    continue
                if changes[0][0] > 0 and changes[0][1] == '+ h' and a[changes[0][0] - 1] == 't':
                    continue
                if changes[0][0] == len(diff) - 1 and changes[0][1] == '+ e' and a[changes[0][0] - 1] == 'n':
                    continue
            elif len(changes) == 2:
                if changes[0][0] == changes[1][0] - 1 and changes[0][1] == '- i' and changes[1][1] == '+ y':
                    continue
            elif len(changes) == 3:
                if changes[0][0] == changes[1][0] - 1 and changes[0][0] == changes[2][0] - 2 and \
                        changes[0][1] == '- f' and changes[1][1] == '+ p' and changes[2][1] == '+ h':
                    continue
            # Check for different word order
            a_parts = sorted([x for x in re.split(r'[, ()]', a) if x])
            b_parts = sorted([x for x in re.split(r'[, ()]', b) if x])
            if ' '.join(a_parts) == ' '.join(b_parts):
                continue
            valid = False
    return valid
