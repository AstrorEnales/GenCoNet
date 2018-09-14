#!/usr/bin/env python3

import io
import csv

omim_umls_map = {}

with io.open('../data/mappings/omim_to_umls.csv', 'r', encoding='utf-8', newline='') as f:
    reader = csv.reader(f, delimiter=',', quotechar='"')
    next(reader, None)
    for row in reader:
        omim_umls_map[row[0]] = row[1:]


def map_omim(omim_id: str) -> [str, str]:
    return omim_umls_map[omim_id] if omim_id in omim_umls_map else [None, None]


def save_mappings():
    with io.open('../data/mappings/omim_to_umls.csv', 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow(['OMIM', 'UMLS', 'UMLS name'])
        for key in sorted(omim_umls_map.keys()):
            writer.writerow([key] + omim_umls_map[key])
