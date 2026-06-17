#!/usr/bin/env python3
# combination_processing.py
# Copyright (c) 2026 Eliška Korbová ORCID 0009-0004-1247-0808
#
# Preparing motifs for Aho-Corasick and combinations for rest of the pipeline
# Output: enzymes_run.csv, combinations.csv
import os
import sys
import csv
import argparse

FAST_COMBOS = [
    ('EcoRI',    'MseI'),
    ('EcoRI',    'MboI'),
    ('EcoRI',    'NlaIII'),
    ('PstI',     'MspI'),
    ('PstI',     'MseI'),
    ('PstI',     'TaqI'),
    ('PstI',     'MboI'),
    ('PstI',     'HpaII'),
    ('SbfI',     'MspI'),
    ('SbfI',     'MseI'),
    ('ApeKI',    'MseI'),
    ('BamHI',    'MseI'),
    ('HindIII',  'MseI'),
    ('HindIII',  'NlaIII'),
    ('NheI',     'MboI'),
    ('NheI',     'MseI'),
    ('XbaI',     'MboI'),
    ('SpeI',     'MboI'),
    ('KpnI',     'MboI'),
    ('NcoI',     'MseI'),
    ('NcoI',     'MboI'),
    ('BglII',    'MseI'),
    ('ClaI',     'MboI'),
]

# loading enzymes.csv with expanded IUPAC characters
def load_enzymes_csv(path) -> dict[str, list[dict]]:
    if not os.path.isfile(path):
        sys.exit(f'[ERROR] enzymes.csv not found: {path}')
    db = {}
    with open(path, newline='') as fh:
        for row in csv.DictReader(fh):
            db.setdefault(row['enzyme_name'], []).append(row)
    return db

# combinations for rest of the pipeline
def load_custom(path) -> list[tuple[str, str]]:
    if not path:
        sys.exit('[ERROR] --combinations custom requires --combinations-file')

    if not os.path.isfile(path):
        sys.exit(f'[ERROR] combinations file not found: {path}')

    combos = []
    with open(path, newline='') as fh:
        reader = csv.reader(fh)
        next(reader, None)  # skip header
        for i, parts in enumerate(reader, 2):
            if len(parts) != 2:
                print(f'  [WARN] Line {i} bad format, skipping: {parts}')
                continue
            combos.append((parts[0].strip(), parts[1].strip()))
    return combos

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workdir', required=True)
    parser.add_argument('--combinations', default='fast', choices=['fast', 'custom'])
    parser.add_argument('--combinations-file', default=None)
    args = parser.parse_args()

    results_dir = os.path.join(args.workdir, 'results')
    os.makedirs(results_dir, exist_ok=True)

    script_dir = os.path.dirname(os.path.abspath(__file__))

    enzymes_csv = os.path.join(script_dir, 'enzymes.csv')
    allowed_path = os.path.join(script_dir, 'allowed_pairs.csv')

    print(f'[1/3] Loading: {enzymes_csv}')
    enzyme_db: dict[str, list[dict]] = load_enzymes_csv(enzymes_csv)
    print(f'      {len(enzyme_db)} enzymes in master CSV')

    if not os.path.isfile(allowed_path):
        sys.exit(f'[ERROR] allowed_pair.csv not found: {allowed_path}')

    allowed_set = set()
    with open(allowed_path, newline='') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            a = row['enzyme_a'].strip()
            b = row['enzyme_b'].strip()
            allowed_set.add(tuple(sorted((a, b))))

    # check - is enzyme in the enzymes.csv and is the pair allowed
    print(f'[2/3] Resolving combinations: {args.combinations}')
    raw_combos = FAST_COMBOS if args.combinations == 'fast' else load_custom(args.combinations_file)

    valid_combos: list[tuple[str, str]] = []  # both enzymes must be in enzyme_db and allowed_set
    for ea, eb in raw_combos:
        missing = [e for e in (ea, eb) if e not in enzyme_db]
        if missing:
            print(f'  [WARN] {missing} not in enzymes.csv, skipping {ea}+{eb}')
            continue

        key = tuple(sorted((ea, eb)))
        if key not in allowed_set:
            print(f'  [WARN] {ea}+{eb} was rejected, cannot simulate whole overlapping recognition sites')
            continue

        valid_combos.append((ea, eb))

    if not valid_combos:
        print('[ERROR] No valid enzyme combinations found.')
        sys.exit(1)
    print(f'      {len(valid_combos)} valid combinations')

    # unique enzymes for the run
    needed = sorted({e for pair in valid_combos for e in pair})
    print(f'      {len(needed)} unique enzymes needed')

    # create enzymes_run.csv for rust
    print('[3/3] Writing run files ...')
    enzymes_run_path = os.path.join(results_dir, 'enzymes_run.csv')
    rows = [row for name in needed for row in enzyme_db[name]]
    with open(enzymes_run_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['enzyme_name', 'expanded_sequence', 'cut_offset'])
        w.writeheader()
        w.writerows(rows)
    print(f'  -> enzymes_run.csv  : {len(rows)} rows ({len(needed)} enzymes)')

    # create combinations.csv for fragment_processing
    combos_path = os.path.join(results_dir, 'combinations.csv')
    with open(combos_path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['enzyme_a', 'enzyme_b'])
        w.writerows(valid_combos)
    print(f'  -> combinations.csv : {len(valid_combos)} pairs')

    print(f'\n[OK] {results_dir}')

if __name__ == '__main__':
    main()
