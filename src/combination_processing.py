#!/usr/bin/env python3
"""
combination_processing.py  -  Prvni krok
Mod renda (fast / custom), kontroluje nazvy enzymu oproti masteru enzymes.csv
vypisuje do <workdir>/results/:
  enzymes_run.csv   - sekvence pouze pro tento beh sloupce: enzyme_name, expanded_sequence, cut_offset primo pro rust
  combinations.csv  - overene kombinace enzymu co existuji sloupce: enzyme_a, enzyme_b (pro rendogbs_pipeline.py)
Implementace (rendogbs.sh):
  python combination_processing.py \\
      --enzymes-csv  enzymes.csv \\
      --workdir      ./run1 \\
      --combinations fast \\
      [--combinations-file my_pairs.txt]
POZOR all zruseno, pak upravim help, to by bylo ohromne, neprehledne grafy a hlavne zbytecne
"""
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
    ('BamHI',    'MboI'),
    ('HindIII',  'MseI'),
    ('HindIII',  'NlaIII'),
    ('NheI',  'MboI'),
    ('NheI',  'MseI'),
    ('XbaI',     'MboI'),
    ('SpeI',  'MboI'),
    ('KpnI',  'MboI'),
    ('NcoI',     'MseI'),
    ('NcoI',     'MboI'),
    ('BglII',    'MboI'),
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

# combinations for rendogbs_pipeline
def load_custom(path) -> list[tuple[str, str]]:
    if not path:
        sys.exit('[ERROR] --combinations custom requires --combinations-file')
    combos = []
    with open(path, newline='') as fh:
        reader = csv.reader(fh)
        next(reader, None) # skip header
        for i, parts in enumerate(reader, 2):
            if len(parts) != 2:
                print(f'  [WARN] Line {i} bad format, skipping: {parts}') 
                continue
            combos.append((parts[0].strip(), parts[1].strip()))
    return combos

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--enzymes-csv', required=True)
    parser.add_argument('--workdir', required=True,)
    parser.add_argument('--combinations', default='fast', choices=['fast', 'custom'])
    parser.add_argument('--combinations-file', default=None)
    args = parser.parse_args()

    results_dir = os.path.join(args.workdir, 'results')
    os.makedirs(results_dir, exist_ok=True)

    print(f'[1/3] Loading: {args.enzymes_csv}')
    enzyme_db: dict[str, list[dict]] = load_enzymes_csv(args.enzymes_csv)
    print(f'      {len(enzyme_db)} enzymes in master CSV')

    # check - is enzyme in the enzymes.csv
    print(f'[2/3] Resolving combinations: {args.combinations}')
    raw_combos = (FAST_COMBOS if args.combinations == 'fast' else load_custom(args.combinations_file))

    valid_combos: list[tuple[str, str]] = [] # both enzymes must be in enzyme_db
    for ea, eb in raw_combos:
        missing = [e for e in (ea, eb) if e not in enzyme_db] # if any is missing discard the pair
        if missing:
            print(f'  [WARN] {missing} not in enzymes.csv, skipping {ea}+{eb}')
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

    # create combinations.csv for fragment_processing.py
    combos_path = os.path.join(results_dir, 'combinations.csv')
    with open(combos_path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['enzyme_a', 'enzyme_b'])
        w.writerows(valid_combos)
    print(f'  -> combinations.csv : {len(valid_combos)} pairs')
    print(f'\n[OK] {results_dir}')

if __name__ == '__main__':
    main()