#!/usr/bin/env python3
# generate_fasta.py
# Copyright (c) 2026 Eliška Korbová ORCID 0009-0004-1247-0808
#
# Creating ddRAD library from filtered.csv
# Usage: python3 generate_fasta.py filtered.csv reference.fasta fragments.fasta
#
# Test results:
# python3 generate_fasta.py /results/HindIII_BfaI/filtered.csv zmays.fna /results/HindIII_BfaI/fragments_zm_HindIII_BfaI.fasta
# tested results (filtered.csv has header = - 1 line)
# grep -c "^>" fragments_zm_HindIII_BfaI.fasta 278320
# wc -l filtered.csv 278321 filtered.csv

from Bio import SeqIO
import csv
import sys

def load_genome(fasta_path):
    genomes = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        genomes[record.id] = record.seq
    return genomes


def make_fragments(csv_path, fasta_path, output_path):
    genomes = load_genome(fasta_path)

    with open(csv_path, newline="") as csvfile, open(output_path, "w") as out_fasta:
        reader = csv.DictReader(csvfile)

        for row in reader:
            acc = row["accession"]
            start = int(row["start_pos"])          # 0-based
            length = int(row["fragment_length"])

            start_enzyme = row["start_enzyme"]
            end_enzyme = row["end_enzyme"]

            if acc not in genomes:
                print(f"WARNING: {acc} not found in FASTA")
                continue

            seq = genomes[acc]

            # 0-based slicing (NO -1 shift)
            fragment = seq[start : start + length]

            header = f">{acc}_{start}_{start_enzyme}_{end_enzyme}_{length}"

            out_fasta.write(header + "\n")
            out_fasta.write(str(fragment) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python make_fragments_fasta.py input.csv ref.fasta output.fasta")
        sys.exit(1)

    csv_path = sys.argv[1]
    fasta_path = sys.argv[2]
    output_path = sys.argv[3]

    make_fragments(csv_path, fasta_path, output_path)