#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from tqdm import tqdm

# We need python-Levenshtein installed for the Levenshtein mode
import Levenshtein  # pip install python-Levenshtein

###############################################################################
# DEFAULTS
###############################################################################
DEFAULT_INPUT_FASTQ = 'Calls.fastq'
DEFAULT_OUTPUT_CSV = 'matched_results.csv'
DEFAULT_DESIGN_CSV = 'design.csv'

# Adapter trimming
DEFAULT_LIB_LENGTH = 64
DEFAULT_FRONT_PRIMER = "ACACGACGCTCTTCCGATCT" #"TAAGAGACAG"
DEFAULT_BACK_PRIMER = "AGATCGGAAGAGCACACGTCT" #"CTGTCTCTTA"
DEFAULT_LENGTH_TOLERANCE = 64

# Barcode matching
DEFAULT_BARCODE_START_IDX = 0
DEFAULT_BARCODE_END_IDX = 10
DEFAULT_MAX_DISTANCE = 2
DEFAULT_DISTANCE_MODE = 'levenshtein'  # or 'hamming'

###############################################################################
# Distance Functions
###############################################################################
def hamming_distance(s1, s2):
    """Compute Hamming distance between two strings of equal length."""
    if len(s1) != len(s2):
        return float('inf')  # treat length mismatch as an impossible match
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def compute_distance(s1, s2, mode='levenshtein'):
    """
    Compute distance between s1 and s2 using either 'levenshtein' or 'hamming'.
    """
    if mode == 'levenshtein':
        return Levenshtein.distance(s1, s2)
    elif mode == 'hamming':
        return hamming_distance(s1, s2)
    else:
        raise ValueError(f"Unknown distance mode: {mode}")

###############################################################################
# Reverse Complement
###############################################################################
def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rc = ''.join(complement.get(base, base) for base in seq)[::-1]
    return rc

###############################################################################
# Adapter/Primer Trimming
###############################################################################
def trim_read(seq, front_primer, back_primer, lib_length, length_tolerance):
    """
    Attempt to trim the read `seq` by finding front_primer and back_primer.
    Return the portion between them (or near them) if it has length of the lib_length
    (plus/minus the length tolerance), else return None.
    """

    seq = seq.strip()

    # Try forward orientation (in case the reads was obtained forwrd)
    front_idx = seq.find(front_primer)
    back_idx = seq.find(back_primer)
    if front_idx != -1 and back_idx != -1:
        insert_seq = seq[front_idx + len(front_primer) : back_idx]
        if abs(len(insert_seq) - lib_length) <= length_tolerance:
            return insert_seq
    elif front_idx != -1:
        start = front_idx + len(front_primer)
        end = min(start + lib_length, len(seq))
        insert_seq = seq[start:end]
        if abs(len(insert_seq) - lib_length) <= length_tolerance:
            return insert_seq
    elif back_idx != -1:
        end = back_idx
        start = max(end - lib_length, 0)
        insert_seq = seq[start:end]
        if abs(len(insert_seq) - lib_length) <= length_tolerance:
            return insert_seq

    # Try reverse orientation
    rc_seq = reverse_complement(seq)
    front_idx = rc_seq.find(front_primer)
    back_idx = rc_seq.find(back_primer)
    if front_idx != -1 and back_idx != -1:
        insert_seq = rc_seq[front_idx + len(front_primer) : back_idx]
        if abs(len(insert_seq) - lib_length) <= length_tolerance:
            return insert_seq
    elif front_idx != -1:
        start = front_idx + len(front_primer)
        end = min(start + lib_length, len(rc_seq))
        insert_seq = rc_seq[start:end]
        if abs(len(insert_seq) - lib_length) <= length_tolerance:
            return insert_seq
    elif back_idx != -1:
        end = back_idx
        start = max(end - lib_length, 0)
        insert_seq = rc_seq[start:end]
        if abs(len(insert_seq) - lib_length) <= length_tolerance:
            return insert_seq

    # If we reach here, no suitable trimming found
    return None

###############################################################################
# Load Designs
###############################################################################
def load_designs(design_csv_path):
    """
    Read barcodes+sequences from a CSV with columns [barcode, sequence] (SOLQC format).
    """
    designs = []
    with open(design_csv_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            barcode = row['barcode'].strip()
            sequence = row['sequence'].strip()
            designs.append((barcode, sequence))
    return designs

###############################################################################
# Best Match to a Design
###############################################################################
def best_match_variant_id(trimmed_seq, designs, start_idx, end_idx,
                          max_distance, distance_mode):
    """
    Extract the substring [start_idx:end_idx] from trimmed_seq as the barcode.
    Compare to each design's barcode using the chosen distance mode.
    Return variant_id if min distance <= max_distance, else -1.
    """
    barcode_region = trimmed_seq[start_idx:end_idx]
    best_dist = None
    best_idx = -1
    for i, (design_barcode, design_seq) in enumerate(designs):
        d = compute_distance(barcode_region, design_barcode, mode=distance_mode)
        if best_dist is None or d < best_dist:
            best_dist = d
            best_idx = i
    if best_dist is not None and best_dist <= max_distance:
        return best_idx
    else:
        return -1

###############################################################################
# Main
###############################################################################
def main():
    parser = argparse.ArgumentParser(
        description="Combine adapter trimming and barcode matching."
    )

    #  Arguments for input (fastq + design file), and output file
    parser.add_argument("--fastq", default=DEFAULT_INPUT_FASTQ,
                        help="Path to input FASTQ (default: %(default)s)")
    parser.add_argument("--design_csv", default=DEFAULT_DESIGN_CSV,
                        help="Path to design CSV (default: %(default)s)")
    parser.add_argument("--output_csv", default=DEFAULT_OUTPUT_CSV,
                        help="Path to output CSV (default: %(default)s)")

    # Arguments for adapter trimming, length tolerance means plus minus the length tolerance
    parser.add_argument("--lib_length", type=int, default=DEFAULT_LIB_LENGTH,
                        help="Expected library length (default: %(default)s)")
    parser.add_argument("--front_primer", default=DEFAULT_FRONT_PRIMER,
                        help="Front primer sequence (default: %(default)s)")
    parser.add_argument("--back_primer", default=DEFAULT_BACK_PRIMER,
                        help="Back primer sequence (default: %(default)s)")
    parser.add_argument("--length_tolerance", type=int,
                        default=DEFAULT_LENGTH_TOLERANCE,
                        help="Allowed deviation from lib_length (default: %(default)s)")

    # Barcode matching
    parser.add_argument("--barcode_start_idx", type=int,
                        default=DEFAULT_BARCODE_START_IDX,
                        help="Barcode region start (default: %(default)s)")
    parser.add_argument("--barcode_end_idx", type=int,
                        default=DEFAULT_BARCODE_END_IDX,
                        help="Barcode region end (default: %(default)s)")
    parser.add_argument("--max_distance", type=int,
                        default=DEFAULT_MAX_DISTANCE,
                        help="Max allowed distance for match (default: %(default)s)")
    parser.add_argument("--distance_mode", choices=["levenshtein", "hamming"],
                        default=DEFAULT_DISTANCE_MODE,
                        help="Distance mode (levenshtein|hamming) (default: %(default)s)")

    args = parser.parse_args()

    # 1) Load designs
    designs = load_designs(args.design_csv)
    # The index of each design is the variant_id.

    # 2) Read all FASTQ lines into memory for a single pass
    with open(args.fastq, 'r') as f:
        lines = f.read().splitlines()

    total_reads = len(lines) // 4

    # 3) Process reads in blocks of 4 lines (fastq format)
    read_counts = defaultdict(int)
    filtered_out_count = 0

    for i in tqdm(range(0, len(lines), 4), desc="Processing", total=total_reads):
        # Each read: ID, SEQ, +, QUAL
        if i+3>len(lines):
            break
        else:
            read_id   = lines[i]
            read_seq  = lines[i+1]
            plus_line = lines[i+2]
            read_qual = lines[i+3]

            # Step A: Adapter Trim
            trimmed = trim_read(
                seq=read_seq,
                front_primer=args.front_primer,
                back_primer=args.back_primer,
                lib_length=args.lib_length,
                length_tolerance=args.length_tolerance
            )

            if trimmed is None:
                filtered_out_count += 1
                continue

            # Step B: Barcode Matching
            variant_id = best_match_variant_id(
                trimmed_seq=trimmed,
                designs=designs,
                start_idx=args.barcode_start_idx,
                end_idx=args.barcode_end_idx,
                max_distance=args.max_distance,
                distance_mode=args.distance_mode
            )

            # If you want to skip unmatched barcodes, you could do:
            # if variant_id == -1:
            #     continue

            # Step C: Aggregate
            read_counts[(trimmed, variant_id)] += 1

    # 4) Write aggregated CSV
    with open(args.output_csv, 'w', newline='') as out_csv:
        writer = csv.writer(out_csv)
        writer.writerow(["enum", "sequence", "count", "variant_id"])
        for idx, ((seq, v_id), count) in enumerate(read_counts.items()):
            writer.writerow([idx, seq, count, v_id])

    print(f"Done! Processed {total_reads} reads.")
    print(f"Filtered out {filtered_out_count} reads (no matching primers).")
    print(f"Results written to {args.output_csv}.")

if __name__ == "__main__":
    main()
