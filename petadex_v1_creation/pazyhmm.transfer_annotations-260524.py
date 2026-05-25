"""
-------------------------------------------------------
Transfer catalytic annotations from sequence to alignment coordinates.
Author: Thomas Quigley
-------------------------------------------------------
Objective: Catalytic residue annotations are stored as sequence positions
(counting only [A-Za-z] characters). HMMalign introduces gap characters ('-')
that shift these positions in the alignment. This script reads a Stockholm
alignment and an annotation CSV, then outputs a new CSV where each sequence
position has been converted to its corresponding alignment column number.

The alignment column number counts ALL characters (letters + '-') from the
start of the alignment, 1-based.

-------------------------------------------------------
Usage:
    python transfer_annotations.py <nr_aligns> <reference_annotations> <component> <nr_annotations_path>

Arguments:
    nr_aligns              : Path to the HMMalign Stockholm output (.sto / .aln)
    reference_annotations  : Path to the input annotation CSV. Format:
                               orf_id,component,aa_1,column_1,aa_2,column_2,aa_3,column_3,aa_4,column_4
                               2544106,c1,S,78,,,,,,
                             column_N values are 1-based sequence positions.
                             orf_id is informational only and is not used for matching.
    component              : Component name to process (e.g. c1, c2, c27)
    nr_annotations_path    : Path to write the output annotation CSV. Format:
                               component,aa_1,column_1,aa_2,column_2,aa_3,column_3,aa_4,column_4
                               c1,S,78,,,,,,
                             column_N values in the output are 1-based alignment columns.
                             When multiple processes write to the same file concurrently,
                             each process appends its row under an exclusive file lock.

Example:
    python transfer_annotations.py \
        /path/to/hmmalign_output.sto \
        /path/to/reference_annotations.csv \
        c1 \
        /path/to/nr_annotations.csv
-------------------------------------------------------
"""
# ─── PACKAGES ─── #
import sys
import os
import re
import csv
import fcntl


# ─── CONFIGURATION ─── #
MAX_CATALYTIC_SITES = 4


def safe_makedirs(path: str):
    """Creates all parent directories for a file path if they do not exist."""
    dir_path = os.path.dirname(os.path.abspath(path))
    os.makedirs(dir_path, exist_ok=True)


def parse_stockholm(stockholm_path: str) -> dict[str, str]:
    """
    Parses a Stockholm alignment file and returns a dict mapping
    seq_id -> full concatenated alignment row (including '-' gap characters).

    Skips all annotation lines (#=GR, #=GC, #=GF, #=GS) and '//' terminators.
    Handles multi-block Stockholm files by concatenating chunks per seq_id.
    """
    sequences: dict[str, str] = {}

    with open(stockholm_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")

            # Skip blank lines, header, terminators, and annotation lines
            if not line.strip():
                continue
            if line.startswith("#") or line.startswith("//"):
                continue

            fields = line.split(None, 1)
            if len(fields) != 2:
                continue

            seq_id, seq_chunk = fields
            sequences[seq_id] = sequences.get(seq_id, "") + seq_chunk.strip()

    return sequences


def seq_pos_to_aln_col(aligned_seq: str, seq_pos: int) -> tuple[int, str]:
    """
    Converts a 1-based sequence position to its HMM column number and returns
    the residue character found at that position.

    Sequence positions count all letter characters [A-Za-z] — both uppercase
    match-state residues and lowercase insertion residues. Gap characters
    ('-' and '.') are not counted.

    The HMM column is derived by:
        1. Finding the alignment column of the target residue (counting all
           letters, upper and lower, as sequence positions).
        2. Counting how many lowercase insertion letters appear up to and
           including that alignment column.
        3. Subtracting that insertion count from the alignment column number.

    This removes insertion columns from the coordinate space, giving a position
    that refers only to HMM match-state columns.

    Parameters
    ----------
    aligned_seq : full alignment row string (letters, '-', and '.' characters)
    seq_pos     : 1-based position in the sequence (counting all letters)

    Returns
    -------
    tuple[int, str] : (1-based HMM column number, residue character at seq_pos)

    Raises
    ------
    ValueError : if seq_pos exceeds the number of residues in aligned_seq
    """
    residue_count   = 0
    insertion_count = 0

    for aln_col, char in enumerate(aligned_seq, start=1):
        if char.isalpha():
            if char.islower():
                insertion_count += 1
            residue_count += 1
            if residue_count == seq_pos:
                hmm_col = aln_col - insertion_count
                return hmm_col, char

    raise ValueError(
        f"Sequence position {seq_pos} exceeds the number of residues "
        f"({residue_count}) in the aligned sequence."
    )


def read_reference_annotations(annotation_path: str, component: str) -> list[tuple[str, int]]:
    """
    Reads the reference annotation CSV and returns the catalytic sites for
    the given component as a list of (amino_acid, seq_position) tuples.

    The orf_id column is present for user readability and is ignored here.
    Only the first matching row for the component is used.

    CSV format:
        orf_id,component,aa_1,column_1,aa_2,column_2,aa_3,column_3,aa_4,column_4
    """
    catalytic_sites: list[tuple[str, int]] = []

    with open(annotation_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["component"] != component:
                continue
            for i in range(1, MAX_CATALYTIC_SITES + 1):
                aa  = row.get(f"aa_{i}",     "").strip()
                col = row.get(f"column_{i}", "").strip()
                if aa and col:
                    catalytic_sites.append((aa, int(col)))
            break  # one row per component

    return catalytic_sites


def write_output_annotations(
    output_path: str,
    component: str,
    catalytic_sites: list[tuple[str, int]],
):
    """
    Appends the output annotation row to a shared CSV file using an exclusive
    file lock, making it safe for multiple processes to write concurrently.

    If the file does not yet exist or is empty at the time this process
    acquires the lock, the CSV header is written first. Otherwise only the
    data row is appended.

    Output format:
        component,aa_1,column_1,aa_2,column_2,aa_3,column_3,aa_4,column_4
    """
    safe_makedirs(output_path)

    # Build a flat row padded to MAX_CATALYTIC_SITES entries
    row: dict[str, str] = {"component": component}
    for i in range(1, MAX_CATALYTIC_SITES + 1):
        if i <= len(catalytic_sites):
            aa, col = catalytic_sites[i - 1]
            row[f"aa_{i}"]     = aa
            row[f"column_{i}"] = str(col)
        else:
            row[f"aa_{i}"]     = ""
            row[f"column_{i}"] = ""

    fieldnames = ["component"] + [
        field
        for i in range(1, MAX_CATALYTIC_SITES + 1)
        for field in (f"aa_{i}", f"column_{i}")
    ]

    # Open in append mode so existing rows are never overwritten.
    # The exclusive lock (LOCK_EX) ensures only one process writes at a time;
    # all others block at flock() until the lock is released.
    with open(output_path, "a", newline="") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            # Write the header only if the file is still empty after we hold
            # the lock (guards against a race where two processes both see an
            # empty file before either has written the header).
            write_header = f.tell() == 0
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            if write_header:
                writer.writeheader()
            writer.writerow(row)
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def main():
    # ─── VALIDATE ARGUMENTS ─── #
    if len(sys.argv) != 5:
        print(
            "Usage: python transfer_annotations.py "
            "<nr_aligns> <reference_annotations> <component> <nr_annotations_path>"
        )
        print("See module docstring for full argument descriptions.")
        sys.exit(1)

    nr_aligns             = sys.argv[1]
    reference_annotations = sys.argv[2]
    component             = sys.argv[3]
    nr_annotations_path   = sys.argv[4]

    # ─── READ REFERENCE ANNOTATIONS ─── #
    print(f"[{component}] Reading reference annotations from {reference_annotations}...")
    catalytic_sites = read_reference_annotations(reference_annotations, component)

    if not catalytic_sites:
        print(f"[{component}] ERROR: No catalytic site annotations found for component '{component}'.")
        sys.exit(1)

    print(f"[{component}] Found {len(catalytic_sites)} catalytic site(s): "
          f"{', '.join(f'{aa}{col}' for aa, col in catalytic_sites)}")

    # ─── PARSE STOCKHOLM ALIGNMENT ─── #
    print(f"[{component}] Parsing Stockholm alignment from {nr_aligns}...")
    sequences = parse_stockholm(nr_aligns)

    if not sequences:
        print(f"[{component}] ERROR: No sequences found in {nr_aligns}.")
        sys.exit(1)

    print(f"[{component}] {len(sequences)} sequence(s) found in alignment.")

    # ─── USE THE FIRST SEQUENCE TO MAP COORDINATES ─── #
    # All sequences share the same alignment column structure, so any one of
    # them can be used to perform the sequence->alignment coordinate mapping.
    ref_seq_id, ref_aligned_seq = next(iter(sequences.items()))
    print(f"[{component}] Using '{ref_seq_id}' as the reference sequence for coordinate mapping.")

    # ─── CONVERT SEQUENCE POSITIONS TO ALIGNMENT COLUMNS ─── #
    print(f"[{component}] Converting sequence positions to alignment columns...")
    converted_sites: list[tuple[str, int]] = []

    for aa, seq_pos in catalytic_sites:
        try:
            aln_col, found_aa = seq_pos_to_aln_col(ref_aligned_seq, seq_pos)
            if found_aa.upper() == aa.upper():
                print(f"[{component}]   {aa}{seq_pos} (seq) -> HMM column {aln_col} "
                      f"[OK: residue '{found_aa}' matches expected '{aa}']")
            else:
                print(f"[{component}] ERROR: residue mismatch at seq pos {seq_pos} "
                      f"(expected '{aa}', found '{found_aa}') -> HMM column {aln_col}")
            converted_sites.append((aa, aln_col))
        except ValueError as e:
            print(f"[{component}] ERROR: {e}")
            sys.exit(1)

    # ─── WRITE OUTPUT ANNOTATIONS ─── #
    print(f"[{component}] Writing output annotations to {nr_annotations_path}...")
    write_output_annotations(nr_annotations_path, component, converted_sites)
    print(f"[{component}] Done.")


if __name__ == "__main__":
    main()