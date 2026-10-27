"""
-------------------------------------------------------
Clean HMMsearch outputs.
Author: Thomas Quigley
-------------------------------------------------------
Objective: The output from HMMsearch is in two forms:
1. The DOMTBLOUT that has e-values and regions hit.
2. The alignments of HMM columns to sequence.

I will parse through the DOMTBLOUT to create a useable
parquet, and I will extract the .afa formatted
alignments for all HMM hits. These will then be used
to select sequences that have a good E-value and
catalytic sites to create the BlastNR HMM.
-------------------------------------------------------
Usage:
    python 260406_clean_hmmsearch.py <results_path> <component> <analysis_path> <annotated_path> <fragments_dir>

Arguments:
    results_path    : Path to the hmmsearch results directory. Expected structure:
                        results_path/
                            aln/        {component}.aln       - Stockholm alignment output from hmmsearch -A
                            domtbl/     {component}.domtbl    - Domain table output from hmmsearch --domtblout
                            e_value_csv/{component}.csv       - [CREATED] All hits with orf_id, hmm_range, i-Evalue
                            aligned_afa/{component}.afa       - [CREATED] Aligned sequences in FASTA format

    component       : Component name matching the HMM and file names (e.g. c1, c2, c27)

    analysis_path   : Path to the analysis output directory. Expected structure:
                        analysis_path/
                            catalytic/  {component}.csv       - [CREATED] orf_id, hmm_range, i_evalue, annotations

    annotated_path  : Path to the catalytic residue annotation CSV. Format:
                        component,aa_1,column_1,aa_2,column_2,aa_3,column_3,aa_4,column_4
                        c1,S,70,D,20,,,,
                        c2,A,30,,,,,,
                        c3,A,10,A,20,A,30,A,40
                      Columns refer to 1-based positions in the consensus (match-state) sequence.

    fragments_dir   : Path to directory where fragment alignment files are written.
                        fragments_dir/
                            {component}.afa   - [CREATED] Aligned sequences below e-value cutoff
                                               that do NOT contain all catalytic residues.

Example:
    python 260406_clean_hmmsearch.py \
        /watson/rnalab/tquigley/projects/logan_paper_revised/redo_hmms/pazy_approach/hmm_results/pazy_search \
        c1 \
        /watson/rnalab/tquigley/projects/logan_paper_revised/redo_hmms/pazy_approach/catalytic_hits/pazy_search \
        /watson/rnalab/tquigley/projects/logan_paper_revised/redo_hmms/pazy_approach/catalytic_annotations.csv \
        /watson/rnalab/tquigley/projects/logan_paper_revised/redo_hmms/pazy_approach/fragments/blastnr
-------------------------------------------------------
"""
# ─── PACKAGES ─── #
import sys
import os
import re
import csv
import pyarrow.parquet as pq
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pa_csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

# ─── CONFIGURATION ─── #
E_VALUE_CUTOFF = 1.0e-5
COVERAGE_CUTOFF = 0.5

# ─── ENSURE PRINT IS HAPPENING ─── #
print("DEBUG: Python script started", flush=True)
print(f"DEBUG: Arguments: {sys.argv}", flush=True)

def safe_makedirs(path: str):
    """Creates all parent directories for a file path if they do not exist."""
    dir_path = os.path.dirname(os.path.abspath(path))
    os.makedirs(dir_path, exist_ok=True)


def parse_domtblout(domtbl_path: str, e_value_csv_path: str) -> pa.Table:
    """
    Parses hmmsearch domtblout into a pyarrow Table.

    Extracts:
        orf_id    : integer before the first '|' in the target name
        hmm_range : 'ali_from-ali_to' built from domtblout fields 17-18
                    (the alignment envelope coordinates, matching the /start-end
                    suffix written into Stockholm seq IDs by hmmsearch -A)
        i_evalue  : field index 12 (domain i-E-value)

    Every individual hit row is kept as-is — no deduplication.
    Writes all results to CSV and returns only rows at or below E_VALUE_CUTOFF.
    """
    orf_ids:    list[int]   = []
    hmm_ranges: list[str]   = []
    i_evalues:  list[float] = []

    with open(domtbl_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue

            fields = line.split()
            if len(fields) < 13:
                continue

            target_name = fields[0]

            # ─── EXTRACT ORF_ID ─── #
            orf_id = int(target_name.split("|")[0])

            # ─── EXTRACT HMM_RANGE ─── #
            # The domtblout does not carry a /start-end suffix in the target name.
            # The alignment envelope coordinates are at fields[17] (ali_from) and
            # fields[18] (ali_to), which match the range written into the Stockholm
            # seq IDs by hmmsearch -A.
            if len(fields) < 19:
                continue
            hmm_range = f"{fields[17]}-{fields[18]}"

            i_evalue = float(fields[12])

            orf_ids.append(orf_id)
            hmm_ranges.append(hmm_range)
            i_evalues.append(i_evalue)

    # ─── BUILD FULL TABLE ─── #
    table = pa.table({
        "orf_id":    pa.array(orf_ids,    type=pa.int64()),
        "hmm_range": pa.array(hmm_ranges, type=pa.string()),
        "i_evalue":  pa.array(i_evalues,  type=pa.float64()),
    })

    # ─── WRITE FULL CSV ─── #
    safe_makedirs(e_value_csv_path)
    pa_csv.write_csv(table, e_value_csv_path)

        # ─── FILTER TO E_VALUE_CUTOFF ─── #
    mask = pc.less_equal(table.column("i_evalue"), E_VALUE_CUTOFF)
    filtered_table = table.filter(mask)

    return filtered_table


def extracted_aligned_sequences(aln_path: str, aligned_afa_path: str) -> pa.Table:
    """
    Parses hmmsearch Stockholm alignment output into a pyarrow Table.

    Extracts:
        orf_id      : integer before the first '|' in the sequence ID
        hmm_range   : the 'start-end' string after the '/' in the sequence ID
        aligned_seq : concatenated aligned sequence chunks

    All rows are kept in the returned table and written to the .afa file.
    Sequences are written with their full alignment including insertion columns
    (lowercase letters and dots are preserved).
    """
    # Each entry is a separate hit row: list of (orf_id, hmm_range, seq_chunks)
    # We need to track per-row seq chunks because Stockholm may split long
    # sequences across multiple lines with the same seq_id.
    # Key: (orf_id, hmm_range) is NOT unique (same orf could hit at same range
    # in a wrapped Stockholm block), so we accumulate chunks per seq_id string.
    seq_chunks: dict[str, str] = {}  # seq_id -> accumulated seq
    row_order:  list[str]     = []  # seq_id insertion order (one entry per unique seq_id)

    with open(aln_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith("#") or line.startswith("//") or not line.strip():
                continue

            fields = line.split(None, 1)
            if len(fields) != 2:
                continue

            seq_id, seq_chunk = fields

            if seq_id not in seq_chunks:
                seq_chunks[seq_id] = ""
                row_order.append(seq_id)

            seq_chunks[seq_id] += seq_chunk.strip()

    # ─── GUARD: NO HITS ─── #
    if not row_order:
        print(f"[WARNING] No aligned sequences found in {aln_path}")
        return pa.table({
            "orf_id":      pa.array([], type=pa.int64()),
            "hmm_range":   pa.array([], type=pa.string()),
            "aligned_seq": pa.array([], type=pa.string()),
        })

    # ─── BUILD FULL TABLE (ALL ROWS) ─── #
    all_orf_ids:    list[int] = []
    all_hmm_ranges: list[str] = []
    all_sequences:  list[str] = []

    for seq_id in row_order:
        orf_id    = int(seq_id.split("|")[0])
        slash_idx = seq_id.rfind("/")
        hmm_range = seq_id[slash_idx + 1:] if slash_idx != -1 else ""
        all_orf_ids.append(orf_id)
        all_hmm_ranges.append(hmm_range)
        all_sequences.append(seq_chunks[seq_id])

    table = pa.table({
        "orf_id":      pa.array(all_orf_ids,    type=pa.int64()),
        "hmm_range":   pa.array(all_hmm_ranges, type=pa.string()),
        "aligned_seq": pa.array(all_sequences,  type=pa.string()),
    })

    # ─── WRITE .AFA FILE (ALL SEQUENCES, INCLUDING INSERTIONS) ─── #
    safe_makedirs(aligned_afa_path)
    with open(aligned_afa_path, "w") as f:
        for orf_id, hmm_range, seq in zip(all_orf_ids, all_hmm_ranges, all_sequences):
            f.write(f">{orf_id}/{hmm_range}\n{seq}\n")

    return table


def contains_catalytic_residues(annotation_path: str, component: str, hmmsearch_output_table: pa.Table, catalytic_csv_path: str) -> pa.Table:
    """
    Reads a CSV containing catalytic residue annotations and checks whether
    each sequence contains the correct amino acid at each catalytic position
    in the consensus sequence (uppercase match-state residues only).

    Annotated CSV format:
        component,aa_1,column_1,aa_2,column_2,aa_3,column_3
        c1,S,70,D,20,,
        c2,A,30,,,,
        c3,A,10,A,20,A,30

    The aa field supports slash-separated alternatives, e.g. "S/C" means
    either serine or cysteine is acceptable at that position.

    Annotations column format:
        "S/C127:TRUE D173:FALSE"
    """
    # ─── READ ANNOTATION CSV FOR THIS COMPONENT ─── #
    # Each catalytic site is (aa_spec, col) where aa_spec may be "S/C" etc.
    catalytic_sites: list[tuple[str, int]] = []

    with open(annotation_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["component"] != component:
                continue
            for i in range(1, 4):
                aa  = row.get(f"aa_{i}",     "").strip()
                col = row.get(f"column_{i}", "").strip()
                if aa and col:
                    catalytic_sites.append((aa, int(col)))
            break

    # ─── EXTRACT CONSENSUS (UPPERCASE ONLY, REMOVE LOWERCASE AND DOTS) ─── #
    def extract_consensus(aligned_seq: str) -> str:
        return re.sub(r'[a-z.]', '', aligned_seq)

    def residue_matches(residue: str, aa_spec: str) -> bool:
        """
        Returns True if residue matches the aa_spec.
        aa_spec may be a single letter ("S") or slash-separated alternatives
        ("S/C"), in which case any listed amino acid is accepted.
        """
        accepted = {a.strip().upper() for a in aa_spec.split("/")}
        return residue.upper() in accepted

    # ─── BUILD ANNOTATION STRING PER ORF ─── #
    orf_ids      = hmmsearch_output_table.column("orf_id").to_pylist()
    hmm_ranges   = hmmsearch_output_table.column("hmm_range").to_pylist()
    aligned_seqs = hmmsearch_output_table.column("aligned_seq").to_pylist()
    i_evalues    = hmmsearch_output_table.column("i_evalue").to_pylist()

    annotations = []
    for aligned_seq in aligned_seqs:
        consensus_seq = extract_consensus(aligned_seq)
        parts = []
        for aa_spec, col_1based in catalytic_sites:
            col_idx = col_1based - 1
            if col_idx < len(consensus_seq):
                is_present = residue_matches(consensus_seq[col_idx], aa_spec)
            else:
                is_present = False
            parts.append(f"{aa_spec}{col_1based}:{'TRUE' if is_present else 'FALSE'}")
        annotations.append(" ".join(parts))

    # ─── BUILD ANNOTATED TABLE ─── #
    annotated_table = pa.table({
        "orf_id":      pa.array(orf_ids,      type=pa.int64()),
        "hmm_range":   pa.array(hmm_ranges,   type=pa.string()),
        "i_evalue":    pa.array(i_evalues,    type=pa.float64()),
        "annotations": pa.array(annotations,  type=pa.string()),
    })

    # ─── WRITE CATALYTIC CSV ─── #
    safe_makedirs(catalytic_csv_path)
    try:
        pa_csv.write_csv(annotated_table, catalytic_csv_path)
    except Exception as e:
        raise RuntimeError(f"Failed to write catalytic CSV to {catalytic_csv_path}: {e}") from e

    return annotated_table


def annotate_stockholm_match_states(aligned_sequences: list[str]) -> str:
    """
    Creates a Stockholm #=GC RF annotation string for hmmbuild --hand.

    Input alignment convention:
        - Insertions are lowercase amino acids or '.'
        - Match states are uppercase amino acids or '-'

    Output RF convention:
        - 'x' = match column
        - '.' = insertion column

    A column is initially treated as a match state if at least one sequence has
    an uppercase amino acid or '-' at that position. If the match-state column
    has coverage below COVERAGE_CUTOFF, it is changed to an insertion state.

    Coverage is calculated as:
        number of uppercase amino acids / number of sequences

    Gaps ('-') mark the column as a match-state column, but they do not count
    toward coverage.
    """
    # ─── GUARD: NO SEQUENCES ─── #
    if not aligned_sequences:
        return ""

    # ─── GUARD: ALIGNMENT LENGTHS MUST MATCH ─── #
    aln_lengths = {len(seq) for seq in aligned_sequences}
    if len(aln_lengths) != 1:
        raise ValueError(f"Aligned sequences have unequal lengths: {sorted(aln_lengths)}")

    n_sequences = len(aligned_sequences)
    aln_length  = len(aligned_sequences[0])
    rf_chars: list[str] = []

    for col_idx in range(aln_length):
        column = [seq[col_idx] for seq in aligned_sequences]

        # ─── DETERMINE ORIGINAL ALIGNMENT STATE ─── #
        # Match columns are represented by uppercase residues or '-'.
        # Insert columns are represented by lowercase residues or '.'.
        is_match_state = any((char.isupper() or char == "-") for char in column)

        if not is_match_state:
            rf_chars.append(".")
            continue

        # ─── CALCULATE MATCH-STATE COVERAGE ─── #
        # Only uppercase residues count as coverage. Gaps mark a match-state
        # column but do not contribute to coverage.
        covered_count = sum(1 for char in column if char.isupper())
        coverage      = covered_count / n_sequences

        if coverage < COVERAGE_CUTOFF:
            rf_chars.append(".")
        else:
            rf_chars.append("x")

    return "".join(rf_chars)

def write_filtered_stockholm(
    annotated_table: pa.Table,
    hmmsearch_output_table: pa.Table,
    filtered_sto_path: str,
) -> tuple[int, list[tuple[str, str]]]:
    """
    Writes a filtered Stockholm alignment to results_path/filtered_sto/{component}.sto.

    A sequence is included only when BOTH conditions are met:
        1. i_evalue <= E_VALUE_CUTOFF
        2. Every catalytic residue annotation is TRUE

    Sequences whose orf_id appears more than once across the passing set are
    excluded from the primary .sto. Instead, their single best (lowest e-value)
    hit is written to repeat_sto_path as repeat_cata_{component}.sto.
    All sequences in both output files are catalytically intact.

    Parameters
    ----------
    annotated_table         : output of contains_catalytic_residues
                              columns: orf_id, hmm_range, i_evalue, annotations
    hmmsearch_output_table  : joined table with aligned_seq
                              columns: orf_id, hmm_range, i_evalue, aligned_seq
    filtered_sto_path       : full path to the primary output .sto file
    repeat_sto_path         : full path to the repeat-hit output .sto file

    Returns
    -------
    tuple[int, int, list[tuple[str, str]]] :
        (sequences written to primary sto,
         sequences written to repeat sto,
         filtered_rows as list of (seq_id, aligned_seq))
    """
    # ─── BUILD ALIGNED_SEQ LOOKUP: (orf_id, hmm_range) -> seq queue ─── #
    seq_queue: dict[tuple[int, str], list[str]] = defaultdict(list)
    for oid, rng, seq in zip(
        hmmsearch_output_table.column("orf_id").to_pylist(),
        hmmsearch_output_table.column("hmm_range").to_pylist(),
        hmmsearch_output_table.column("aligned_seq").to_pylist(),
    ):
        seq_queue[(oid, rng)].append(seq)

    # ─── APPLY FILTERS ─── #
    def all_catalytic_true(annotation_str: str) -> bool:
        if not annotation_str:
            return False
        return all(token.endswith(":TRUE") for token in annotation_str.split())

    # Include i_evalue in passing so write_repeat_catalytic_stockholm can rank hits
    passing: list[tuple[int, str, str, float]] = []  # (orf_id, hmm_range, aligned_seq, i_evalue)
    passing_orf_counts: dict[int, int] = {}

    for oid, rng, ev, ann in zip(
        annotated_table.column("orf_id").to_pylist(),
        annotated_table.column("hmm_range").to_pylist(),
        annotated_table.column("i_evalue").to_pylist(),
        annotated_table.column("annotations").to_pylist(),
    ):
        if ev > E_VALUE_CUTOFF:
            continue
        if not all_catalytic_true(ann):
            continue

        key = (oid, rng)
        seq = seq_queue[key].pop(0) if seq_queue[key] else ""
        if not seq:
            continue

        passing.append((oid, rng, seq, ev))
        passing_orf_counts[oid] = passing_orf_counts.get(oid, 0) + 1

    # ─── SPLIT INTO UNIQUE AND DUPLICATE ORF_IDS ─── #
    duplicate_orf_ids: set[int] = {
        oid for oid, count in passing_orf_counts.items() if count > 1
    }
    if duplicate_orf_ids:
        print(f"[filtered_sto] {len(duplicate_orf_ids)} orf_id(s) with >1 catalytic hit — "
              f"best hit per orf written to repeat sto: {sorted(duplicate_orf_ids)}")

    # ─── BUILD BEST HIT PER DUPLICATE FOR PRIMARY STO ─── #
    best_duplicate_hit: dict[int, tuple[str, str, float]] = {}  # orf_id -> (hmm_range, seq, evalue)
    for oid, rng, seq, ev in passing:
        if oid not in duplicate_orf_ids:
            continue
        if oid not in best_duplicate_hit or ev < best_duplicate_hit[oid][2]:
            best_duplicate_hit[oid] = (rng, seq, ev)

    # Primary sto: unique orf_ids + best hit for each duplicate
    filtered_rows: list[tuple[str, str]] = [
        (f"{oid}/{rng}", seq)
        for oid, rng, seq, ev in passing
        if oid not in duplicate_orf_ids
    ] + [
        (f"{oid}/{rng}", seq)
        for oid, (rng, seq, _) in sorted(best_duplicate_hit.items())
    ]

    # ─── WRITE PRIMARY STOCKHOLM ─── #
    safe_makedirs(filtered_sto_path)
    with open(filtered_sto_path, "w") as f:
        f.write("# STOCKHOLM 1.0\n")

        if filtered_rows:
            seq_ids = [seq_id for seq_id, _ in filtered_rows]
            seqs    = [seq    for _, seq    in filtered_rows]
            rf      = annotate_stockholm_match_states(seqs)

            max_name_len = max(len(seq_id) for seq_id in seq_ids + ["#=GC RF"])
            for seq_id, seq in filtered_rows:
                f.write(f"{seq_id:<{max_name_len}}  {seq}\n")
            f.write(f"{'#=GC RF':<{max_name_len}}  {rf}\n")

        f.write("//\n")

    return len(filtered_rows), filtered_rows

def _count_trues(annotation_str: str) -> int:
    """Returns the number of ':TRUE' tokens in an annotation string."""
    if not annotation_str:
        return 0
    return sum(1 for token in annotation_str.split() if token.endswith(":TRUE"))


def _annotation_sort_key(annotation_str: str) -> int:
    """Sort key: descending TRUE count (negate for ascending sort calls)."""
    return -_count_trues(annotation_str)


def print_summary(
    component:              str,
    e_value_table:          pa.Table,
    aligned_table:          pa.Table,
    annotated_table:        pa.Table,
    hmmsearch_output_table: pa.Table,
    summary_dir:            str,
):
    """
    Prints a textual and graphical summary of the HMMsearch results.

    Outputs
    -------
    Printed to stdout:
        - Lowest and highest 10 i-evalue hits (from the full, unfiltered domtbl table)
        - Number of unique orf_ids present in the aligned Stockholm output
        - Top 3 orf_ids by number of HMM domain hits

    Saved to summary_dir/{component}_annotations.png:
        - Bar plot of counts for each unique annotation string, sorted left
          (most TRUEs) to right (least TRUEs). Includes total sequences
          searched and sequences represented by the plot in the title.

    Saved to summary_dir/{component}_hits_per_orf.png:
        - Histogram of number of HMM domain hits per orf_id (from the full domtbl table)

    Saved to summary_dir/{component}_multi_hit_matrix.png:
        - Matrix plot for orf_ids with >1 HMM hit showing the annotation of
          their lowest e-value hit (x-axis) vs second-lowest e-value hit
          (y-axis). Both axes are sorted most TRUEs (left/top) to least
          TRUEs (right/bottom). Cell values are counts of orf_ids.

    Parameters
    ----------
    component              : component name, used for titles and filenames
    e_value_table          : full domtbl table (all hits, no e-value filter)
                             columns: orf_id, hmm_range, i_evalue
    aligned_table          : Stockholm-parsed table (all aligned sequences)
                             columns: orf_id, hmm_range, aligned_seq
    annotated_table        : catalytic annotation table
                             columns: orf_id, hmm_range, i_evalue, annotations
    hmmsearch_output_table : joined table (e-value rows + aligned_seq)
                             columns: orf_id, hmm_range, i_evalue, aligned_seq
    summary_dir            : directory to write plot images into
    """

    os.makedirs(summary_dir, exist_ok=True)

    SEP = "─" * 60

    # ─── E-VALUE EXTREMES (pre-filtered, all domtbl hits) ─── #
    all_orf_ids    = e_value_table.column("orf_id").to_pylist()
    all_hmm_ranges = e_value_table.column("hmm_range").to_pylist()
    all_evalues    = e_value_table.column("i_evalue").to_pylist()

    rows = sorted(zip(all_evalues, all_orf_ids, all_hmm_ranges))

    # ─── BUILD ALIGNED_SEQ LOOKUP: (orf_id, hmm_range) -> seq ─── #
    # Used to print the aligned sequence alongside each orf_id where available.
    # Queue-based so multiple hits at the same key are each shown once.
    seq_lookup: dict[tuple[int, str], list[str]] = defaultdict(list)
    for oid, rng, seq in zip(
        hmmsearch_output_table.column("orf_id").to_pylist(),
        hmmsearch_output_table.column("hmm_range").to_pylist(),
        hmmsearch_output_table.column("aligned_seq").to_pylist(),
    ):
        if seq:
            seq_lookup[(oid, rng)].append(seq)

    def get_seq(oid: int, rng: str) -> str:
        """Peek at the first available aligned sequence for this (orf_id, hmm_range), or ''."""
        seqs = seq_lookup.get((oid, rng), [])
        return seqs[0] if seqs else ""

    def print_hit(ev: float, oid: int, rng: str) -> None:
        seq = get_seq(oid, rng)
        print(f"  orf_id={oid}  hmm_range={rng}  i_evalue={ev:.3e}")
        if seq:
            print(f"    seq: {seq}")
        else:
            print(f"    seq: (no alignment)")

    print(SEP)
    print(f"[{component}] SUMMARY")
    print(SEP)

    print(f"\n{'─'*30}")
    print(f"  LOWEST 10 E-VALUE HITS  (best hits)")
    print(f"{'─'*30}")
    for ev, oid, rng in rows[:10]:
        print_hit(ev, oid, rng)

    print(f"\n{'─'*30}")
    print(f"  HIGHEST 10 E-VALUE HITS  (worst hits)")
    print(f"{'─'*30}")
    for ev, oid, rng in rows[-10:]:
        print_hit(ev, oid, rng)

    # ─── ALIGNED ORF COUNT ─── #
    unique_aln_orfs = len(set(aligned_table.column("orf_id").to_pylist()))
    print(f"\n{SEP}")
    print(f"  Unique orf_ids in HMMsearch aligned output: {unique_aln_orfs}")
    print(SEP)

    # ─── BAR PLOT: ANNOTATION VALUE COUNTS, SORTED BY TRUE COUNT ─── #
    # Sorted left (most TRUEs) to right (least TRUEs).
    # Title includes total sequences searched and sequences represented.
    #
    # For each orf_id, determine whether it has exactly 1 or >1 HMM domain hits
    # (from the full unfiltered domtbl e_value_table, not the filtered table).
    orf_hit_counts: Counter = Counter(all_orf_ids)  # all_orf_ids from e_value_table above

    ann_orf_ids     = annotated_table.column("orf_id").to_pylist()
    annotation_list = annotated_table.column("annotations").to_pylist()

    single_counts: Counter = Counter()
    multi_counts:  Counter = Counter()
    for oid, ann in zip(ann_orf_ids, annotation_list):
        if orf_hit_counts[oid] == 1:
            single_counts[ann] += 1
        else:
            multi_counts[ann] += 1

    # ─── SORT LABELS: MOST TRUEs LEFT, LEAST TRUEs RIGHT ─── #
    all_ann_labels = sorted(
        set(single_counts) | set(multi_counts),
        key=_annotation_sort_key,
    )

    single_vals = [single_counts[l] for l in all_ann_labels]
    multi_vals  = [multi_counts[l]  for l in all_ann_labels]

    # ─── COUNT TOTALS FOR TITLE ─── #
    # Sequences searched = all rows in the annotated table (one per domtbl hit row).
    # Sequences represented = unique orf_ids present in annotated table.
    n_searched_bar     = len(ann_orf_ids)
    n_represented_bar  = len(set(ann_orf_ids))

    bar_width = 0.4
    x_pos = range(len(all_ann_labels))

    fig, ax = plt.subplots(figsize=(max(8, len(all_ann_labels) * 1.4), 6))

    bars_single = ax.bar(
        [x - bar_width / 2 for x in x_pos], single_vals,
        width=bar_width, label="= 1 HMM hit",
        color="#4C72B0", edgecolor="white", linewidth=0.6,
    )
    bars_multi = ax.bar(
        [x + bar_width / 2 for x in x_pos], multi_vals,
        width=bar_width, label="> 1 HMM hit",
        color="#DD8452", edgecolor="white", linewidth=0.6,
    )

    # Annotate bar tops
    all_vals = single_vals + multi_vals
    top_val  = max(all_vals) if all_vals else 1
    for bar, val in list(zip(bars_single, single_vals)) + list(zip(bars_multi, multi_vals)):
        if val > 0:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + top_val * 0.01,
                str(val),
                ha="center", va="bottom", fontsize=7,
            )

    ax.set_xticks(list(x_pos))
    ax.set_xticklabels(all_ann_labels, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title(
        f"{component} — Catalytic annotation counts  (= 1 vs > 1 HMM hit)\n"
        f"Sequences searched: {n_searched_bar:,}  |  Unique orf_ids represented: {n_represented_bar:,}",
        fontsize=13,
    )
    ax.legend(fontsize=10)

    plt.tight_layout()
    ann_plot_path = os.path.join(summary_dir, f"{component}_annotations.png")
    fig.savefig(ann_plot_path, dpi=150)
    plt.close(fig)
    print(f"[{component}] Annotation bar plot saved to {ann_plot_path}")

    # ─── MATRIX PLOT: MULTI-HIT ORF_IDS ─── #
    # For every orf_id that has >1 HMM hit, identify its two lowest e-value
    # hits from the annotated table and cross-tabulate their annotation strings.
    #
    # X-axis: annotation of the lowest e-value hit
    # Y-axis: annotation of the second-lowest e-value hit
    # Both axes sorted most TRUEs (left/top) to least TRUEs (right/bottom).
    # Cell values: count of orf_ids matching that (x, y) annotation pair.

    # Collect all (orf_id, i_evalue, annotation) rows from annotated_table
    ann_evalues = annotated_table.column("i_evalue").to_pylist()

    # Group hits by orf_id, keeping track of evalue and annotation
    orf_hits: dict[int, list[tuple[float, str]]] = defaultdict(list)
    for oid, ev, ann in zip(ann_orf_ids, ann_evalues, annotation_list):
        orf_hits[oid].append((ev, ann))

    # Extract the two lowest e-value annotations for each multi-hit orf_id
    matrix_pairs: list[tuple[str, str]] = []  # (best_ann, second_ann)
    for oid, hits in orf_hits.items():
        if len(hits) < 2:
            continue
        sorted_hits = sorted(hits, key=lambda h: h[0])  # ascending evalue
        best_ann   = sorted_hits[0][1]
        second_ann = sorted_hits[1][1]
        matrix_pairs.append((best_ann, second_ann))

    if matrix_pairs:
        # Determine unique annotation labels present in the matrix data,
        # sorted most TRUEs (index 0) to least TRUEs (last index).
        x_labels_set = {p[0] for p in matrix_pairs}
        y_labels_set = {p[1] for p in matrix_pairs}
        all_matrix_labels = sorted(
            x_labels_set | y_labels_set,
            key=_annotation_sort_key,
        )

        x_labels = sorted(x_labels_set, key=_annotation_sort_key)
        y_labels = sorted(y_labels_set, key=_annotation_sort_key)

        x_idx = {lbl: i for i, lbl in enumerate(x_labels)}
        y_idx = {lbl: i for i, lbl in enumerate(y_labels)}

        matrix = np.zeros((len(y_labels), len(x_labels)), dtype=int)
        for best_ann, second_ann in matrix_pairs:
            matrix[y_idx[second_ann], x_idx[best_ann]] += 1

        # ─── COUNTS FOR TITLE ─── #
        n_multi_orf      = len(matrix_pairs)          # unique orf_ids in the plot
        n_searched_mat   = sum(                        # total domtbl rows for multi-hit orfs
            len(v) for k, v in orf_hits.items() if len(v) >= 2
        )

        fig, ax = plt.subplots(figsize=(max(6, len(x_labels) * 1.6), max(5, len(y_labels) * 1.2)))
        im = ax.imshow(matrix, aspect="auto", cmap="Blues")

        # Annotate each cell with its count
        for row_i in range(len(y_labels)):
            for col_i in range(len(x_labels)):
                val = matrix[row_i, col_i]
                ax.text(
                    col_i, row_i, str(val),
                    ha="center", va="center",
                    fontsize=9,
                    color="white" if val > matrix.max() * 0.6 else "black",
                )

        ax.set_xticks(range(len(x_labels)))
        ax.set_xticklabels(x_labels, rotation=40, ha="right", fontsize=8)
        ax.set_yticks(range(len(y_labels)))
        ax.set_yticklabels(y_labels, fontsize=8)
        ax.set_xlabel("Annotation — lowest e-value HMM hit", fontsize=11)
        ax.set_ylabel("Annotation — 2nd lowest e-value HMM hit", fontsize=11)
        ax.set_title(
            f"{component} — Multi-hit orf_id annotation matrix\n"
            f"Sequences searched: {n_searched_mat:,}  |  orf_ids represented: {n_multi_orf:,}",
            fontsize=12,
        )
        plt.colorbar(im, ax=ax, label="orf_id count")
        plt.tight_layout()
        matrix_plot_path = os.path.join(summary_dir, f"{component}_multi_hit_matrix.png")
        fig.savefig(matrix_plot_path, dpi=150)
        plt.close(fig)
        print(f"[{component}] Multi-hit matrix plot saved to {matrix_plot_path}")
    else:
        print(f"[{component}] No orf_ids with >1 HMM hit found; matrix plot skipped.")


def main():
    # ─── DIRECTORY FORMAT ─── #
    #
    # results_path/
    #   aln/            {component}.aln
    #   domtbl/         {component}.domtbl
    #   e_value_csv/    {component}.csv       [CREATED]
    #   aligned_afa/    {component}.afa       [CREATED]
    #   filtered_afa/   {component}.afa       [CREATED]
    #   filtered_sto/   {component}.sto       [CREATED]
    #   fastaa/         {component}.fa        [CREATED]
    #   summary/        {component}_annotations.png      [CREATED]
    #   summary/        {component}_hits_per_orf.png     [CREATED]
    #   summary/        {component}_multi_hit_matrix.png [CREATED]
    # analysis_path/
    #   catalytic/      {component}.csv       [CREATED]
    # fragments_dir/
    #   {component}.afa                       [CREATED]

    # ─── VALIDATE ARGUMENTS ─── #
    if len(sys.argv) != 6:
        print("Usage: python 260406_clean_hmmsearch.py <results_path> <component> <analysis_path> <annotated_path> <fragments_dir>")
        print("See module docstring for full argument descriptions.")
        sys.exit(1)

    results_path   = sys.argv[1]
    component      = sys.argv[2]
    analysis_path  = sys.argv[3]
    annotated_path = sys.argv[4]
    fragments_dir  = sys.argv[5]

    aln_path              = os.path.join(results_path,  "aln",         component + ".aln")
    domtbl_path           = os.path.join(results_path,  "domtbl",      component + ".domtbl")
    e_value_csv_path      = os.path.join(results_path,  "e_value_csv", component + ".csv")
    aligned_afa_path      = os.path.join(results_path,  "aligned_afa", component + ".afa")
    catalytic_csv         = os.path.join(analysis_path, "catalytic",   component + ".csv")

    # ─── PARSE HMMSEARCH OUTPUTS ─── #
    print(f"[{component}] Parsing domtblout...")
    e_value_table = parse_domtblout(domtbl_path, e_value_csv_path)
    print(f"[{component}] {len(e_value_table)} hits below E-value cutoff {E_VALUE_CUTOFF}")

    print(f"[{component}] Extracting aligned sequences...")
    aligned_table = extracted_aligned_sequences(aln_path, aligned_afa_path)
    print(f"[{component}] {len(aligned_table)} aligned sequences extracted")

    # ─── JOIN ALIGNED SEQUENCES ONTO E-VALUE ROWS FOR CATALYTIC CHECKING ─── #
    # The catalytic CSV is built from ALL e-value rows (every individual domtbl hit,
    # including duplicate orf_ids). aligned_seq is looked up as a best-effort for
    # the residue check; e-value rows with no matching alignment get an empty string.
    print(f"[{component}] Joining aligned sequences onto e-value rows...")

    # (orf_id, hmm_range) -> queue of aligned_seqs (consumed one per matching e-value row)
    aligned_seq_queue: dict[tuple[int, str], list[str]] = defaultdict(list)
    for oid, rng, seq in zip(
        aligned_table.column("orf_id").to_pylist(),
        aligned_table.column("hmm_range").to_pylist(),
        aligned_table.column("aligned_seq").to_pylist(),
    ):
        aligned_seq_queue[(oid, rng)].append(seq)

    # Iterate every e-value row; pop one aligned_seq per matching key if available.
    all_orf_ids:    list[int]   = e_value_table.column("orf_id").to_pylist()
    all_hmm_ranges: list[str]   = e_value_table.column("hmm_range").to_pylist()
    all_i_evalues:  list[float] = e_value_table.column("i_evalue").to_pylist()
    all_seqs:       list[str]   = []

    for oid, rng in zip(all_orf_ids, all_hmm_ranges):
        key = (oid, rng)
        if aligned_seq_queue[key]:
            all_seqs.append(aligned_seq_queue[key].pop(0))
        else:
            all_seqs.append("")

    hmmsearch_output_table = pa.table({
        "orf_id":      pa.array(all_orf_ids,    type=pa.int64()),
        "hmm_range":   pa.array(all_hmm_ranges, type=pa.string()),
        "i_evalue":    pa.array(all_i_evalues,  type=pa.float64()),
        "aligned_seq": pa.array(all_seqs,       type=pa.string()),
    })
    print(f"[{component}] {len(hmmsearch_output_table)} sequences in output table")

    # ─── CHECK FOR CATALYTIC RESIDUES ─── #
    print(f"[{component}] Checking catalytic residues...")
    annotated_table = contains_catalytic_residues(annotated_path, component, hmmsearch_output_table, catalytic_csv)
    print(f"[{component}] Catalytic annotations written to {catalytic_csv}")

    # ─── WRITE FILTERED STOCKHOLM ─── #
    filtered_sto_path = os.path.join(results_path, "filtered_sto", component + ".sto")
    print(f"[{component}] Writing filtered Stockholm alignment for hmmbuild --hand...")
    n_written, n_repeat = write_filtered_stockholm(
        annotated_table, hmmsearch_output_table, filtered_sto_path
    )
    print(f"[{component}] {n_written} sequences written to {filtered_sto_path}")
    # ─── SUMMARY ─── #
    summary_dir = os.path.join(results_path, "summary")
    print_summary(
        component              = component,
        e_value_table          = e_value_table,
        aligned_table          = aligned_table,
        annotated_table        = annotated_table,
        hmmsearch_output_table = hmmsearch_output_table,
        summary_dir            = summary_dir,
    )

    print(f"[{component}] Done.")


if __name__ == "__main__":
    main()
