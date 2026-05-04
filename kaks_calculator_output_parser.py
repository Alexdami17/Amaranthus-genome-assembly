#!/usr/bin/env python3
"""
Parse KaKs Calculator output and optionally integrate with a GFF3 file.

Modes
-----
1) Without GFF (simple KaKs table):

    python kaks_gff_sort.py \
        --kaks kaks_output.txt \
        --out kaks_parsed.tsv

    Output columns:
        gene1_id  gene2_id  Ka  Ks  Ka_Ks  Pvalue

2) With GFF (sorted by gene start coordinate):

    python kaks_gff_sort.py \
        --kaks kaks_output.txt \
        --gff genes.gff3 \
        --out kaks_sorted_by_start.tsv

    Output columns:
        gene_id  start  Ka  Ks  Ka_Ks  Pvalue

3) With GFF + extra KaKs table:

    python kaks_gff_sort.py \
        --kaks kaks_output.txt \
        --gff genes.gff3 \
        --out kaks_sorted_by_start.tsv \
        --extra-kaks kaks_parsed.tsv
"""

import argparse
import csv
import logging
import re
from typing import Dict, List, Tuple, Optional, Set


LOG = logging.getLogger("kaks_gff")


# ----------------------------------------------------------------------
# GFF parsing
# ----------------------------------------------------------------------
def parse_gff(gff_path: str):
    """
    Parse a GFF3 file.

    Returns
    -------
    genes : dict
        gene_id -> (seqid, start)
    transcript_to_gene : dict
        transcript_id -> gene_id
    """
    genes: Dict[str, Tuple[str, int]] = {}
    transcript_to_gene: Dict[str, str] = {}

    with open(gff_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs = parts

            # Parse attributes into a dict
            attr_dict = {}
            for attr in attrs.split(";"):
                if not attr:
                    continue
                if "=" in attr:
                    k, v = attr.split("=", 1)
                    attr_dict[k] = v

            if ftype == "gene":
                gene_id = attr_dict.get("ID")
                if gene_id is None:
                    continue
                try:
                    start_pos = int(start)
                except ValueError:
                    continue
                genes[gene_id] = (seqid, start_pos)

            elif ftype in ("mRNA", "transcript"):
                tid = attr_dict.get("ID")
                parent = attr_dict.get("Parent")
                if tid and parent:
                    # Some formats have comma-separated parents; take first
                    parent_gene = parent.split(",")[0]
                    transcript_to_gene[tid] = parent_gene

    LOG.info(
        "Parsed GFF: %d genes, %d transcripts",
        len(genes),
        len(transcript_to_gene),
    )
    return genes, transcript_to_gene


# ----------------------------------------------------------------------
# KaKs parsing
# ----------------------------------------------------------------------
def split_sequence_field(
    seq_field: str,
    known_ids: Optional[Set[str]] = None,
) -> Tuple[str, str]:
    """
    Split the KaKs 'Sequence' field into gene1_id and gene2_id.

    If known_ids is provided (from GFF: transcript + gene IDs), we try to find
    a split such that both parts are in known_ids. That avoids relying on
    naming conventions.

    Fallbacks:
      - regex suited to "-mRNA-<n>" pattern
      - split on first '-'
    """
    seq_field = seq_field.strip()

    # 1) Smart split using known IDs from GFF
    if known_ids:
        for i in range(1, len(seq_field)):
            left = seq_field[:i]
            right = seq_field[i:]
            if left in known_ids and right in known_ids:
                return left, right
        # If we reach here, we didn't find a clean split using known IDs
        LOG.debug(
            "Could not split '%s' using known IDs; falling back to heuristics",
            seq_field,
        )

    # 2) Heuristic: typical "-mRNA-<n>" style
    m = re.match(r"(.+?-mRNA-\d+)-(.*)", seq_field)
    if m:
        return m.group(1), m.group(2)

    # 3) Fallback: split at first '-'
    parts = seq_field.split("-", 1)
    if len(parts) == 2:
        return parts[0], parts[1]

    # 4) Last resort: treat everything as gene1
    LOG.warning("Could not confidently split Sequence field '%s'", seq_field)
    return seq_field, ""


def parse_kaks_file(
    kaks_path: str,
    known_ids_for_split: Optional[Set[str]] = None,
) -> List[dict]:
    """
    Parse KaKs Calculator output.

    Removes any row where Ka, Ks, or Ka_Ks is NA (case-insensitive).
    Returns
    -------
    rows_for_output : list of dicts with keys:
        gene1_id, gene2_id, Ka, Ks, Ka_Ks, Pvalue
    """
    rows_for_output: List[dict] = []
    total_lines = 0
    removed_na = 0

    with open(kaks_path, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue
            total_lines += 1

            if row[0].startswith("#"):
                continue
            # Skip header line
            if row[0].strip().lower().startswith("sequence"):
                continue
            # Require at least 6 columns: Sequence, Method, Ka, Ks, Ka/Ks, P-Value(Fisher)
            if len(row) < 6:
                LOG.debug("Skipping short KaKs row (len=%d): %s", len(row), row)
                continue

            seq_field = row[0].strip()
            ka = row[2].strip()
            ks = row[3].strip()
            ka_ks = row[4].strip()
            pval = row[5].strip().rstrip(">")  # strip trailing '>' if present

            # NA filter (exact 'NA', case-insensitive)
            if ka.upper() == "NA" or ks.upper() == "NA" or ka_ks.upper() == "NA":
                removed_na += 1
                continue

            gene1_id, gene2_id = split_sequence_field(seq_field, known_ids_for_split)

            row_dict = {
                "gene1_id": gene1_id,
                "gene2_id": gene2_id,
                "Ka": ka,
                "Ks": ks,
                "Ka_Ks": ka_ks,
                "Pvalue": pval,
            }
            rows_for_output.append(row_dict)

    LOG.info(
        "Parsed KaKs file: %d data rows (from %d lines, %d removed due to NA)",
        len(rows_for_output),
        total_lines,
        removed_na,
    )
    return rows_for_output


def write_kaks_tsv(rows: List[dict], out_path: str):
    """Write the parsed KaKs data to a TSV file."""
    fieldnames = ["gene1_id", "gene2_id", "Ka", "Ks", "Ka_Ks", "Pvalue"]
    with open(out_path, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    LOG.info("Wrote KaKs table to %s (%d rows)", out_path, len(rows))


# ----------------------------------------------------------------------
# Mapping KaKs -> genes via GFF
# ----------------------------------------------------------------------
def build_gene_to_kaks(
    kaks_rows: List[dict],
    genes: Dict[str, Tuple[str, int]],
    transcript_to_gene: Dict[str, str],
) -> Dict[str, List[dict]]:
    """
    Build mapping: gene_id -> list of KaKs rows, based on gene1_id.

    Strategy:
      1) If gene1_id is a transcript in transcript_to_gene, map to its gene.
      2) Else if gene1_id itself is a gene ID in 'genes', use it directly.
      3) Else try some simple suffix-strip heuristics (e.g. '-mRNA-1', '-T1').
      4) If still not found, skip with a warning.

    Returns
    -------
    gene_to_kaks : dict
        gene_id -> list of KaKs row dicts
    """
    gene_to_kaks: Dict[str, List[dict]] = {}
    skipped = 0
    mapped = 0

    for rec in kaks_rows:
        tid = rec["gene1_id"]
        gene_id: Optional[str] = None

        # 1) Transcript -> gene via GFF mapping
        if tid in transcript_to_gene:
            gene_id = transcript_to_gene[tid]

        # 2) Exact gene ID
        elif tid in genes:
            gene_id = tid

        else:
            # 3) Heuristic: strip common transcript suffixes
            base = re.sub(r"(-mRNA-\d+|-T\d+|-t\d+)$", "", tid)
            if base in genes:
                gene_id = base

        if gene_id is None:
            skipped += 1
            LOG.warning(
                "Could not map KaKs gene1_id '%s' to any gene in GFF; skipping this pair for GFF-based output",
                tid,
            )
            continue

        if gene_id not in genes:
            skipped += 1
            LOG.warning(
                "Mapped gene_id '%s' from '%s' not found as a gene feature in GFF; skipping",
                gene_id,
                tid,
            )
            continue

        mapped += 1
        gene_to_kaks.setdefault(gene_id, []).append(rec)

    LOG.info(
        "Mapped KaKs -> GFF genes: %d mapped, %d skipped (total KaKs rows: %d)",
        mapped,
        skipped,
        len(kaks_rows),
    )
    return gene_to_kaks


def get_sorted_gene_list(genes: Dict[str, Tuple[str, int]]) -> List[Tuple[str, str, int]]:
    """
    Convert genes dict to a sorted list.

    Returns
    -------
    list of (gene_id, seqid, start), sorted by seqid then start.
    """
    out: List[Tuple[str, str, int]] = []
    for gid, (seqid, start) in genes.items():
        out.append((gid, seqid, start))
    out.sort(key=lambda x: (x[1], x[2]))
    return out


def write_sorted_by_gff(
    genes_sorted: List[Tuple[str, str, int]],
    gene_to_kaks: Dict[str, List[dict]],
    out_path: str,
):
    """
    Create a TSV where rows are KaKs entries whose gene1 base ID
    maps to a 'gene' feature in the GFF, sorted by gene start.

    Output columns:
        gene_id  start  Ka  Ks  Ka_Ks  Pvalue
    """
    written = 0
    with open(out_path, "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(["gene_id", "start", "Ka", "Ks", "Ka_Ks", "Pvalue"])

        for gene_id, seqid, start in genes_sorted:
            if gene_id not in gene_to_kaks:
                continue
            for rec in gene_to_kaks[gene_id]:
                writer.writerow([
                    gene_id,
                    start,
                    rec["Ka"],
                    rec["Ks"],
                    rec["Ka_Ks"],
                    rec["Pvalue"],
                ])
                written += 1

    LOG.info(
        "Wrote GFF-sorted KaKs table to %s (%d rows, %d genes with KaKs)",
        out_path,
        written,
        len(gene_to_kaks),
    )


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description=(
            "Parse KaKs Calculator output and optionally integrate with a GFF3 file.\n\n"
            "Without --gff: --out is a KaKs table.\n"
            "With --gff:    --out is a GFF-sorted table by gene start.\n"
            "Optional --extra-kaks always writes the raw KaKs table."
        )
    )
    parser.add_argument(
        "-k", "--kaks",
        required=True,
        help="KaKs Calculator output file (tab-delimited).",
    )
    parser.add_argument(
        "-g", "--gff",
        required=False,
        help="GFF3 file containing gene/mRNA annotations (optional).",
    )
    parser.add_argument(
        "-o", "--out",
        required=True,
        help="Main output TSV file. "
             "Without --gff: KaKs table. With --gff: GFF-sorted table.",
    )
    parser.add_argument(
        "--extra-kaks",
        required=False,
        help="Optional additional TSV file for the KaKs table "
             "(gene1_id, gene2_id, Ka, Ks, Ka_Ks, Pvalue), "
             "written in all modes.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO).",
    )

    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    # If GFF is supplied, we parse it first to get ID sets for smarter splitting
    if args.gff:
        genes, transcript_to_gene = parse_gff(args.gff)
        known_ids_for_split = set(genes.keys()) | set(transcript_to_gene.keys())
        kaks_rows = parse_kaks_file(args.kaks, known_ids_for_split=known_ids_for_split)

        # Extra KaKs table if requested
        if args.extra_kaks:
            write_kaks_tsv(kaks_rows, args.extra_kaks)

        # Build gene -> KaKs mapping and write sorted output
        gene_to_kaks = build_gene_to_kaks(kaks_rows, genes, transcript_to_gene)
        genes_sorted = get_sorted_gene_list(genes)
        write_sorted_by_gff(genes_sorted, gene_to_kaks, args.out)

    else:
        # No GFF: simple KaKs parsing and write to --out
        kaks_rows = parse_kaks_file(args.kaks, known_ids_for_split=None)
        write_kaks_tsv(kaks_rows, args.out)
        # Also write to extra_kaks if requested (same content)
        if args.extra_kaks and args.extra_kaks != args.out:
            write_kaks_tsv(kaks_rows, args.extra_kaks)


if __name__ == "__main__":
    main()
