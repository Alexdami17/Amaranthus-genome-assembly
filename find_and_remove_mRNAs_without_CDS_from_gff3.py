#!/usr/bin/env python3

"""
find_and_remove_mrna_without_cds.py

Purpose
-------
Single tool to:
  (1) scan a GFF3 (or .gz) to find mRNA/transcript features lacking any CDS;
  (2) log a summary (counts + IDs);
  (3) if any are missing CDS, remove those mRNAs/transcripts AND any genes that
      become orphaned (i.e., all their mRNAs removed), plus any child features
      whose Parent includes a removed mRNA or dropped gene.

Behavior
--------
- Treats both "mRNA" and "transcript" as mRNA-like.
- Uses tqdm progress bars if installed; otherwise runs quietly.
- If no mRNA-like features are missing CDS: writes the log and exits with 0,
  and DOES NOT write an output unless --always-write is provided.
- If removals are needed but --output is not provided, logs the summary and exits with code 2.

Usage
-----
  python find_and_remove_mrna_without_cds.py \
      -i input.gff3[.gz] \
      -o filtered.gff3 \
      -l input.missingCDS.log

  # Disable progress bars even if tqdm is installed
  python find_and_remove_mrna_without_cds.py -i input.gff3 -o filtered.gff3 --no-progress

  # If nothing needs removal, still write a verbatim copy
  python find_and_remove_mrna_without_cds.py -i input.gff3 -o copy.gff3 --always-write
"""

from __future__ import annotations
from collections import defaultdict
import argparse
import gzip
import os
import re
import sys
from typing import Dict, Iterable, List, Set, Tuple

# ────────────────────────────────────────────────────────────────────────────────
# Optional tqdm
# ────────────────────────────────────────────────────────────────────────────────
try:
    from tqdm import tqdm as _tqdm  # type: ignore

    def tqdm_wrap(it: Iterable[str], **kw):
        return _tqdm(it, **kw)

    TQDM = True
except Exception:
    TQDM = False

    def tqdm_wrap(it: Iterable[str], **kw):
        return it


# ────────────────────────────────────────────────────────────────────────────────
# I/O helpers
# ────────────────────────────────────────────────────────────────────────────────
def open_any(path: str, mode: str = "rt"):
    """Open plain or gzipped files with the same call."""
    if path.endswith(".gz"):
        return gzip.open(path, mode)  # "rt"/"wt" supported
    return open(path, mode)


def parse_attrs(attr_str: str) -> Dict[str, str]:
    """Parse a GFF3 attributes column into a dict (minimal, robust)."""
    out: Dict[str, str] = {}
    if not attr_str or attr_str == ".":
        return out
    for part in attr_str.split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        else:
            # flag-like attributes are rare; keep as true-ish string
            out[part] = "true"
    return out


_ATTR_RE_CACHE: Dict[str, re.Pattern] = {}


def get_attr(attrs: str, key: str) -> str | None:
    """Fast-ish regex getter for a single attribute value."""
    pat = _ATTR_RE_CACHE.get(key)
    if pat is None:
        pat = re.compile(rf"(?:^|;){re.escape(key)}=([^;]+)")
        _ATTR_RE_CACHE[key] = pat
    m = pat.search(attrs)
    return m.group(1) if m else None


# ────────────────────────────────────────────────────────────────────────────────
# Pass 1: collect transcript + CDS info
# ────────────────────────────────────────────────────────────────────────────────
def pass1_collect(
    input_path: str,
    show_progress: bool,
    max_ids_in_log: int | None = None,
) -> Tuple[Set[str], Set[str], Dict[str, Set[str]], List[str]]:
    """
    Returns:
      - mrnas: set of all mRNA-like IDs
      - missing_mrnas: set of mRNA-like IDs with ZERO CDS
      - gene_to_mrnas: map of gene ID -> set of its mRNA-like child IDs
      - log_lines: ready-to-write summary lines
    """
    mrnas: Set[str] = set()
    mrnas_with_cds: Set[str] = set()
    gene_to_mrnas: Dict[str, Set[str]] = defaultdict(set)
    log_lines: List[str] = []

    f = open_any(input_path, "rt")
    try:
        iterator = tqdm_wrap(f, desc="Pass1 scan", unit="line") if show_progress else f
        for ln in iterator:
            if not ln.strip() or ln.startswith("#"):
                continue
            t = ln.rstrip("\n").split("\t")
            if len(t) != 9:
                continue
            ftype, attrs = t[2], t[8]

            if ftype in ("mRNA", "transcript"):
                mid = get_attr(attrs, "ID")
                if mid:
                    mrnas.add(mid)
                    par = get_attr(attrs, "Parent")
                    if par:
                        g = par.split(",")[0].strip()
                        if g:
                            gene_to_mrnas[g].add(mid)

            elif ftype == "CDS":
                parents = get_attr(attrs, "Parent")
                if parents:
                    for pid in parents.split(","):
                        pid = pid.strip()
                        if pid:
                            mrnas_with_cds.add(pid)
    finally:
        f.close()

    missing_list = sorted([m for m in mrnas if m not in mrnas_with_cds])
    missing = set(missing_list)

    # Build log
    log_lines.append(f"Total mRNA-like features: {len(mrnas)}")
    log_lines.append(f"mRNA-like with >=1 CDS:   {len(mrnas_with_cds & mrnas)}")
    log_lines.append(f"mRNA-like missing CDS:    {len(missing)}")
    if missing:
        log_lines.append("IDs missing CDS:")
        if max_ids_in_log is None:
            for mid in missing_list:
                log_lines.append(f"  {mid}")
        else:
            for mid in missing_list[: max_ids_in_log]:
                log_lines.append(f"  {mid}")
            if len(missing_list) > max_ids_in_log:
                log_lines.append(f"  ... and {len(missing_list) - max_ids_in_log} more")

    else:
        log_lines.append("No mRNA-like features missing CDS were found.")

    return mrnas, missing, gene_to_mrnas, log_lines


def compute_drop_genes(missing_mrnas: Set[str], gene_to_mrnas: Dict[str, Set[str]]) -> Set[str]:
    """Genes to drop are those for which all mRNA-like children are removed."""
    drop: Set[str] = set()
    for gid, mids in gene_to_mrnas.items():
        remaining = [m for m in mids if m not in missing_mrnas]
        if len(mids) > 0 and len(remaining) == 0:
            drop.add(gid)
    return drop


# ────────────────────────────────────────────────────────────────────────────────
# Pass 2: write filtered GFF3
# ────────────────────────────────────────────────────────────────────────────────
def pass2_filter(
    input_path: str,
    output_path: str,
    missing_mrnas: Set[str],
    drop_genes: Set[str],
    show_progress: bool,
) -> Tuple[int, int, int]:
    """
    Writes filtered GFF3 to output_path.
    Returns (dropped_mrnas, dropped_genes, dropped_children).
    """
    d_mrna = d_gene = d_child = 0

    inp = open_any(input_path, "rt")
    out = open_any(output_path, "wt")
    try:
        iterator = tqdm_wrap(inp, desc="Pass2 filter", unit="line") if show_progress else inp
        for ln in iterator:
            if not ln.strip():
                out.write(ln)
                continue
            if ln.startswith("#"):
                out.write(ln)
                continue

            t = ln.rstrip("\n").split("\t")
            if len(t) != 9:
                out.write(ln)
                continue

            ftype, attr_str = t[2], t[8]
            attrs = parse_attrs(attr_str)

            if ftype in ("mRNA", "transcript"):
                mid = attrs.get("ID", "")
                if mid in missing_mrnas:
                    d_mrna += 1
                    continue

            elif ftype == "gene":
                gid = attrs.get("ID", "")
                if gid in drop_genes:
                    d_gene += 1
                    continue

            else:
                parents = attrs.get("Parent", "")
                if parents:
                    par_list = [p.strip() for p in parents.split(",") if p.strip()]
                    if any((p in missing_mrnas) or (p in drop_genes) for p in par_list):
                        d_child += 1
                        continue

            out.write(ln)
    finally:
        out.close()
        inp.close()

    return d_mrna, d_gene, d_child


# ────────────────────────────────────────────────────────────────────────────────
# CLI
# ────────────────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description="Find mRNA/transcript features missing CDS, log summary, and remove them + orphaned genes."
    )
    ap.add_argument("-i", "--input", required=True, help="Input GFF3 (.gff3 or .gff3.gz)")
    ap.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output GFF3 path (written only if removals occur unless --always-write)",
    )
    ap.add_argument(
        "-l",
        "--log",
        default=None,
        help="Log file path (default: <input>.missingCDS.log)",
    )
    ap.add_argument(
        "--always-write",
        action="store_true",
        help="If nothing to remove, still copy input to --output (required).",
    )
    ap.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bars even if tqdm is installed.",
    )
    ap.add_argument(
        "--max-ids-in-log",
        type=int,
        default=500,
        help="Max IDs to list explicitly in the log for missing mRNAs (default: 500, set 0 for none, -1 for unlimited).",
    )
    args = ap.parse_args()

    input_path = args.input
    output_path = args.output
    show_progress = (not args.no_progress) and TQDM

    # Resolve log path
    if args.log:
        log_path = args.log
    else:
        base = os.path.basename(input_path)
        log_path = os.path.join(os.path.dirname(input_path), f"{base}.missingCDS.log")

    # Normalize max-ids-in-log
    if args.max_ids_in_log is None or args.max_ids_in_log < 0:
        max_ids = None  # unlimited
    else:
        max_ids = args.max_ids_in_log

    # Pass 1
    mrnas, missing_mrnas, gene_to_mrnas, log_lines = pass1_collect(
        input_path=input_path,
        show_progress=show_progress,
        max_ids_in_log=max_ids,
    )

    if missing_mrnas:
        drop_genes = compute_drop_genes(missing_mrnas, gene_to_mrnas)
        log_lines.append(f"Genes that will be dropped (fully orphaned): {len(drop_genes)}")

        if not output_path:
            # Need an output to write the filtered file
            with open(log_path, "w") as lf:
                lf.write("\n".join(log_lines) + "\n")
            sys.stderr.write(
                f"[ERROR] Removals required but no --output provided. Log written to {log_path}\n"
            )
            sys.exit(2)

        # Pass 2 filter
        d_mrna, d_gene, d_child = pass2_filter(
            input_path=input_path,
            output_path=output_path,
            missing_mrnas=missing_mrnas,
            drop_genes=drop_genes,
            show_progress=show_progress,
        )
        log_lines.append(f"Dropped mRNAs/transcripts: {d_mrna}")
        log_lines.append(f"Dropped genes:             {d_gene}")
        log_lines.append(f"Dropped child features:    {d_child}")
        log_lines.append(f"Output written to:         {output_path}")

        with open(log_path, "w") as lf:
            lf.write("\n".join(log_lines) + "\n")
        sys.stderr.write(
            f"[DONE] Removed {d_mrna} mRNAs, {d_gene} genes, {d_child} children. Log: {log_path}\n"
        )
        sys.exit(0)

    # No missing CDS; optionally copy through
    if output_path and args.always_write:
        with open_any(input_path, "rt") as src, open_any(output_path, "wt") as dst:
            for ln in src:
                dst.write(ln)
        log_lines.append(
            "No mRNA-like features missing CDS were found; input copied to output (--always-write)."
        )
        log_lines.append(f"Output copied to: {output_path}")
    else:
        log_lines.append("No mRNA-like features missing CDS were found; no output written.")

    with open(log_path, "w") as lf:
        lf.write("\n".join(log_lines) + "\n")
    sys.stderr.write(f"[DONE] Nothing to remove. Log: {log_path}\n")
    sys.exit(0)


if __name__ == "__main__":
    main()
