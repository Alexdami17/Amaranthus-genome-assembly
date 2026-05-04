#!/usr/bin/env python3
"""
filter_repeatmasker_align.py

Filter RepeatMasker .align files by query sequence (chromosome/contig).

An "alignment record" in .align starts with a header line like:
  876 5.89 0.00 1.68 Chr01X 2 122 (53302002) C rnd-4_family-930#Unknown ...

and includes subsequent alignment/matrix/kimura lines until the next header.

This script streams line-by-line, preserves formatting, and supports:
- exact chromosome list match (--chr)
- regex match (--chr-regex)
- invert match (--invert)
- optional split output into one file per chr (--split-outdir)
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from typing import Dict, Optional, TextIO, Tuple


HEADER_RE = re.compile(
    r"""^\s*
    (?P<score>\d+)\s+
    (?P<pdiv>\d+(?:\.\d+)?)\s+
    (?P<pdel>\d+(?:\.\d+)?)\s+
    (?P<pins>\d+(?:\.\d+)?)\s+
    (?P<qseq>\S+)
    (?:\s+|$)
    """,
    re.VERBOSE,
)


def parse_header_qseq(line: str) -> Optional[str]:
    """
    Return query sequence/chr name if line is a header, else None.
    """
    m = HEADER_RE.match(line)
    if not m:
        return None
    return m.group("qseq")


def decide_keep(qseq: str, chr_set: Optional[set], chr_regex: Optional[re.Pattern], invert: bool) -> bool:
    matched = False
    if chr_set is not None:
        matched = qseq in chr_set
    if chr_regex is not None:
        matched = bool(chr_regex.search(qseq)) or matched

    return (not matched) if invert else matched


def get_writer(
    qseq: str,
    split_outdir: Optional[str],
    writers: Dict[str, TextIO],
    default_out: Optional[TextIO],
) -> TextIO:
    """
    Choose output handle:
    - if split_outdir: one file per qseq
    - else: default_out must be set
    """
    if split_outdir:
        if qseq not in writers:
            os.makedirs(split_outdir, exist_ok=True)
            outpath = os.path.join(split_outdir, f"{qseq}.align")
            writers[qseq] = open(outpath, "w", encoding="utf-8")
        return writers[qseq]
    if default_out is None:
        raise RuntimeError("Internal error: default_out is None but split_outdir not set.")
    return default_out


def flush_block(block_lines: list[str], keep: bool, qseq: Optional[str],
                split_outdir: Optional[str], writers: Dict[str, TextIO], out_fh: Optional[TextIO]) -> None:
    if not block_lines:
        return
    if not keep:
        return
    if qseq is None:
        # Should not happen if blocks begin with a header, but be safe.
        writer = out_fh
        if writer is None:
            return
        writer.writelines(block_lines)
        return

    writer = get_writer(qseq, split_outdir, writers, out_fh)
    writer.writelines(block_lines)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Filter RepeatMasker .align records by query sequence (chromosome/contig)."
    )
    ap.add_argument("-i", "--input", required=True, help="Input RepeatMasker .align file")
    ap.add_argument("-o", "--output", help="Output file (default: stdout). Ignored if --split-outdir is used.")
    ap.add_argument("--chr", nargs="+", help="Exact query sequence names to keep (e.g., Chr01X Chr02)")
    ap.add_argument("--chr-regex", help=r"Regex to match query sequence names (e.g., '^Chr0(1|2)')")
    ap.add_argument("--invert", action="store_true", help="Invert match (exclude specified chromosomes)")
    ap.add_argument("--split-outdir", help="Write one output .align per query sequence into this directory")
    args = ap.parse_args()

    if not args.chr and not args.chr_regex:
        ap.error("You must provide at least one of --chr or --chr-regex")

    if args.output and args.split_outdir:
        print("Note: --output is ignored when --split-outdir is used.", file=sys.stderr)

    chr_set = set(args.chr) if args.chr else None
    chr_regex = re.compile(args.chr_regex) if args.chr_regex else None

    out_fh: Optional[TextIO] = None
    if not args.split_outdir:
        if args.output:
            out_fh = open(args.output, "w", encoding="utf-8")
        else:
            out_fh = sys.stdout

    writers: Dict[str, TextIO] = {}

    current_block: list[str] = []
    current_qseq: Optional[str] = None
    current_keep: bool = False
    seen_first_header = False

    with open(args.input, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            qseq = parse_header_qseq(line)

            if qseq is not None:
                # New record begins. Flush previous record.
                if seen_first_header:
                    flush_block(current_block, current_keep, current_qseq, args.split_outdir, writers, out_fh)
                else:
                    # If there is any preamble before first header, we drop it by default.
                    seen_first_header = True

                # Start new record
                current_block = [line]
                current_qseq = qseq
                current_keep = decide_keep(qseq, chr_set, chr_regex, args.invert)
            else:
                # Continuation of current record (or preamble before first header)
                if seen_first_header:
                    current_block.append(line)

        # EOF flush
        if seen_first_header:
            flush_block(current_block, current_keep, current_qseq, args.split_outdir, writers, out_fh)

    # Close handles
    for w in writers.values():
        w.close()
    if out_fh not in (None, sys.stdout):
        out_fh.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
