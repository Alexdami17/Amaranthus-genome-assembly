#!/usr/bin/env python3

"""
Combine two coordinate-sorted GFF3 files:
- Insert file B (gff3 from tRNAscan-SE-2.0.12) into file A (gff3) while preserving A's original line order (no resorting of A).
- Transform file B:
  * col3 'pseudogene' -> 'tRNA'
  * Prefix all ID/Parent values with a detected prefix from file A (or --prefix)
  * Also prefix the leading token in Name= (before the first '-') to match IDs
  * Remove exactly one stray trailing ';' from ID/Parent/Name values (not separators)
  * Rename exonN -> exon.N inside ID/Parent (and defensively elsewhere)
  * Model tRNA as transcript-level features under a gene container
  * Set column 2 (source) to match A's dominant source by default
- Preserve A's pragma/comment lines; write a single ##gff-version 3 header.

Progress:
- Always shows progress. Uses tqdm if installed, else a lightweight stderr counter.
"""

import argparse
import sys
import re
from collections import Counter, defaultdict
from typing import Optional, List, Tuple, Iterator

GFF_VERSION_LINE = "##gff-version 3"

# ------------------------ Progress Helpers ------------------------

_warned_tqdm = False

class _FallbackPbar:
    def __init__(self, total: int, desc: str):
        self.total = total
        self.desc = desc
        self.n = 0
        self._step = max(1, total // 100) if total else 10000

    def update(self, n=1):
        self.n += n
        if self.n % self._step == 0 or self.n >= self.total:
            pct = f"{(self.n / self.total * 100):5.1f}%" if self.total else "N/A"
            sys.stderr.write(f"\r[{self.desc}] {self.n}/{self.total} ({pct})")
            sys.stderr.flush()

    def close(self):
        sys.stderr.write("\n")
        sys.stderr.flush()

def _new_pbar(total: int, desc: str):
    global _warned_tqdm
    try:
        from tqdm import tqdm  # type: ignore
        return tqdm(total=total, desc=desc, unit="lines", leave=False)
    except Exception:
        if not _warned_tqdm:
            sys.stderr.write("Note: install tqdm (`pip install tqdm`) for a nicer progress bar.\n")
            _warned_tqdm = True
        return _FallbackPbar(total, desc)

def _count_data_lines(path: str) -> int:
    c = 0
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            c += 1
    return c

# ------------------------ CLI & Core ------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Insert (merge) a coordinate-sorted GFF3 from file B (tRNAscan-SE-2.0.12) into file A (gff3) with specific renames/fixes."
    )
    p.add_argument("-a", "--gffA", required=True, help="Primary GFF3, already coordinate-sorted")
    p.add_argument("-b", "--gffB", required=True, help="Secondary GFF3 (to be transformed) and merged into A")
    p.add_argument("-o", "--out", required=True, help="Output GFF3 path")
    p.add_argument("--prefix", default=None,
                   help="Prefix to add to all ID/Parent values from file B (e.g., 'AmaTu_'). "
                        "If omitted, the script tries to infer it from file A.")
    p.add_argument("--rename-pseudogene-to", default="tRNA",
                   help='For file B: rename feature type "pseudogene" in column 3 to this (default: tRNA)')
    p.add_argument("--infer-source-from", choices=["A", "B", "auto"], default="auto",
                   help="Column 2 for B rows: match A’s most common source (auto/A) or keep B’s (B). Default: auto")
    return p.parse_args()

def is_header(line: str) -> bool:
    return line.startswith("#")

# ---- ADDED: detect_common_source_A ----
def detect_common_source_A(pathA: str):
    """
    Return the most common value from column 2 (source) in a GFF3 file,
    ignoring header/pragma lines and malformed rows.
    """
    counter = Counter()
    with open(pathA, "r", encoding="utf-8") as f:
        for line in f:
            if is_header(line):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 9:
                counter[parts[1]] += 1
    return counter.most_common(1)[0][0] if counter else None
# ---- END ADDED ----

def parse_gff_line(line: str):
    parts = line.rstrip("\n").split("\t")
    if len(parts) != 9:
        raise ValueError(f"Invalid GFF3 (not 9 columns): {line.strip()}")
    parts[3] = int(parts[3])
    parts[4] = int(parts[4])
    return parts  # [seqid, source, type, start, end, score, strand, phase, attrs]

def split_attrs(attr_str: str) -> List[Tuple[str, str]]:
    if attr_str.strip() in (".", ""):
        return []
    items = attr_str.split(";")
    out = []
    for item in items:
        if item == "":
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out.append((k, v))
        else:
            out.append((item, ""))
    return out

def join_attrs(kvs: List[Tuple[str, str]]) -> str:
    if not kvs:
        return "."
    return ";".join(f"{k}={v}" if v != "" else k for k, v in kvs)

def remove_one_trailing_semicolon(value: str) -> str:
    return value[:-1] if value.endswith(";") else value

def _fix_exon_suffix(val: str) -> str:
    # Normalize "exon1" -> "exon.1" and ".exon1" -> ".exon.1"
    val = re.sub(r"\.exon(\d+)", r".exon.\1", val)
    val = re.sub(r"(^|[^.])exon(\d+)", lambda m: m.group(1) + f"exon.{m.group(2)}", val)
    return val

def fix_name_token(val: str, prefix: str) -> str:
    """Prefix the leading token in Name= (before the first '-') if not already prefixed."""
    v = remove_one_trailing_semicolon(val)
    if not prefix:
        return v
    parts = v.split("-", 1)
    lead = parts[0]
    tail = parts[1] if len(parts) > 1 else ""
    if not lead.startswith(prefix):
        lead = prefix + lead
    return lead if tail == "" else f"{lead}-{tail}"

def fix_id_like_token(val: str, prefix: str) -> str:
    v = remove_one_trailing_semicolon(val)
    if prefix and not v.startswith(prefix):
        v = prefix + v
    return _fix_exon_suffix(v)

def transform_attrs_for_B(attr_str: str, prefix: str) -> str:
    """Transform general attributes for B lines (ID/Parent prefixed; Name leading token prefixed; exonN normalized)."""
    kvs = split_attrs(attr_str)
    new = []
    for k, v in kvs:
        v = remove_one_trailing_semicolon(v)
        if k == "ID":
            v = fix_id_like_token(v, prefix)
            new.append((k, v))
        elif k == "Parent":
            parents = [p.strip() for p in v.split(",") if p.strip() != ""]
            parents = [fix_id_like_token(p, prefix) for p in parents]
            new.append((k, ",".join(parents)))
        elif k == "Name":
            new.append((k, fix_name_token(v, prefix)))
        else:
            new.append((k, _fix_exon_suffix(v)))
    return join_attrs(new)

def attrs_to_kv_list(attr_str: str) -> List[Tuple[str, str]]:
    kvs = split_attrs(attr_str)
    return [(k, remove_one_trailing_semicolon(v)) for k, v in kvs]

def attrs_get(kvs: List[Tuple[str, str]], key: str) -> Optional[str]:
    for k, v in kvs:
        if k == key:
            return v
    return None

def attrs_set(kvs: List[Tuple[str, str]], key: str, value: str):
    for i, (k, _) in enumerate(kvs):
        if k == key:
            kvs[i] = (key, value)
            return
    kvs.append((key, value))

# ---------- Helpers for ID pattern ----------
_trna_fallback_counter = defaultdict(int)   # per-seqid fallback numbering

def _extract_trna_number(tid: str, name: Optional[str]) -> Optional[int]:
    """Find N from '*trnaN' / '*tRNA N' etc in ID or Name (case-insensitive)."""
    for s in (tid, name or ""):
        m = re.search(r'(?i)(?:^|[._-])trna\s*(\d+)\b', s)
        if m:
            try:
                return int(m.group(1))
            except Exception:
                pass
    return None

def _build_gene_name(existing_name: Optional[str], prefix: str, seqid: str, n: int, kv: List[Tuple[str,str]]) -> str:
    """Gene Name leading token should be '<prefix><seqid>.tRNA<n>' then keep any suffix (e.g. '-HisGTG')."""
    base = f"{prefix}{seqid}.tRNA{n}"
    if existing_name:
        parts = existing_name.split("-", 1)
        tail = parts[1] if len(parts) > 1 else ""
        return base if not tail else f"{base}-{tail}"
    # If no Name provided, try to build a nice one from isotype/anticodon
    iso = attrs_get(kv, "isotype") or ""
    ac  = attrs_get(kv, "anticodon") or ""
    suffix = f"-{iso}{ac}" if (iso or ac) else ""
    return base + suffix

def make_gene_and_trna_attrs_from_trna_line(attrs_str: str, prefix: str, seqid: str) -> Tuple[str, str, str]:
    """
    From a tRNA line's (already prefix-transformed) attrs, build:
    - gene_attrs: ID=<prefix><seqid>tRNA<n>, Name=<prefix><seqid>.tRNA<n>-<isotype><anticodon> (if present), keep isotype/anticodon/gene_biotype; drop Parent
    - trna_attrs: ID=<prefix><seqid>.trna<n>, Parent=<geneID> (minimal set)
    Returns tuple (gene_id, gene_attrs_str, trna_attrs_str)
    """
    kv = attrs_to_kv_list(attrs_str)
    tid = attrs_get(kv, "ID")
    if tid is None:
        raise ValueError("tRNA line missing ID attribute after transformation.")

    name_in = attrs_get(kv, "Name")
    n = _extract_trna_number(tid, name_in)
    if n is None:
        # fallback per seqid
        _trna_fallback_counter[seqid] += 1
        n = _trna_fallback_counter[seqid]

    gene_id = f"{prefix}{seqid}tRNA{n}"   # no dot, tRNA caps
    trna_id = f"{prefix}{seqid}.trna{n}"  # with dot, trna lower

    # Build gene attrs
    gene_name = _build_gene_name(name_in, prefix, seqid, n, kv)
    gene_kv: List[Tuple[str, str]] = [("ID", gene_id), ("Name", gene_name)]
    for k, v in kv:
        if k in ("ID", "Parent", "Name"):
            continue
        gene_kv.append((k, v))
    gene_attrs = join_attrs(gene_kv)

    # Build tRNA (transcript) attrs: keep it minimal and clean
    trna_kv: List[Tuple[str, str]] = [("ID", trna_id), ("Parent", gene_id)]
    trna_attrs = join_attrs(trna_kv)

    return gene_id, gene_attrs, trna_attrs

# ---- Helpers to infer prefix & collect pragmas ----
def _first_id_from_attrs(attrs: str) -> Optional[str]:
    for k, v in split_attrs(attrs):
        if k == "ID":
            return v
    return None

def detect_prefix_from_A(pathA: str) -> Optional[str]:
    """
    Infer a prefix by finding the substring of an ID that precedes the seqid,
    e.g. ID='AmaAk_Chr01Xg000010' with seqid='Chr01X' -> prefix 'AmaAk_'.
    Returns None if no stable prefix is detected.
    """
    candidates = Counter()
    with open(pathA, "r", encoding="utf-8") as f:
        for line in f:
            if is_header(line):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, _, _, _, _, _, _, _, attrs = parts
            idv = _first_id_from_attrs(attrs)
            if not idv:
                continue
            idv = remove_one_trailing_semicolon(idv)
            pos = idv.find(seqid)
            if pos > 0:
                prefix = idv[:pos]
                if 1 <= len(prefix) <= 20:
                    candidates[prefix] += 1
    if not candidates:
        return None
    prefix, count = candidates.most_common(1)[0]
    if len(candidates) > 1:
        total = sum(candidates.values())
        if count / total < 0.7:
            return None
    return prefix

def read_pragmas_until_fasta(pathA: str) -> List[str]:
    pragmas = []
    with open(pathA, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith("##") and not line.startswith(GFF_VERSION_LINE):
                pragmas.append(line.rstrip("\n"))
    return pragmas
# ---- END ADDED ----

# ------------------------ Streaming B (transformed) ------------------------

def iter_transformed_B(path: str, source_for_B: Optional[str], rename_pseudogene_to: str, prefix: str,
                       total_lines: int) -> Iterator[Tuple[str, int, str]]:
    """
    Yields transformed B lines as (seqid, start, line) in B's order,
    expanding each tRNA into: gene (.) then tRNA (Parent=gene). Exons stay as-is
    (with ID/Parent/Name fixes) and follow wherever they appear in B.
    """
    pbar = _new_pbar(total_lines, "Reading/transforming B")
    # exon numbering per tRNA transcript ID (so exon IDs become <trnaID>.exon.N)
    exon_count = defaultdict(int)

    try:
        with open(path, "r", encoding="utf-8") as f:
            for raw in f:
                if is_header(raw):
                    continue
                seqid, source, ftype, start, end, score, strand, phase, attrs = parse_gff_line(raw)

                # column-3 rename for B
                if ftype == "pseudogene":
                    ftype = rename_pseudogene_to

                # harmonize source (col 2)
                if source_for_B is not None:
                    source = source_for_B

                # transform attributes: ID/Parent/Name prefixing + exonN normalization
                attrs = transform_attrs_for_B(attrs, prefix=prefix or "")
                kv = attrs_to_kv_list(attrs)

                if ftype == "tRNA":
                    gene_id, gene_attrs, trna_attrs = make_gene_and_trna_attrs_from_trna_line(attrs, prefix or "", seqid)
                    gene_line = "\t".join(map(str, [seqid, source, "gene", start, end, ".", strand, ".", gene_attrs]))
                    yield (seqid, start, gene_line)

                    trna_line = "\t".join(map(str, [seqid, source, "tRNA", start, end, score if score != "." else ".", strand, ".", trna_attrs]))
                    yield (seqid, start, trna_line)

                elif ftype == "exon":
                    # Re-ID exon as <trnaID>.exon.N where trnaID is the (single) Parent we expect for tRNAs
                    parent_val = attrs_get(kv, "Parent") or ""
                    parents = [p.strip() for p in parent_val.split(",") if p.strip()]
                    if parents:
                        p0 = parents[0]
                        n = _extract_trna_number(p0, attrs_get(kv, "Name"))
                        if n is not None:
                            trna_id = f"{prefix}{seqid}.trna{n}"
                        else:
                            trna_id = p0  # fallback: use as-is
                        exon_count[trna_id] += 1
                        exon_id = f"{trna_id}.exon.{exon_count[trna_id]}"
                        attrs_set(kv, "ID", exon_id)
                        parents[0] = trna_id
                        attrs_set(kv, "Parent", ",".join(parents))
                        attrs = join_attrs(kv)
                    out_line = "\t".join(map(str, [seqid, source, ftype, start, end, score, strand, phase, attrs]))
                    yield (seqid, start, out_line)

                else:
                    out_line = "\t".join(map(str, [seqid, source, ftype, start, end, score, strand, phase, attrs]))
                    yield (seqid, start, out_line)

                pbar.update(1)
    finally:
        pbar.close()

# ------------------------ Main (streaming merge, A-stable; no warnings) ------------------------

def main():
    args = parse_args()

    # Determine source for B
    source_for_B = None
    if args.infer_source_from in ("A", "auto"):
        srcA = detect_common_source_A(args.gffA)
        if srcA:
            source_for_B = srcA

    # Determine prefix to apply to B's ID/Parent/Name
    prefix = args.prefix
    if prefix is None:
        prefix = detect_prefix_from_A(args.gffA)
        if prefix is None:
            sys.stderr.write(
                "ERROR: Could not infer an ID prefix from file A. "
                "Please supply one explicitly with --prefix (e.g., --prefix AmaTu_).\n"
            )
            sys.exit(2)

    # Read pragmas/comments from A (except version line), up to FASTA
    pragmas = read_pragmas_until_fasta(args.gffA)

    # Count lines for progress bars (always shown)
    total_A = _count_data_lines(args.gffA)
    total_B_raw = _count_data_lines(args.gffB)

    # Prepare transformed B iterator (preserves B order)
    b_iter = iter_transformed_B(args.gffB, source_for_B, args.rename_pseudogene_to, prefix, total_B_raw)
    b_next = next(b_iter, None)  # (seqid, start, line) or None

    with open(args.out, "w", encoding="utf-8") as w:
        # Header + pragmas
        w.write(f"{GFF_VERSION_LINE}\n")
        for line in pragmas:
            w.write(line + "\n")

        # Merge: write A as-is, inserting B lines whose (seqid, start) are strictly before the current A line
        pbarA = _new_pbar(total_A, "Writing A (+inserting B)")
        try:
            with open(args.gffA, "r", encoding="utf-8") as fa:
                for rawA in fa:
                    if is_header(rawA):
                        continue
                    a_parts = parse_gff_line(rawA)
                    a_seqid, _, _, a_start, _, _, _, _, _ = a_parts

                    while b_next is not None and (
                        (b_next[0] < a_seqid) or (b_next[0] == a_seqid and b_next[1] < a_start)
                    ):
                        w.write(b_next[2].rstrip("\n") + "\n")
                        b_next = next(b_iter, None)

                    w.write(rawA.rstrip("\n") + "\n")
                    pbarA.update(1)

            # After all A lines, flush remaining B
            pbarFlush = _new_pbar(0, "Flushing remaining B")
            while b_next is not None:
                w.write(b_next[2].rstrip("\n") + "\n")
                b_next = next(b_iter, None)
                pbarFlush.update(1)
            pbarFlush.close()

        finally:
            pbarA.close()

if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        try:
            sys.stderr.close()
        except Exception:
            pass
        sys.exit(0)
