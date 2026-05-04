#!/usr/bin/env python3

import argparse
import csv
import logging
from collections import defaultdict
import re
import sys  

logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

NULL_TOKENS = {"", "na", "nan", "none", "-", ".", "null"}

# ------------------------ Progress Helpers ------------------------
_WARNED_TQDM = False

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
    """Use tqdm if available; otherwise warn once and switch to a lightweight stderr counter."""
    global _WARNED_TQDM
    try:
        from tqdm import tqdm  # type: ignore
        return tqdm(total=total, desc=desc, unit="lines", leave=False)
    except Exception:
        if not _WARNED_TQDM:
            sys.stderr.write("Warning: install tqdm (`pip install tqdm`) for a nicer progress bar.\n")
            _WARNED_TQDM = True
        return _FallbackPbar(total, desc)

def _count_data_lines_gff(path: str) -> int:
    c = 0
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c += 1
    return c

def _count_rows_tsv(path: str) -> int:
    # count total lines minus header
    n = 0
    with open(path, "r", newline="") as f:
        for _ in f:
            n += 1
    return max(0, n - 1)
# ---------------------- END: Progress Helpers ----------------------


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate enriched GFF3 from EnTAP TSV and PASA-filtered GFF3.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    req = parser.add_argument_group("required arguments")
    req.add_argument("-i", "--input", required=True, help="Input GFF3 file")
    req.add_argument("-t", "--tsv", required=True, help="TSV file with attributes")
    req.add_argument("-o", "--output", required=True, help="Output GFF3 file")
    req.add_argument("-m", "--map", required=True, help="Output mapping file")
    req.add_argument("--id", required=True, help="Base name for gene IDs, e.g., AmaTu")

    opt = parser.add_argument_group("optional arguments")
    opt.add_argument("-s", "--source",default=".",help="Source string for column 2 of GFF3 output",
    )

    return parser.parse_args()

def safe_get(row, idx, default=""):
    if idx < 0:
        return default
    try:
        val = row[idx].strip()
    except IndexError:
        return default
    return val

def is_null(val: str) -> bool:
    return val.strip().lower() in NULL_TOKENS

def load_tsv(tsv_file):
    """
    Return dict keyed by original mRNA ID with:
      dbxref  -> accession (first token before space) OR "" (treated as missing)
      note    -> remainder of that field after accession (can be "")
      go      -> list of GO terms (one per GO column, first GO:* only, de-duped)
      interpro-> raw string from InterPro column (column 52)
    """
    mapping = {}

    # TSV progress
    total_rows = _count_rows_tsv(tsv_file)
    pbar = _new_pbar(total_rows, "Reading TSV")

    with open(tsv_file, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, [])
        col_orig_id = 0
        col_dbxref = 12
        col_go_terms = [43, 44, 45]
        col_interpro = 51

        for row in reader:
            orig_id = safe_get(row, col_orig_id)
            if is_null(orig_id):
                pbar.update(1)
                continue

            dbxref_field = safe_get(row, col_dbxref)
            dbxref = ""
            note_from_dbxref = ""
            if not is_null(dbxref_field):
                toks = dbxref_field.split()
                dbxref = toks[0]
                if len(toks) > 1:
                    note_from_dbxref = " ".join(toks[1:])

            go_terms = []
            seen_go = set()
            for i in col_go_terms:
                colval = safe_get(row, i)
                if is_null(colval):
                    continue
                first = None
                for piece in re.split(r"[,\s]+", colval):
                    if piece.startswith("GO:"):
                        first = piece
                        break
                if first and first not in seen_go:
                    go_terms.append(first)
                    seen_go.add(first)

            interpro_raw = safe_get(row, col_interpro)
            interpro_raw = "" if is_null(interpro_raw) else interpro_raw

            mapping[orig_id] = {
                "dbxref": dbxref,                    # "" if missing
                "note": note_from_dbxref.strip(),    # may be ""
                "go": go_terms,                      # list
                "interpro": interpro_raw,            # raw; may be ""
            }

            pbar.update(1)

    pbar.close()
    return mapping

def parse_attributes(attr_string):
    attrs = {}
    if not attr_string:
        return attrs
    for part in attr_string.split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key] = value
    return attrs

def format_attributes(attrs):
    # Normalize key order for readability; include any extra keys at the end
    order = ["ID", "Name", "Parent", "Dbxref", "Ontology_term", "Note"]
    keys = [k for k in order if k in attrs] + [k for k in attrs.keys() if k not in order]
    return ";".join(f"{k}={attrs[k]}" for k in keys)

def interpro_from_raw(raw: str):
    """
    Return:
      ':(InterPro|IPRxxxxx|<desc>)' using the first IPR found and everything inside its matching ()
      ':(InterPro|Unknown|Unknown)' only if the entire cell is exactly '-(-)'
      None if cell is NaN/blank or no IPR is present
    Handles commas and inner parentheses inside the description.
    """
    if not raw or is_null(raw):
        return None

    s = raw.strip()
    if s == "-(-)":
        return ":(InterPro|Unknown|Unknown)"

    # Find first IPR id followed by '('
    m = re.search(r"(IPR\d+)\s*\(", s)
    if not m:
        return None

    ipr = m.group(1)
    open_i = s.find("(", m.start())  # index of the '(' after the IPR
    # Match the corresponding closing ')' with a small parenthesis counter
    depth = 0
    close_i = None
    for i, ch in enumerate(s[open_i:], start=open_i):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                close_i = i
                break

    # Fallback if we never find a balanced ')'
    if close_i is None:
        j = s.find(")", open_i + 1)
        desc = s[open_i + 1 : (j if j != -1 else len(s))].strip()
    else:
        desc = s[open_i + 1 : close_i].strip()

    if not desc:
        desc = "Unknown"

    return f":(InterPro|{ipr}|{desc})"

def enrich_gene_attrs_from_tsv(attrs, tsv_row):
    """
    Apply:
      - Dbxref if present (non-NaN). If missing/NaN -> omit Dbxref AND set Note=Unknown Protein.
      - Ontology_term from GO list (comma-joined)
      - Note: start with existing Note or dbxref remainder; if Dbxref missing, ensure it starts with 'Unknown Protein'.
        Append InterPro suffix per rules (None/-(-)/IPRxxx(desc)).
    """
    existing_note = attrs.get("Note", "").strip()
    dbx = tsv_row.get("dbxref", "").strip()
    dbxref_present = bool(dbx)

    # GO terms
    if tsv_row.get("go"):
        attrs["Ontology_term"] = ",".join(tsv_row["go"])

    # Dbxref handling + base note
    if dbxref_present:
        attrs["Dbxref"] = dbx
        base_note = tsv_row.get("note", "").strip() or existing_note
    else:
        # Omit Dbxref entirely and force Unknown Protein
        if "Dbxref" in attrs:
            del attrs["Dbxref"]
        base_note = "Unknown Protein"

    # InterPro suffix
    ip_suffix = interpro_from_raw(tsv_row.get("interpro", ""))  # None / ':(InterPro|...)'
    if ip_suffix:
        note_full = (base_note + ip_suffix) if base_note else ip_suffix.lstrip(":")
    else:
        note_full = base_note

    if note_full:
        attrs["Note"] = note_full

    return attrs

def main():
    args = parse_args()
    tsv_map = load_tsv(args.tsv)

    # Collect features grouped by chromosome, preserving input order
    features_by_chrom = defaultdict(list)

    # GFF3 read progress
    total_gff = _count_data_lines_gff(args.input)
    pbar_in = _new_pbar(total_gff, "Reading GFF3")

    with open(args.input) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                logging.warning("Skipping malformed GFF3 line (expected 9 cols): %s", line.strip())
                continue
            features_by_chrom[parts[0]].append(parts)
            pbar_in.update(1)

    pbar_in.close()

    # Gene IDs increment by 10
    global_gene_counter = 10

    # Order by numeric part (Chr01, Chr02, …) regardless of trailing letters; then ChrC, ChrM; then others.
    all_seqids = list(features_by_chrom.keys())

    # Group seqids by numeric part (any width), keep first-seen order within each numeric group
    first_seen_index = {s: i for i, s in enumerate(all_seqids)}
    from collections import defaultdict as _dd
    numeric_groups = _dd(list)

    for s in all_seqids:
        m = re.fullmatch(r"Chr(\d+)([A-Za-z].*)?\Z", s)
        if m:
            num = int(m.group(1))
            numeric_groups[num].append(s)

    preferred_order = []
    for num in sorted(numeric_groups):
        # FIX: use the actual dict name and key
        preferred_order.extend(sorted(numeric_groups[num], key=lambda x: first_seen_index[x]))

    # Add specials if present
    for special in ("ChrC", "ChrM"):
        if special in features_by_chrom:
            preferred_order.append(special)

    # Anything not covered above goes to extras; sort extras by numeric part if any, else by name
    def _nat_key(ch):
        m = re.fullmatch(r"Chr(\d+)([A-Za-z].*)?\Z", ch)
        if m:
            return (0, int(m.group(1)))
        if ch == "ChrC": return (1, 0)
        if ch == "ChrM": return (1, 1)
        return (2, ch)

    extras = sorted([c for c in features_by_chrom if c not in preferred_order], key=_nat_key)
    ordered_seqids = preferred_order + extras
    
    new_features = []
    mrna_map = {}             # original mRNA ID -> new mRNA ID
    mrna_seen_for_gene = {}   # gene_id -> bool

    for chrom in ordered_seqids:
        features = features_by_chrom[chrom]
        last_gene_id = None # <-- initialize per seqid

        for feature in features:
            seqid, source, ftype, start, end, score, strand, phase, attr_string = feature
            attrs = parse_attributes(attr_string)

            if ftype == "gene":
                gene_id = f"{args.id}_{chrom}g{global_gene_counter:06d}"
                global_gene_counter += 10
                last_gene_id = gene_id
                mrna_seen_for_gene[gene_id] = False
                attrs["ID"] = gene_id
                attrs["Name"] = gene_id
                # defer enrichment until we see mRNA (we'll patch this record in-place then)
                new_features.append((seqid, args.source, ftype, start, end, score, strand, phase, format_attributes(attrs)))
                continue

            if ftype in ("mRNA", "transcript"):
                # If mRNA appears before gene, synthesize gene using mRNA coords
                if not last_gene_id:
                    gene_id = f"{args.id}_{chrom}g{global_gene_counter:06d}"
                    global_gene_counter += 10
                    last_gene_id = gene_id
                    mrna_seen_for_gene[gene_id] = False
                    synth_gene_attrs = {"ID": gene_id, "Name": gene_id}
                    orig_mrna_id_tmp = attrs.get("ID", "")
                    if orig_mrna_id_tmp in tsv_map:
                        synth_gene_attrs = enrich_gene_attrs_from_tsv(synth_gene_attrs, tsv_map[orig_mrna_id_tmp])
                    new_features.append((seqid, args.source, "gene", start, end, score, strand, ".", format_attributes(synth_gene_attrs)))

                # Enforce single isoform per gene
                if mrna_seen_for_gene.get(last_gene_id, False):
                    logging.warning("Extra mRNA under gene %s on %s; dropping extra mRNA line.", last_gene_id, chrom)
                    orig_mrna_id = attrs.get("ID")
                    if not orig_mrna_id:
                        raise ValueError(
                            f"Input GFF3 {ftype} on {seqid}:{start}-{end} lacks an ID; cannot join to TSV col 1. "
                            "Please ensure mRNA/transcript features have IDs matching TSV column 1."
                        )
                    if orig_mrna_id:
                        mrna_map[orig_mrna_id] = f"{last_gene_id}-mRNA-1"
                    continue

                mrna_id = f"{last_gene_id}-mRNA-1"
                orig_mrna_id = attrs.get("ID")
                if orig_mrna_id:
                    mrna_map[orig_mrna_id] = mrna_id

                # Clean / set attributes
                attrs["ID"] = mrna_id
                attrs["Parent"] = last_gene_id
                if "Name" in attrs:
                    del attrs["Name"]

                # Enrich the most recent gene now that we know the mRNA’s TSV row
                if orig_mrna_id and orig_mrna_id in tsv_map:
                    for i in range(len(new_features) - 1, -1, -1):
                        if new_features[i][0] == chrom and new_features[i][2] == "gene":
                            g_seqid, g_source, g_ftype, g_start, g_end, g_score, g_strand, g_phase, g_attr = new_features[i]
                            g_attrs = parse_attributes(g_attr)
                            g_attrs = enrich_gene_attrs_from_tsv(g_attrs, tsv_map[orig_mrna_id])
                            new_features[i] = (g_seqid, g_source, g_ftype, g_start, g_end, g_score, g_strand, g_phase, format_attributes(g_attrs))
                            break

                new_features.append((seqid, args.source, ftype, start, end, score, strand, phase, format_attributes(attrs)))
                mrna_seen_for_gene[last_gene_id] = True
                continue

            # Other features (exon, CDS, UTRs, etc.)
            if "Parent" in attrs:
                parents = attrs["Parent"].split(",")
                new_parents = []
                for p in parents:
                    p = p.strip()
                    new_parents.append(mrna_map.get(p, p))
                attrs["Parent"] = ",".join(new_parents)

            if "ID" in attrs and last_gene_id:
                old_id = attrs["ID"]
                attrs["ID"] = re.sub(r"^[^\.]+", last_gene_id, old_id)

            new_features.append((seqid, args.source, ftype, start, end, score, strand, phase, format_attributes(attrs)))

    # ---- Strand-aware renumbering + grouped reordering under each mRNA ----
    mrna_strand = {}
    children_by_mrna = defaultdict(lambda: {
        "exon": [], "CDS": [], "five_prime_UTR": [], "three_prime_UTR": [], "UTR": []
    })

    for idx, (seqid, source, ftype, start, end, score, strand, phase, attr_string) in enumerate(new_features):
        attrs = parse_attributes(attr_string)

        if ftype in ("mRNA", "transcript"):
            mrna_id = attrs.get("ID")
            if mrna_id:
                mrna_strand[mrna_id] = strand

        elif ftype in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR", "UTR"):
            parent = attrs.get("Parent", "")
            if not parent:
                continue
            mrna_id = parent.split(",")[0].strip()
            if not mrna_id:
                continue
            try:
                s = int(start)
            except ValueError:
                s = 0
            children_by_mrna[mrna_id][ftype if ftype in ("five_prime_UTR","three_prime_UTR") else ftype].append((s, idx))

    for mrna_id, groups in children_by_mrna.items():
        strand = mrna_strand.get(mrna_id, "+")
        rev = (strand == "-")

        groups["exon"].sort(key=lambda x: x[0], reverse=rev)
        for rank, (_, feat_idx) in enumerate(groups["exon"], start=1):
            seqid, source, ftype0, start, end, score, strand0, phase, attr_string = new_features[feat_idx]
            attrs = parse_attributes(attr_string)
            attrs["ID"] = f"{mrna_id}.exon.{rank}"
            new_features[feat_idx] = (seqid, source, ftype0, start, end, score, strand0, phase, format_attributes(attrs))

        groups["CDS"].sort(key=lambda x: x[0], reverse=rev)
        for rank, (_, feat_idx) in enumerate(groups["CDS"], start=1):
            seqid, source, ftype0, start, end, score, strand0, phase, attr_string = new_features[feat_idx]
            attrs = parse_attributes(attr_string)
            attrs["ID"] = f"{mrna_id}.cds.{rank}"
            new_features[feat_idx] = (seqid, source, ftype0, start, end, score, strand0, phase, format_attributes(attrs))

        groups["five_prime_UTR"].sort(key=lambda x: x[0], reverse=rev)
        groups["three_prime_UTR"].sort(key=lambda x: x[0], reverse=rev)
        groups["UTR"].sort(key=lambda x: x[0], reverse=rev)

    consumed = set()
    rebuilt = []
    child_types = {"exon", "CDS", "five_prime_UTR", "three_prime_UTR", "UTR"}

    for idx, rec in enumerate(new_features):
        if idx in consumed:
            continue

        seqid, source, ftype, start, end, score, strand, phase, attr_string = rec

        if ftype in child_types:
            continue

        rebuilt.append(rec)

        if ftype in ("mRNA", "transcript"):
            attrs = parse_attributes(attr_string)
            mrna_id = attrs.get("ID", "")
            if mrna_id and mrna_id in children_by_mrna:
                groups = children_by_mrna[mrna_id]
                for key in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR", "UTR"):
                    for _, cidx in groups[key]:
                        rebuilt.append(new_features[cidx])
                        consumed.add(cidx)

    for idx, rec in enumerate(new_features):
        if idx not in consumed and rec not in rebuilt:
            rebuilt.append(rec)

    new_features = rebuilt
    # ---- end renumbering + grouped reordering ----

    known = set(ordered_seqids)
    final_features = []
    for chrom in ordered_seqids:
        for feat in new_features:
            if feat[0] == chrom:
                final_features.append(feat)

    for feat in new_features:
        if feat[0] not in known:
            final_features.append(feat)

    # Writing progress
    pbar_out = _new_pbar(len(final_features), "Writing GFF3")

    with open(args.output, "w") as out:
        out.write("##gff-version 3\n")
        for feat in final_features:
            out.write("\t".join(feat) + "\n")
            pbar_out.update(1)

    pbar_out.close()

    with open(args.map, "w", newline="") as map_out:
        writer = csv.writer(map_out, delimiter="\t")
        writer.writerow(["Original_ID", "New_ID"])
        for orig, new in sorted(mrna_map.items(), key=lambda x: x[1]):
            writer.writerow([orig, new])

if __name__ == "__main__":
    main()
