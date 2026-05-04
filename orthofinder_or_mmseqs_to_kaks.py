#!/usr/bin/env python3
"""
orthofinder_kaks_pipeline.py

Pipeline to compute Ka/Ks for 1:1 pairs from either:

1. OrthoFinder TSV (orthogroup + two gene IDs), *or*
2. MMseqs RBH TSV (query + target IDs), with optional same-chromosome filtering.

Workflow:

1. Read input TSV:
   - OrthoFinder mode: orthogroup + two gene IDs, keep only 1:1 groups.
   - MMseqs mode: query + target IDs, keep only strict 1-to-1 RBH pairs,
     optionally requiring same chromosome (Chr01 == Chr01, etc.).
2. Optionally restrict to one or more chromosomes.
3. Extract protein/CDS sequences for each pair.
4. Align proteins with MAFFT.
5. Generate codon alignments with PAL2NAL.
6. Build quasi-AXT file (gene1-gene2 + two codon-aligned sequences).
7. Run KaKs_Calculator3.0 on all pairs.

Examples
--------

OrthoFinder mode (as before):

    python orthofinder_kaks_pipeline.py \
        --input-mode orthofinder \
        -i orthofinder_pairs.tsv \
        -p speciesA_prot.faa speciesB_prot.faa \
        -c speciesA_cds.fna  speciesB_cds.fna \
        -o my_kaks_out \
        --KaKs /path/to/KaKs \
        --KaKs_opts -m YN -c 1

MMseqs RBH mode (RBH TSV from `mmseqs easy-rbh` or `convertalis`):

    mmseqs easy-rbh X.faa Y.faa XY_rbh tmp \
        --min-seq-id 0.25 \
        -e 1e-5 \
        -c 0.4

    mmseqs convertalis X.faa Y.faa XY_rbh XY_rbh.tsv \
        --format-output "query,target,fident,alnLen,qlen,tlen,qstart,qend,tstart,tend,evalue,bits"

    python orthofinder_kaks_pipeline.py \
        --input-mode mmseqs \
        -i XY_rbh.tsv \
        -p X_prot.faa Y_prot.faa \
        -c X_cds.fna  Y_cds.fna \
        -o XY_kaks_out \
        --chr-regex "(Chr\\d+)" \
        --KaKs /path/to/KaKs \
        --KaKs_opts -m YN -c 1

Chromosome-related options
--------------------------

- `--chr-regex` (default: "(Chr\\d+)"):
    Regex with one capturing group for chromosome name, used in MMseqs mode
    (and optionally when filtering FASTAs).
    Example: "AmaPaChr01g010480" -> "Chr01".

- `-C / --chromosome`:
    One or more chromosome names to keep (e.g., -C Chr01 Chr05 Chr16).
    In MMseqs mode:
        *only* pairs where chr(query) == chr(target) and that chromosome
        is in this list are kept.

In OrthoFinder mode, behaviour is unchanged except that -C is interpreted
as "both IDs must contain the SAME substring from the list."
"""

import argparse
import shutil
import subprocess
import logging
import json
import re
import sys
from pathlib import Path
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    import yaml
except ImportError:
    yaml = None


def setup_logging(log_dir):
    """Configure logging to file and console, writing pipeline.log in log_dir."""
    log_file = Path(log_dir) / "pipeline.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger("ortholog_pipeline")


def normalize_id(text, id_regex=None):
    """Normalize an ID string using regex or first token fallback."""
    text = text.strip()
    if id_regex:
        m = re.search(id_regex, text)
        if m:
            return m.group(0)
    return text.split()[0]


def load_fasta_sequences(fasta_files, id_regex=None, chromosomes=None, chr_regex=None):
    """
    Load sequences from multiple FASTA files into a dict keyed by normalized IDs.

    If 'chr_regex' is provided, it is used to parse the chromosome (first capturing group).
    If 'chromosomes' is provided (list of chromosome names), only keep sequences
    whose parsed chromosome name is in this list.

    If chr_regex is None but chromosomes is provided, we fall back on substring matching.
    """
    seq_dict = {}
    chr_pattern = re.compile(chr_regex) if chr_regex else None

    for fasta in fasta_files or []:
        for record in SeqIO.parse(fasta, "fasta"):
            clean_id = normalize_id(record.description, id_regex)

            if chromosomes:
                keep = False
                if chr_pattern:
                    m = chr_pattern.search(clean_id)
                    if m and m.group(1) in chromosomes:
                        keep = True
                else:
                    # fallback: any chromosome substring contained in ID
                    if any(ch in clean_id for ch in chromosomes):
                        keep = True
                if not keep:
                    continue

            seq_dict[clean_id] = str(record.seq)

    return seq_dict


def run_mafft(mafft_path, infile, outfile, extra_opts, log_file):
    """Run MAFFT alignment, logging verbose output to log_file."""
    cmd = [mafft_path] + (extra_opts or []) + [infile]
    with open(outfile, "w") as out, open(log_file, "w") as lf:
        subprocess.run(cmd, stdout=out, stderr=lf, check=True)


def run_pal2nal(pal2nal_path, prot_aln, cds_file, outfile, extra_opts, log_file):
    """
    Run PAL2NAL codon alignment, logging verbose output to log_file.

    We use -output fasta so we can directly convert to quasi-AXT in Python.
    """
    cmd = [pal2nal_path, prot_aln, cds_file, "-output", "fasta"]
    if extra_opts:
        cmd += extra_opts
    with open(outfile, "w") as out, open(log_file, "w") as lf:
        subprocess.run(cmd, stdout=out, stderr=lf, check=True)


def run_KaKs(KaKs_path, infile, outfile, extra_opts, log_file):
    """Run KaKs_Calculator3.0 on the quasi-AXT file."""
    cmd = [KaKs_path, "-i", infile, "-o", outfile] + (extra_opts or [])
    with open(log_file, "w") as lf:
        subprocess.run(cmd, stderr=lf, check=True)


def summarize_logs(log_dir, summary_file, logger):
    """Scan all .log files and collect those containing errors/warnings."""
    issues = []
    for log in Path(log_dir).glob("*.log"):
        try:
            with open(log) as f:
                content = f.read()
                if "ERROR" in content or "fail" in content.lower() or "warning" in content.lower():
                    issues.append(log.name)
        except Exception as e:
            logger.warning(f"Could not read {log}: {e}")
    with open(summary_file, "w") as out:
        for logname in issues:
            out.write(f"{logname}\n")
    logger.info(f"Log summary written: {summary_file} ({len(issues)} problematic logs)")
    return issues


def process_orthogroup(og_id, seq_ids, prot_dict, cds_dict,
                       prot_cds_dir, prot_aln_dir, cds_aln_dir, log_dir,
                       mafft_path, mafft_opts, pal2nal_path,
                       pal2nal_opts, dry_run=False):
    """
    Process one orthogroup/pair: write FASTA, run MAFFT, run PAL2NAL.

    NOTE: Logger is obtained inside to avoid pickling issues with ProcessPoolExecutor.
    """
    logger = logging.getLogger("ortholog_pipeline")
    try:
        # Write protein FASTA into prot_cds_seq
        og_prot_fasta = prot_cds_dir / f"{og_id}.prot.fa"
        with open(og_prot_fasta, "w") as out:
            for seq_id in seq_ids:
                if seq_id in prot_dict:
                    out.write(f">{seq_id}\n{prot_dict[seq_id]}\n")
                else:
                    logger.warning(f"{og_id}: Protein {seq_id} not found")

        logger.info(f"{og_id}: Protein FASTA written")

        if dry_run:
            logger.info(f"{og_id}: Dry-run mode, skipping MAFFT/PAL2NAL")
            return (og_id, True)

        # Run MAFFT on proteins
        prot_aln_file = prot_aln_dir / f"{og_id}.prot_aln.fa"
        if mafft_path:
            mafft_log = log_dir / f"{og_id}.mafft.log"
            logger.info(f"{og_id}: Running MAFFT, see {mafft_log} for details...")
            run_mafft(mafft_path, str(og_prot_fasta), str(prot_aln_file),
                      mafft_opts, str(mafft_log))
            logger.info(f"{og_id}: MAFFT alignment complete")

            # Run PAL2NAL if CDS provided
            if cds_dict and pal2nal_path:
                # Write CDS FASTA
                og_cds_fasta = prot_cds_dir / f"{og_id}.cds.fa"
                with open(og_cds_fasta, "w") as out:
                    for seq_id in seq_ids:
                        if seq_id in cds_dict:
                            out.write(f">{seq_id}\n{cds_dict[seq_id]}\n")
                        else:
                            logger.warning(f"{og_id}: CDS {seq_id} not found")

                cds_aln_file = cds_aln_dir / f"{og_id}.cds_aln.fasta"
                pal2nal_log = log_dir / f"{og_id}.pal2nal.log"
                logger.info(f"{og_id}: Running PAL2NAL, see {pal2nal_log} for details...")
                run_pal2nal(pal2nal_path, str(prot_aln_file), str(og_cds_fasta),
                            str(cds_aln_file), pal2nal_opts, str(pal2nal_log))
                logger.info(f"{og_id}: PAL2NAL codon alignment complete")

        return (og_id, True)
    except Exception as e:
        logger.error(f"{og_id}: FAILED with error {e}")
        return (og_id, False)


def fasta_to_axt_block(fasta_path: Path, handle, logger=None):
    """
    Convert a 2-sequence codon FASTA alignment into one quasi-AXT block.

    Header format:
        gene1-gene2
    """
    if logger is None:
        logger = logging.getLogger("ortholog_pipeline")

    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(records) != 2:
        logger.warning(f"{fasta_path.name}: expected 2 sequences, found {len(records)}; skipping")
        return False

    r1, r2 = records

    # Force the AXT header to be geneID1-geneID2
    header = f"{r1.id}-{r2.id}"

    handle.write(f"{header}\n")
    handle.write(str(r1.seq) + "\n")
    handle.write(str(r2.seq) + "\n\n")

    return True


# ------------------------------
#  Pair parsing helpers
# ------------------------------

def parse_pairs_orthofinder(tsv_path, id_regex, chromosomes):
    """
    Parse OrthoFinder TSV: orthogroup + two IDs, keep only 1:1 single-copy groups.

    Chromosome behaviour is unchanged from your original script:
    if 'chromosomes' is provided (list of substrings), we keep
    only pairs where BOTH IDs contain the SAME substring from this list.
    """
    og_counts, row_map = {}, {}
    with open(tsv_path) as infile:
        header = next(infile, None)  # skip header if present
        for line in infile:
            # This assumes multi-copy orthogroups have commas in the IDs columns
            if "," in line:
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            og_id = parts[0].strip()
            seq1 = normalize_id(parts[1], id_regex)
            seq2 = normalize_id(parts[2], id_regex)
            og_counts[og_id] = og_counts.get(og_id, 0) + 1
            row_map[og_id] = (seq1, seq2)

    tasks = []
    for og_id, count in og_counts.items():
        if count != 1:
            continue
        seq1, seq2 = row_map[og_id]

        if chromosomes:
            # require that BOTH IDs share at least one chromosome substring
            if not any((ch in seq1) and (ch in seq2) for ch in chromosomes):
                continue

        tasks.append((og_id, [seq1, seq2]))
    return tasks


def parse_pairs_mmseqs(tsv_path, id_regex, chromosomes, chr_regex):
    """
    Parse MMseqs TSV (e.g. from easy-rbh / convertalis):

      col0: query ID
      col1: target ID
      (extra columns ignored here)

    Logic:
      1) Optionally extract chromosome from IDs via chr_regex (first group).
         - If chr_regex is given:
             - require chr(query) and chr(target) both defined and equal.
             - if 'chromosomes' list is given, require chr in this list.
         - If chr_regex is None:
             - if 'chromosomes' list is given, require that
               some substring in chromosomes appears in both IDs.
      2) Count how many times each ID appears as query/target.
      3) Keep only strict 1-to-1 pairs, i.e. genes that occur exactly once
         across all retained candidate pairs.
    """
    chr_pattern = re.compile(chr_regex) if chr_regex else None
    candidate_rows = []
    q_counts = {}
    t_counts = {}

    with open(tsv_path) as infile:
        for line in infile:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                continue

            qid = normalize_id(parts[0], id_regex)
            tid = normalize_id(parts[1], id_regex)

            # Same-chromosome + chromosome subset logic
            if chr_pattern:
                mq = chr_pattern.search(qid)
                mt = chr_pattern.search(tid)
                if not mq or not mt:
                    continue
                chr_q = mq.group(1)
                chr_t = mt.group(1)
                if chr_q != chr_t:
                    continue
                if chromosomes and chr_q not in chromosomes:
                    continue
            else:
                if chromosomes:
                    # fallback: require any of the substrings to appear in BOTH IDs
                    if not any((ch in qid) and (ch in tid) for ch in chromosomes):
                        continue

            candidate_rows.append((qid, tid))
            q_counts[qid] = q_counts.get(qid, 0) + 1
            t_counts[tid] = t_counts.get(tid, 0) + 1

    # Enforce strict 1-to-1
    tasks = []
    pair_idx = 0
    for qid, tid in candidate_rows:
        if q_counts[qid] == 1 and t_counts[tid] == 1:
            og_id = f"pair_{pair_idx}"
            pair_idx += 1
            tasks.append((og_id, [qid, tid]))

    return tasks


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Wrapper: parse OrthoFinder or MMseqs pairs, align with MAFFT+PAL2NAL, "
            "compute Ka/Ks with KaKs_Calculator."
        )
    )
    parser.add_argument("--config", help="Optional YAML/JSON config file")
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate inputs and setup directories without running MAFFT/PAL2NAL/KaKs")
    parser.add_argument("-i", "--tsv", help="Path to TSV file (OrthoFinder or MMseqs)")
    parser.add_argument("-p", "--prot", nargs="+", help="Protein FASTA files")
    parser.add_argument("-c", "--cds", nargs="+", help="CDS FASTA files")
    parser.add_argument("-o", "--outdir", default="ortholog_sequences", help="Output directory")
    parser.add_argument("-m", "--mafft", default=shutil.which("mafft"), help="Path to MAFFT")
    parser.add_argument("--mafft_opts", nargs="*", default=["--genafpair", "--maxiterate", "1000"])
    parser.add_argument("--pal2nal", default=shutil.which("pal2nal.pl"), help="Path to PAL2NAL")
    parser.add_argument("--pal2nal_opts", nargs="*", default=[], help="Extra options for pal2nal.pl")
    parser.add_argument("--KaKs", default=shutil.which("KaKs"), help="Path to KaKs_Calculator3.0")
    parser.add_argument("--KaKs_opts", nargs="*", default=[],
                        help="Extra options for KaKs_Calculator3.0 (e.g., -c, -m, -d)")
    parser.add_argument("--id-regex", help="Regex pattern to extract IDs from FASTA/TSV")

    parser.add_argument(
        "-C", "--chromosome",
        nargs="+",
        help=(
            "One or more chromosome names to keep (e.g., -C Chr01 Chr05 Chr16). "
            "In OrthoFinder mode: a pair is kept only if BOTH IDs contain the SAME substring. "
            "In MMseqs mode: a pair is kept only if chr(query) == chr(target) and that "
            "chromosome is in this list (if --chr-regex is provided)."
        )
    )

    parser.add_argument(
        "--chr-regex",
        default=r"(Chr\d+)",
        help="Regex with one capturing group for chromosome name (used mainly in MMseqs mode)."
    )

    parser.add_argument(
        "--input-mode",
        choices=["orthofinder", "mmseqs"],
        default="orthofinder",
        help="How to interpret the TSV input: 'orthofinder' (orthogroup + 2 IDs) or 'mmseqs' (query + target)."
    )

    parser.add_argument("-t", "--threads", type=int, default=4)

    args = parser.parse_args()

    # Load config file if provided (config keys override defaults if CLI didn't set them)
    if args.config:
        cfg = Path(args.config)
        with open(cfg) as f:
            if cfg.suffix in [".yaml", ".yml"] and yaml:
                config_data = yaml.safe_load(f)
            else:
                config_data = json.load(f)
        for k, v in config_data.items():
            if hasattr(args, k) and getattr(args, k) in [None, parser.get_default(k)]:
                setattr(args, k, v)

    # Normalize chromosome argument
    if args.chromosome is None:
        chromosomes = None
    elif isinstance(args.chromosome, str):
        chromosomes = [args.chromosome]
    else:
        chromosomes = list(args.chromosome)

    # Create directories
    outdir = Path(args.outdir)
    prot_cds_dir = outdir / "prot_cds_seq"
    prot_aln_dir = outdir / "prot_aln_seq"
    cds_aln_dir = outdir / "cds_aln_seq"
    log_dir = outdir / "log"

    for d in [outdir, log_dir, prot_cds_dir, prot_aln_dir]:
        d.mkdir(parents=True, exist_ok=True)
    if args.cds:
        cds_aln_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    logger = setup_logging(log_dir)
    logger.info("Pipeline started")

    # === Sanity checks on required arguments ===
    if not args.tsv:
        logger.error("No TSV file provided. Use -i / --tsv.")
        sys.exit(1)

    if not args.prot:
        logger.error("No protein FASTA provided. Use -p / --prot.")
        sys.exit(1)

    if not args.cds:
        logger.error("No CDS FASTA provided. Use -c / --cds.")
        sys.exit(1)

    # === Check external tools ===
    if not args.mafft:
        logger.error("MAFFT not found in PATH and no --mafft path provided.")
        sys.exit(1)

    if not args.pal2nal:
        logger.error("pal2nal.pl not found in PATH and no --pal2nal path provided.")
        sys.exit(1)

    if not args.KaKs:
        logger.error("KaKs_Calculator3.0 not found in PATH and no --KaKs path provided.")
        sys.exit(1)

    # Load sequences (optionally filtered by chromosome names)
    prot_dict = load_fasta_sequences(args.prot, id_regex=args.id_regex,
                                     chromosomes=chromosomes, chr_regex=args.chr_regex)
    cds_dict = load_fasta_sequences(args.cds, id_regex=args.id_regex,
                                    chromosomes=chromosomes, chr_regex=args.chr_regex) if args.cds else {}

    # Parse pairs depending on input mode
    if args.input_mode == "orthofinder":
        tasks = parse_pairs_orthofinder(args.tsv, args.id_regex, chromosomes)
    else:
        tasks = parse_pairs_mmseqs(args.tsv, args.id_regex, chromosomes, args.chr_regex)

    # Track missing IDs
    missing_ids = []
    for og_id, seq_ids in tasks:
        for sid in seq_ids:
            if sid not in prot_dict and sid not in cds_dict:
                missing_ids.append((og_id, sid))

    if missing_ids:
        mis_file = log_dir / "missing_ids.txt"
        with open(mis_file, "w") as f:
            for og, sid in missing_ids:
                f.write(f"{og}\t{sid}\n")
        logger.warning(f"Missing IDs written to {mis_file}")

    # Dry-run shortcut
    if args.dry_run:
        logger.info("Dry-run mode: skipping MAFFT/PAL2NAL/KaKs stages")
        with open(log_dir / "summary.tsv", "w") as f:
            pass
        with open(log_dir / "summary.json", "w") as f:
            json.dump([], f)
        summarize_logs(log_dir, outdir / "log_summary.txt", logger)
        logger.info("Pipeline finished (dry-run).")
        return

    # Run MAFFT/PAL2NAL in parallel
    summary = []
    if tasks:
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            futures = [executor.submit(process_orthogroup, og_id, seq_ids,
                                       prot_dict, cds_dict,
                                       prot_cds_dir, prot_aln_dir, cds_aln_dir, log_dir,
                                       args.mafft, args.mafft_opts,
                                       args.pal2nal, args.pal2nal_opts,
                                       args.dry_run)
                       for og_id, seq_ids in tasks]
            for fut in as_completed(futures):
                og_id, ok = fut.result()
                summary.append({
                    "orthogroup": og_id,
                    "status": "success" if ok else "fail"
                })

    # Structured summary logs into log dir
    with open(log_dir / "summary.tsv", "w") as f:
        for row in summary:
            f.write(f"{row['orthogroup']}\t{row['status']}\n")
    with open(log_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    logger.info("Alignment stages complete. Starting AXT/KaKs post-processing...")

    # If no CDS or PAL2NAL, we can't do Ka/Ks (should already be caught, but keep guard)
    if not args.cds or not cds_dict:
        logger.warning("No CDS sequences loaded; skipping Ka/Ks stage.")
        summarize_logs(log_dir, outdir / "log_summary.txt", logger)
        logger.info(
            f"Pipeline finished. Successes: {sum(1 for r in summary if r['status']=='success')}, "
            f"Failures: {sum(1 for r in summary if r['status']=='fail')}"
        )
        return

    # Build master quasi-AXT from all PAL2NAL FASTA alignments
    all_axt = outdir / "all_pairs.axt"
    num_blocks = 0
    with open(all_axt, "w") as master:
        for aln in sorted(cds_aln_dir.glob("*.cds_aln.fasta")):
            ok = fasta_to_axt_block(aln, master, logger=logger)
            if ok:
                num_blocks += 1

    if num_blocks == 0:
        logger.warning("No valid codon FASTA alignments found; skipping Ka/Ks stage.")
    else:
        logger.info(f"Quasi-AXT file written: {all_axt} ({num_blocks} blocks)")

        # KaKs calculation on all_pairs.axt
        if args.KaKs and all_axt.exists():
            KaKs_out = outdir / "results.KaKs"
            KaKs_log = log_dir / "KaKs.log"
            try:
                run_KaKs(args.KaKs, str(all_axt), str(KaKs_out), args.KaKs_opts, str(KaKs_log))
                logger.info(f"Ka/Ks analysis complete: {KaKs_out}")
            except subprocess.CalledProcessError as e:
                logger.error(f"KaKs_Calculator failed: {e}")
        else:
            if not args.KaKs:
                logger.warning("KaKs_Calculator path not provided; skipping Ka/Ks stage")
            elif not all_axt.exists():
                logger.warning("Master AXT file not found; skipping Ka/Ks stage")

    logger.info("Summary logs written (log/summary.tsv, log/summary.json)")
    logger.info(
        f"Pipeline finished. Successes: {sum(1 for r in summary if r['status']=='success')}, "
        f"Failures: {sum(1 for r in summary if r['status']=='fail')}"
    )

    # Summarize all per-orthogroup logs
    summarize_logs(log_dir, outdir / "log_summary.txt", logger)


if __name__ == "__main__":
    main()
