# Amaranthus genome assembly

A collection of standalone scripts used in the assembly, annotation, and
comparative analysis of *Amaranthus* genomes. Each script is self-contained and
exposes a command-line interface (run with `-h`/`--help` for full options).

## Repeat analysis

- **`create_repeat_landscape.py`** - Builds a manuscript-ready repeat landscape
  plot from a RepeatMasker `*.divsum` file (output of
  `calcDivergenceFromAlign.pl`). Produces spaced stacked bars of repeat
  divergence with hard-coded, colorblind-safe (Okabe-Ito) colors grouped by
  repeat class. Outputs SVG/PDF/PNG with editable text; supports class selection
  via `--features`/`--classes-file` and genome-size normalisation.
- **`filter_repeatmasker_align_file.py`** - Filters a RepeatMasker `.align` file
  by query sequence (chromosome/contig). Streams the file record-by-record,
  preserving original formatting, and supports exact name lists (`--chr`), regex
  matching (`--chr-regex`), inverted selection (`--invert`), and optional
  splitting into one file per chromosome.

## Gene annotation / GFF3 processing

- **`generate_gff3_from_entap_and_pasa.py`** - Generates a GFF3 annotation by
  combining EnTAP functional annotation with PASA structural annotation.
- **`merge_tRNAScan-SE_gff3_with_entap_derived_gff3.py`** - Inserts a
  coordinate-sorted tRNAscan-SE GFF3 (file B) into an existing GFF3 (file A)
  while preserving A's line order. Normalises file B's features (`pseudogene` →
  `tRNA`, ID/Parent prefixing to match A, `exonN` → `exon.N`, tidy trailing
  separators), models tRNAs as transcript-level features under a gene container,
  and writes a single `##gff-version 3` header.
- **`find_and_remove_mRNAs_without_CDS_from_gff3.py`** - Scans a GFF3 (or `.gz`)
  for mRNA/transcript features that lack any CDS, logs a summary, and optionally
  rewrites the file with those transcripts removed along with newly orphaned
  genes and any child features that referenced them.

## Ka/Ks (selection) analysis

- **`orthofinder_or_mmseqs_to_kaks.py`** - End-to-end pipeline to compute Ka/Ks
  for 1:1 pairs from either an OrthoFinder TSV or an MMseqs reciprocal-best-hit
  TSV (with optional same-chromosome filtering). Extracts protein/CDS sequences,
  aligns proteins with MAFFT, builds codon alignments with PAL2NAL, assembles AXT
  files, and runs KaKs_Calculator 3.0.
- **`kaks_calculator_output_parser.py`** - Parses KaKs_Calculator output into a
  tidy TSV (`gene1`, `gene2`, `Ka`, `Ks`, `Ka/Ks`, `Pvalue`) and can optionally
  join it against a GFF3 to sort results by gene start coordinate.

## Pangenome

- **`orthofinder_to_pangenome_5groups.py`** - Converts an OrthoFinder gene-count
  matrix into a pangenome classification across five occupancy groups and
  produces editable vector plots.

## Strata / changepoint analysis (R)

- **`fit_mcp_strata_analysis.R`** - Driver script for analysing Ks along a
  chromosome: QC with LOESS/GAM smoothing and Rosner outlier detection, then
  Bayesian changepoint ("strata") modelling. Sources the helper functions below
  and saves QC/GAM plots.
- **`fit_mcp_strata_helper_scripts.R`** - Helper function library used by the
  driver: `analyze_ks()` (Rosner outliers + LOESS/GAM), `fit_mcp_strata()`
  (Rosner outliers + `mcp` changepoint model + LOO comparison), and
  `plot_mcp_cp_posterior()` (changepoint posterior density plotting).
