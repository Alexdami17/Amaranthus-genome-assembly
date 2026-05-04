#!/usr/bin/env python3
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from typing import Dict, List, Tuple

# ----------------------------
# Plot style helpers
# ----------------------------
def apply_editable_vector_font_settings(fmt: str) -> None:
    """
    Ensure editable text in SVG and reasonable font embedding in PDF.
    - SVG: keep text as text (not paths) so Inkscape can edit it.
    - PDF: use TrueType (fonttype=42) so Illustrator/Inkscape typically edits text better.
    """
    fmt = fmt.lower()
    if fmt == "svg":
        mpl.rcParams["svg.fonttype"] = "none"
    elif fmt == "pdf":
        mpl.rcParams["pdf.fonttype"] = 42
        mpl.rcParams["ps.fonttype"] = 42

def classic_axes(ax):
    """Remove top/right spines; keep left/bottom (classic axes)."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")
    return ax

def save_figure(outfile: str, fmt: str, dpi: int = 200) -> None:
    fmt = fmt.lower()
    if fmt == "png":
        plt.savefig(outfile, dpi=dpi)
    else:
        plt.savefig(outfile)  # vector formats

# ----------------------------
# Parsing helpers
# ----------------------------
def parse_accession(colname: str) -> str:
    # acanthoHap1 -> acantho ; cruentus -> cruentus
    m = re.match(r"^(.*?)(Hap[12])?$", colname)
    return m.group(1) if m else colname

def read_genecount(genecount_tsv: str) -> Tuple[pd.DataFrame, str]:
    df = pd.read_csv(genecount_tsv, sep="\t", dtype=str)
    og_col = df.columns[0]
    return df, og_col

def to_numeric_matrix(df: pd.DataFrame, og_col: str) -> Tuple[np.ndarray, List[str], np.ndarray]:
    # Drop OrthoFinder "Total" column if present; also drop OG id col
    cols = [c for c in df.columns if c not in {og_col, "Total"}]
    mat = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(dtype=float)
    ogs = df[og_col].astype(str).to_numpy()
    return mat, cols, ogs

def collapse_haplotypes(pa: np.ndarray, cols: List[str]) -> Tuple[np.ndarray, List[str]]:
    """
    OR-collapse Hap1/Hap2 columns to accession level.
    Presence(accession) = Presence(Hap1) OR Presence(Hap2) OR ... (if only one exists, that one).
    """
    accs = [parse_accession(c) for c in cols]
    uniq_accs: List[str] = []
    acc_to_idxs: Dict[str, List[int]] = {}

    for j, a in enumerate(accs):
        acc_to_idxs.setdefault(a, []).append(j)
        if a not in uniq_accs:
            uniq_accs.append(a)

    collapsed = np.zeros((pa.shape[0], len(uniq_accs)), dtype=bool)
    for k, a in enumerate(uniq_accs):
        idxs = acc_to_idxs[a]
        collapsed[:, k] = np.any(pa[:, idxs], axis=1)
    return collapsed, uniq_accs

# ----------------------------
# Compartment labeling (5 categories)
# ----------------------------
def label_5cat_by_int(occ_n: np.ndarray,
                      n_genomes: int,
                      soft_core_min: int,
                      shell_min: int,
                      cloud_max: int) -> pd.Categorical:
    """
    5 categories: Core, Soft-core, Shell, Cloud, Private
      Core:      occ_n == n_genomes
      Soft-core: soft_core_min .. n_genomes-1
      Private:   occ_n == 1
      Cloud:     2..cloud_max
      Shell:     shell_min .. (soft_core_min-1) excluding Cloud/Private
    """
    lab = np.full(len(occ_n), "Shell", dtype=object)
    lab[occ_n == n_genomes] = "Core"
    lab[(occ_n >= soft_core_min) & (occ_n < n_genomes)] = "Soft-core"
    lab[occ_n == 1] = "Private"
    lab[(occ_n >= 2) & (occ_n <= cloud_max)] = "Cloud"

    # Ensure very low occupancy doesn't remain Shell accidentally:
    lab[occ_n < shell_min] = np.where(lab[occ_n < shell_min] == "Shell", "Cloud", lab[occ_n < shell_min])

    cats = ["Core", "Soft-core", "Shell", "Cloud", "Private"]
    return pd.Categorical(lab, categories=cats, ordered=True)

# ----------------------------
# TSV summary output
# ----------------------------
def write_compartment_summary(prefix: str, comp_table: pd.DataFrame) -> None:
    """
    Writes:
      - total orthogroups
      - counts per compartment
      - percentages per compartment
    """
    categories = ["Core", "Soft-core", "Shell", "Cloud", "Private"]
    total = int(comp_table.shape[0])
    counts = comp_table["compartment"].value_counts().reindex(categories).fillna(0).astype(int)

    rows = [("total_orthogroups", total)]
    for c in categories:
        rows.append((f"n_{c}", int(counts[c])))
    for c in categories:
        pct = (counts[c] / total * 100.0) if total > 0 else 0.0
        rows.append((f"pct_{c}", float(pct)))

    out = pd.DataFrame(rows, columns=["metric", "value"])
    out.to_csv(f"{prefix}_pangenome_summary.tsv", sep="\t", index=False)

# ----------------------------
# Plotting / summaries
# ----------------------------
def summarize(prefix: str,
              ogs: np.ndarray,
              pa: np.ndarray,
              cols: List[str],
              comp: pd.Categorical,
              fmt: str) -> pd.DataFrame:
    occ_n = pa.sum(axis=1)
    occ_frac = occ_n / pa.shape[1]

    out = pd.DataFrame({
        "Orthogroup": ogs,
        "occ_n": occ_n,
        "occ_frac": occ_frac,
        "compartment": comp
    })
    out.to_csv(f"{prefix}_orthogroup_compartments.tsv", sep="\t", index=False)

    counts = out["compartment"].value_counts().reindex(comp.categories).fillna(0).astype(int)
    counts_df = counts.reset_index()
    counts_df.columns = ["compartment", "n_orthogroups"]
    counts_df.to_csv(f"{prefix}_compartment_counts.tsv", sep="\t", index=False)

    write_compartment_summary(prefix, out)

    # Bar chart (keep default color cycle -> blue)
    fig, ax = plt.subplots()
    ax.bar(counts_df["compartment"], counts_df["n_orthogroups"])
    ax.set_title(f"{prefix}: pangenome compartments (OrthoFinder orthogroups)")
    ax.set_ylabel("Number of orthogroups")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    classic_axes(ax)
    fig.tight_layout()
    save_figure(f"{prefix}_compartments_bar.{fmt}", fmt=fmt, dpi=200)
    plt.close(fig)

    return out

def _pick_center_and_band(pan: np.ndarray, core: np.ndarray, band: str):
    """
    band:
      - '95'  => mean line + 2.5–97.5% envelope
      - 'iqr' => median line + 25–75% envelope
    Returns: center_pan, center_core, lo_q, hi_q
    """
    band = band.lower()
    if band == "iqr":
        center_pan = np.nanmedian(pan, axis=0)
        center_core = np.nanmedian(core, axis=0)
        lo_q, hi_q = 0.25, 0.75
    elif band == "95":
        center_pan = np.nanmean(pan, axis=0)
        center_core = np.nanmean(core, axis=0)
        lo_q, hi_q = 0.025, 0.975
    else:
        raise ValueError(f"Invalid --band value: {band}. Use 'iqr' or '95'.")
    return center_pan, center_core, lo_q, hi_q

def accumulation(prefix: str,
                 pa: np.ndarray,
                 fmt: str,
                 band: str,
                 n_perm: int = 500,
                 seed: int = 123) -> None:
    """
    Produces ONE plot:
      1) Absolute (orthogroup counts):   {prefix}_accumulation_curves.{fmt}

    Ribbon is computed across permutations (order effects).
    Center line depends on --band:
      - band=95  => mean line; ribbon = 2.5–97.5%
      - band=iqr => median line; ribbon = 25–75%
    """
    rng = np.random.default_rng(seed)
    n_og, n_gen = pa.shape

    present_lists = [np.flatnonzero(pa[:, j]) for j in range(n_gen)]
    pan = np.empty((n_perm, n_gen), dtype=float)
    core = np.empty((n_perm, n_gen), dtype=float)

    for p in range(n_perm):
        order = rng.permutation(n_gen)
        union_set = set()
        inter_set = None
        for k, j in enumerate(order):
            idx = present_lists[j].tolist()
            union_set.update(idx)
            if inter_set is None:
                inter_set = set(idx)
            else:
                inter_set.intersection_update(idx)

            pan[p, k] = len(union_set)
            core[p, k] = len(inter_set)

    def q(arr, prob): return np.nanquantile(arr, prob, axis=0)
    x = np.arange(1, n_gen + 1)

    pan_center, core_center, lo_q, hi_q = _pick_center_and_band(pan, core, band)
    pan_lo, pan_hi = q(pan, lo_q), q(pan, hi_q)
    core_lo, core_hi = q(core, lo_q), q(core, hi_q)

    # Colors: Pangenome=green (C2), Core=orange (C1)
    fig, ax = plt.subplots()
    ax.plot(x, pan_center, linewidth=2, color="C2", label="Pangenome")
    ax.plot(x, core_center, linewidth=2, color="C1", label="Core genome")
    ax.fill_between(x, pan_lo, pan_hi, color="C2", alpha=0.15)
    ax.fill_between(x, core_lo, core_hi, color="C1", alpha=0.15)

    ax.set_title(f"{prefix}: accumulation curves (n={n_gen}, {n_perm} permutations)")
    ax.set_xlabel("Number of genomes added")
    ax.set_ylabel("Number of orthogroups")
    classic_axes(ax)
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(f"{prefix}_accumulation_curves.{fmt}", fmt=fmt, dpi=200)
    plt.close(fig)

# ----------------------------
# Gene list export from Orthogroups.tsv
# ----------------------------
def export_gene_lists(orthogroups_tsv: str, comp_table: pd.DataFrame, prefix: str) -> None:
    """
    Exports:
      - {prefix}_genes_{Compartment}.txt (unique gene IDs)
      - {prefix}_gene_to_orthogroup_compartment.tsv
    """
    og_to_comp = dict(zip(comp_table["Orthogroup"].astype(str), comp_table["compartment"].astype(str)))

    og_df = pd.read_csv(orthogroups_tsv, sep="\t", dtype=str).fillna("")
    og_col = og_df.columns[0]
    genome_cols = [c for c in og_df.columns if c != og_col and c != "Total"]

    comp_to_genes: Dict[str, set] = {c: set() for c in ["Core", "Soft-core", "Shell", "Cloud", "Private"]}
    gene_map_rows = []

    for _, row in og_df.iterrows():
        og = str(row[og_col])
        c = og_to_comp.get(og, None)
        if c is None:
            continue

        for gc in genome_cols:
            cell = row[gc].strip()
            if not cell:
                continue
            parts = [p.strip() for p in cell.split(",") if p.strip()]
            for g in parts:
                comp_to_genes[c].add(g)
                gene_map_rows.append((g, og, c, gc))

    for c, genes in comp_to_genes.items():
        Path(f"{prefix}_genes_{c}.txt").write_text("\n".join(sorted(genes)) + "\n")

    gene_map = pd.DataFrame(gene_map_rows, columns=["gene_id", "Orthogroup", "compartment", "source_column"])
    gene_map.drop_duplicates().to_csv(f"{prefix}_gene_to_orthogroup_compartment.tsv", sep="\t", index=False)

# ----------------------------
# CLI
# ----------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Gene-centric pangenome from OrthoFinder outputs (Option A assemblies; Option B haplotype-collapsed accessions)."
    )
    p.add_argument("genecount_tsv", help="OrthoFinder Orthogroups.GeneCount.tsv")
    p.add_argument("orthogroups_tsv", nargs="?", default=None,
                   help="Optional OrthoFinder Orthogroups.tsv (enables per-compartment gene list export)")

    fmt_group = p.add_mutually_exclusive_group()
    fmt_group.add_argument("--png", action="store_true", help="Write PNG figures (default)")
    fmt_group.add_argument("--pdf", action="store_true", help="Write PDF figures")
    fmt_group.add_argument("--svg", action="store_true", help="Write SVG figures (editable text for Inkscape)")

    p.add_argument("--band", choices=["iqr", "95"], default="95",
                   help="Permutation ribbon: '95' (2.5–97.5%%) or 'iqr' (25–75%%). "
                        "Note: 'iqr' uses median center line automatically. Default: 95")

    p.add_argument("--n-perm", type=int, default=500, help="Permutations for accumulation curves (default: 500)")
    p.add_argument("--seed", type=int, default=123, help="Random seed (default: 123)")

    p.add_argument("--cloud-max", type=int, default=2, help="Cloud occupancy = 2..cloud_max (default: 2)")
    p.add_argument("--shell-min", type=int, default=3, help="Shell minimum occupancy (default: 3)")
    p.add_argument("--softcore-missing", type=int, default=1,
                   help="Soft-core allows this many missing genomes: soft_core_min = n_genomes - softcore_missing (default: 1)")

    return p.parse_args()

# ----------------------------
# Main
# ----------------------------
def main():
    args = parse_args()

    fmt = "png"
    if args.pdf:
        fmt = "pdf"
    elif args.svg:
        fmt = "svg"
    elif args.png:
        fmt = "png"

    apply_editable_vector_font_settings(fmt)

    df, og_col = read_genecount(args.genecount_tsv)
    mat, cols, ogs = to_numeric_matrix(df, og_col)
    pa_A = mat > 0

    # ---- Option A (assemblies)
    nA = pa_A.shape[1]
    soft_core_min_A = max(nA - args.softcore_missing, 2)

    compA = label_5cat_by_int(
        occ_n=pa_A.sum(axis=1),
        n_genomes=nA,
        soft_core_min=soft_core_min_A,
        shell_min=args.shell_min,
        cloud_max=args.cloud_max
    )
    compA_table = summarize("A", ogs, pa_A, cols, compA, fmt=fmt)
    accumulation("A", pa_A, fmt=fmt, band=args.band, n_perm=args.n_perm, seed=args.seed)

    # ---- Option B (collapsed accessions)
    pa_B, acc_cols = collapse_haplotypes(pa_A, cols)
    nB = pa_B.shape[1]
    soft_core_min_B = max(nB - args.softcore_missing, 2)

    compB = label_5cat_by_int(
        occ_n=pa_B.sum(axis=1),
        n_genomes=nB,
        soft_core_min=soft_core_min_B,
        shell_min=args.shell_min,
        cloud_max=args.cloud_max
    )
    compB_table = summarize("B", ogs, pa_B, acc_cols, compB, fmt=fmt)
    accumulation("B", pa_B, fmt=fmt, band=args.band, n_perm=args.n_perm, seed=args.seed)

    # ---- Gene list export
    if args.orthogroups_tsv:
        export_gene_lists(args.orthogroups_tsv, compA_table, "A")
        export_gene_lists(args.orthogroups_tsv, compB_table, "B")

    print("Done.")
    print(f"Figures written as .{fmt} (A_*, B_*).")
    print(f"Ribbon band: {args.band}  (iqr => median+25–75%; 95 => mean+2.5–97.5%)")
    print("Summaries: A_pangenome_summary.tsv, B_pangenome_summary.tsv")

if __name__ == "__main__":
    main()
