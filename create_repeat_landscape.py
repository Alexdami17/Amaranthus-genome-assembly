#!/usr/bin/env python3
"""
create_repeat_landscape_v3_4.py

Landscape-only plot from RepeatMasker *.divsum (SVG/PDF/PNG), manuscript-ready.

What this version does:
  - NO palette argument; ALL colors are hard-coded for reproducibility.
  - Colors are CVD-friendly (built from Okabe–Ito colorblind-safe palette + deterministic variants),
    while preserving semantic grouping:
      * DNA/*      -> warm (orange/vermillion/gold) with strong within-family separation
      * LTR/*      -> green/teal (Copia vs Gypsy robust under CVD)
      * LINE/*     -> blues with lightness separation
      * SINE/*     -> purple/pink with lightness separation
      * RC/Helitron-> brown
      * Other      -> light neutral
      * Unknown    -> gray (ONLY Unknown is gray)
  - Select classes via --features and/or --classes-file (TXT)
  - Spaced stacked bars (less histogram-like)
  - Proper axis spines
  - Legend does NOT squeeze plot (margin or panel modes)
  - Legend sorted alphabetically (without changing stacking order)
  - Matplotlib compatibility kept (not used for colors, but harmless)

Input: *.divsum from calcDivergenceFromAlign.pl

Examples:
  ./create_repeat_landscape_v3_4.py --div species.divsum --genome-size 812345678 \
    --out landscape.svg --ylim 0 3.5

  ./create_repeat_landscape_v3_4.py --div species.divsum --twoBit genome.2bit \
    --classes-file classes.txt --out focused.svg --legend-mode panel
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["svg.fonttype"] = "none"   # keep text as <text>, not paths
mpl.rcParams["pdf.fonttype"] = 42       # embed TrueType fonts (Type 42)
mpl.rcParams["ps.fonttype"]  = 42

# Optional but recommended: stable, widely available font
mpl.rcParams["font.family"] = "DejaVu Sans"

from matplotlib.gridspec import GridSpec

# Matplotlib colormap compatibility (>=3.7 vs older) - kept for completeness
# (Not used for colors, since we hardcode all colors.)
try:
    import matplotlib.colormaps as cmaps  # Matplotlib >= 3.7

    def get_cmap(name: str):
        return cmaps.get_cmap(name)

except ImportError:

    def get_cmap(name: str):
        return plt.get_cmap(name)


# -----------------------------
# Canonical label ordering used for stacking (from Perl script)
# -----------------------------
GRAPH_LABELS_ORDER: List[str] = [
    "Unknown",
    "Other",
    "DNA/Academ",
    "DNA/CMC",
    "DNA/Crypton",
    "DNA/Ginger",
    "DNA/Harbinger",
    "DNA/hAT",
    "DNA/Kolobok",
    "DNA/Maverick",
    "DNA",
    "DNA/Merlin",
    "DNA/MULE",
    "DNA/P",
    "DNA/PiggyBac",
    "DNA/Sola",
    "DNA/TcMar",
    "DNA/Transib",
    "DNA/Zator",
    "DNA/Dada",
    "RC/Helitron",
    "LTR/DIRS",
    "LTR/Ngaro",
    "LTR/Pao",
    "LTR/Copia",
    "LTR/Gypsy",
    "LTR/ERVL",
    "LTR",
    "LTR/ERV1",
    "LTR/ERV",
    "LTR/ERVK",
    "LINE/L1",
    "LINE",
    "LINE/RTE",
    "LINE/CR1",
    "LINE/Rex-Babar",
    "LINE/L2",
    "LINE/Proto2",
    "LINE/LOA",
    "LINE/R1",
    "LINE/Jockey-I",
    "LINE/Dong-R4",
    "LINE/R2",
    "LINE/CRE",
    "PLE",
    "Retroposon/SVA",
    "SINE",
    "SINE/5S",
    "SINE/7SL",
    "SINE/Alu",
    "SINE/tRNA",
    "SINE/tRNA-Alu",
    "SINE/tRNA-RTE",
    "SINE/RTE",
    "SINE/Deu",
    "SINE/tRNA-V",
    "SINE/MIR",
    "SINE/U",
    "SINE/tRNA-7SL",
    "SINE/tRNA-CR1",
]

# -----------------------------
# Name normalization (ported from Perl fixName map)
# -----------------------------
NAME_MAP: Dict[str, str] = {
    "DNA/Chompy": "DNA",
    "DNA/CMC-Chapaev": "DNA/CMC",
    "DNA/CMC-Chapaev-3": "DNA/CMC",
    "DNA/CMC-EnSpm": "DNA/CMC",
    "DNA/CMC-Transib": "DNA/CMC",
    "DNA/En-Spm": "DNA/CMC",
    "DNA/Crypton-3": "DNA/Crypton",
    "DNA/PIF-Harbinger": "DNA/Harbinger",
    "DNA/PIF-ISL2EU": "DNA/Harbinger",
    "DNA/Tourist": "DNA/Harbinger",
    "DNA/AcHobo": "DNA/hAT",
    "DNA/Charlie": "DNA/hAT",
    "DNA/Chompy1": "DNA/hAT",
    "DNA/MER1_type": "DNA/hAT",
    "DNA/Tip100": "DNA/hAT",
    "DNA/hAT-Ac": "DNA/hAT",
    "DNA/hAT-Blackjack": "DNA/hAT",
    "DNA/hAT-Charlie": "DNA/hAT",
    "DNA/hAT-Tag1": "DNA/hAT",
    "DNA/hAT-Tip100": "DNA/hAT",
    "DNA/hAT-hATw": "DNA/hAT",
    "DNA/hAT-hobo": "DNA/hAT",
    "DNA/hAT_Tol2": "DNA/hAT",
    "DNA/Kolobok-IS4EU": "DNA/Kolobok",
    "DNA/Kolobok-T2": "DNA/Kolobok",
    "DNA/T2": "DNA/Kolobok",
    "DNA/MULE-MuDR": "DNA/MULE",
    "DNA/MULE-NOF": "DNA/MULE",
    "DNA/MuDR": "DNA/MULE",
    "DNA/piggyBac": "DNA/PiggyBac",
    "DNA/MER2_type": "DNA/TcMar",
    "DNA/Mariner": "DNA/TcMar",
    "DNA/Pogo": "DNA/TcMar",
    "DNA/Stowaway": "DNA/TcMar",
    "DNA/Tc1": "DNA/TcMar",
    "DNA/Tc2": "DNA/TcMar",
    "DNA/Tc4": "DNA/TcMar",
    "DNA/TcMar-Fot1": "DNA/TcMar",
    "DNA/TcMar-ISRm11": "DNA/TcMar",
    "DNA/TcMar-Mariner": "DNA/TcMar",
    "DNA/TcMar-Pogo": "DNA/TcMar",
    "DNA/TcMar-Tc1": "DNA/TcMar",
    "DNA/TcMar-Tc2": "DNA/TcMar",
    "DNA/TcMar-Tigger": "DNA/TcMar",
    "DNA/Tigger": "DNA/TcMar",
    "DNA/Helitron": "RC/Helitron",
    "LTR/DIRS1": "LTR/DIRS",
    "LTR/ERV-Foamy": "LTR/ERVL",
    "LTR/ERV-Lenti": "LTR/ERV",
    "LTR/ERVL-MaLR": "LTR/ERVL",
    "LTR/Gypsy-Troyka": "LTR/Gypsy",
    "LTR/MaLR": "LTR/ERVL",
    "LINE/CR1-Zenon": "LINE/CR1",
    "LINE/I": "LINE/Jockey-I",
    "LINE/Jockey": "LINE/Jockey-I",
    "LINE/L1-Tx1": "LINE/L1",
    "LINE/R2-Hero": "LINE/R2",
    "LINE/RTE-BovB": "LINE/RTE",
    "LINE/RTE-RTE": "LINE/RTE",
    "LINE/RTE-X": "LINE/RTE",
    "LINE/telomeric": "LINE/Jockey-I",
    "LINE/Penlope": "PLE",
    "PLE/Hydra": "PLE",
    "PLE/Naiad": "PLE",
    "PLE/Athena": "PLE",
    "PLE/Poseidon": "PLE",
    "PLE/Nematis": "PLE",
    "PLE/Neptune": "PLE",
    "PLE/Coprina": "PLE",
    "SINE/B2": "SINE/tRNA",
    "SINE/B4": "SINE/tRNA-Alu",
    "SINE/BovA": "SINE/tRNA-RTE",
    "SINE/C": "SINE/tRNA",
    "SINE/Core": "SINE",
    "SINE/ID": "SINE/tRNA",
    "SINE/Lys": "SINE/tRNA",
    "SINE/MERMAID": "SINE/tRNA-V",
    "SINE/RTE-BovB": "SINE/RTE",
    "SINE/tRNA-Glu": "SINE/tRNA",
    "SINE/tRNA-Lys": "SINE/tRNA",
    "SINE/V": "SINE/tRNA-V",
    "Unknown/Y-chromosome": "Unknown",
    "DNA/CMC-3": "DNA/CMC",
    "DNA/CMC-Mirage": "DNA/CMC",
    "DNA/Dada": "DNA/Dada",
    "DNA/Kolobok-Hydra": "DNA/Kolobok",
    "DNA/MULE-F": "DNA/MULE",
    "DNA/TcMar-Stowaway": "DNA/TcMar",
    "DNA/TcMar-Tc4": "DNA/TcMar",
    "DNA/TcMar-m44": "DNA/TcMar",
    "DNA/hAT-Pegasus": "DNA/hAT",
    "DNA/hAT-Tol2": "DNA/hAT",
    "DNA/hAT-hAT1": "DNA/hAT",
    "DNA/hAT-hAT5": "DNA/hAT",
    "DNA/hAT-hAT6": "DNA/hAT",
    "DNA/hAT-hAT19": "DNA/hAT",
    "LINE/Jockey-I-I": "LINE/Jockey-I",
    "LINE/I-Jockey": "LINE/Jockey-I",
    "LINE/R2-NeSL": "LINE/R2",
    "LTR/Copia(Xen1)": "LTR/Copia",
    "LTR/ERV4": "LTR/ERV",
    "LTR/Ginger": "LTR/Gypsy",
    "Retroposon": "Other",
    "Other/Composite": "Other",
    "SINE/5S-Deu-L2": "SINE/5S",
    "SINE/5S-Sauria-RTE": "SINE/5S",
    "SINE/L2": "SINE",
    "SINE/Mermaid": "SINE/tRNA",
    "SINE/R2": "SINE",
    "SINE/U": "SINE/U",
    "SINE/Core-RTE": "SINE/RTE",
    "SINE/tRNA-C": "SINE/tRNA",
    "SINE/tRNA-Core-RTE": "SINE/tRNA",
    "SINE/tRNA-Core": "SINE/tRNA",
    "SINE/tRNA-Deu-CR1": "SINE/Deu",
    "SINE/tRNA-Deu-L2": "SINE/Deu",
    "SINE/tRNA-Deu": "SINE/Deu",
    "SINE/tRNA-Jockey": "SINE/tRNA",
    "SINE/tRNA-L2": "SINE/tRNA",
    "SINE/tRNA-Rex": "SINE/tRNA",
    "SINE/tRNA-Sauria-L2": "SINE/tRNA",
    "SINE/tRNA-Sauria-RTE": "SINE/tRNA",
    "SINE/tRNA-Sauria": "SINE/tRNA",
    "SINE/Sauria": "SINE/tRNA",
    "SINE/tRNA-V-Core-L2": "SINE/tRNA",
    "SINE/tRNAore-RTE": "SINE/tRNA",
    "SINE/tRNAore": "SINE/tRNA",
    "tRNA": "Structural_RNA",
    "scRNA": "Structural_RNA",
    "RNA": "Structural_RNA",
    "snRNA": "Structural_RNA",
    "rRNA": "Structural_RNA",
    "Satellite/centromeric": "Satellite",
    "Satellite/telomeric": "Satellite",
    "Satellite/acromeric": "Satellite",
}

CANONICAL_ALIASES: Dict[str, str] = {
    "Rc/Helitron": "RC/Helitron",
    "RC/Helitron": "RC/Helitron",
}


def fix_name(name: str) -> str:
    """Normalize to canonical RepeatMasker-like class names used in the plot."""
    n = name.replace("?", "").strip()
    n = CANONICAL_ALIASES.get(n, n)
    return NAME_MAP.get(n, n)


# -----------------------------
# Hard-coded CVD-friendly colors (RGBA)
# Based on Okabe–Ito colorblind-safe palette + deterministic variants.
# Only "Unknown" is gray.
# -----------------------------
def _hex_rgba(h: str, a: float = 1.0) -> Tuple[float, float, float, float]:
    h = h.lstrip("#")
    return (int(h[0:2], 16) / 255.0, int(h[2:4], 16) / 255.0, int(h[4:6], 16) / 255.0, a)


def _mix(
    c1: Tuple[float, float, float, float],
    c2: Tuple[float, float, float, float],
    t: float,
) -> Tuple[float, float, float, float]:
    """Linear mix between RGBA colors c1 and c2. t=0 -> c1, t=1 -> c2."""
    return (
        c1[0] * (1 - t) + c2[0] * t,
        c1[1] * (1 - t) + c2[1] * t,
        c1[2] * (1 - t) + c2[2] * t,
        1.0,
    )


WHITE = (1.0, 1.0, 1.0, 1.0)
BLACK = (0.0, 0.0, 0.0, 1.0)

# Okabe–Ito palette (CVD-friendly)
OI_ORANGE = _hex_rgba("#E69F00")
OI_SKYBLUE = _hex_rgba("#56B4E9")
OI_GREEN = _hex_rgba("#009E73")   # bluish green
OI_YELLOW = _hex_rgba("#F0E442")
OI_BLUE = _hex_rgba("#0072B2")
OI_VERMIL = _hex_rgba("#D55E00")  # vermillion
OI_PURPLE = _hex_rgba("#CC79A7")  # reddish purple
OI_GRAY = _hex_rgba("#999999")

UNKNOWN_COLOR = OI_GRAY                    # Only Unknown is gray
OTHER_COLOR = _mix(OI_GRAY, WHITE, 0.45)   # light neutral for "Other"
FALLBACK_SLATE = _hex_rgba("#7A8A99")      # non-gray fallback (rare)

# --- DNA: warm family (orange/vermillion/gold), strong within-family separation ---
DNA_BASE = OI_VERMIL
DNA_GOLD = _mix(OI_ORANGE, OI_YELLOW, 0.35)

DNA_SHADES: Dict[str, Tuple[float, float, float, float]] = {
    "DNA": DNA_BASE,
    "DNA/CMC": _mix(OI_ORANGE, WHITE, 0.15),        # amber
    "DNA/TcMar": DNA_BASE,                          # vermillion
    "DNA/MULE": _mix(DNA_BASE, BLACK, 0.10),        # deeper vermillion
    "DNA/Harbinger": DNA_GOLD,                      # gold (very distinct)
    "DNA/hAT": _mix(OI_ORANGE, BLACK, 0.05),        # orange
    "DNA/PiggyBac": _mix(DNA_BASE, WHITE, 0.20),    # lighter warm
    "DNA/Kolobok": _mix(DNA_BASE, BLACK, 0.15),     # darker warm
    "DNA/Crypton": _mix(DNA_BASE, BLACK, 0.22),
    "DNA/Transib": _mix(DNA_BASE, BLACK, 0.28),
    "DNA/Merlin": _mix(OI_ORANGE, WHITE, 0.30),
    "DNA/Sola": _mix(OI_ORANGE, WHITE, 0.05),
    "DNA/Academ": _mix(DNA_BASE, BLACK, 0.35),
    "DNA/Maverick": _mix(DNA_BASE, BLACK, 0.45),
    "DNA/Dada": _mix(DNA_BASE, WHITE, 0.28),
    "DNA/Ginger": _mix(OI_ORANGE, BLACK, 0.18),
    "DNA/P": _mix(DNA_BASE, WHITE, 0.10),
}

# --- LTR: green family; Copia vs Gypsy made CVD-robust via green vs teal/skyblue mix ---
LTR_SHADES: Dict[str, Tuple[float, float, float, float]] = {
    "LTR": OI_GREEN,
    "LTR/Copia": OI_GREEN,                           # bluish-green
    "LTR/Gypsy": _mix(OI_GREEN, OI_SKYBLUE, 0.55),   # teal/blue-green (distinct)
    "LTR/DIRS": _mix(OI_GREEN, WHITE, 0.25),
    "LTR/Ngaro": _mix(OI_GREEN, OI_SKYBLUE, 0.35),
    "LTR/Pao": _mix(OI_GREEN, BLACK, 0.10),
    "LTR/ERVL": _mix(OI_GREEN, WHITE, 0.40),
    "LTR/ERV1": _mix(OI_GREEN, WHITE, 0.18),
    "LTR/ERV": _mix(OI_GREEN, BLACK, 0.06),
    "LTR/ERVK": _mix(OI_GREEN, BLACK, 0.16),
}

# --- LINE: blue family; separated by lightness ---
LINE_SHADES: Dict[str, Tuple[float, float, float, float]] = {
    "LINE": OI_BLUE,
    "LINE/L1": OI_BLUE,                          # deep blue
    "LINE/RTE": _mix(OI_BLUE, WHITE, 0.35),      # lighter blue
    "LINE/CR1": _mix(OI_BLUE, WHITE, 0.18),
    "LINE/L2": _mix(OI_BLUE, BLACK, 0.10),
    "LINE/R2": _mix(OI_BLUE, BLACK, 0.18),
    "LINE/Jockey-I": _mix(OI_BLUE, BLACK, 0.14),
    "LINE/CRE": _mix(OI_BLUE, WHITE, 0.45),
    "LINE/R1": _mix(OI_BLUE, WHITE, 0.12),
    "LINE/Dong-R4": _mix(OI_BLUE, WHITE, 0.28),
    "LINE/LOA": _mix(OI_BLUE, WHITE, 0.22),
    "LINE/Proto2": _mix(OI_BLUE, WHITE, 0.30),
    "LINE/Rex-Babar": _mix(OI_BLUE, WHITE, 0.40),
}

# --- SINE: purple/pink family (NOT gray), separated by lightness ---
SINE_BASE = OI_PURPLE
SINE_SHADES: Dict[str, Tuple[float, float, float, float]] = {
    "SINE": SINE_BASE,
    "SINE/tRNA": _mix(SINE_BASE, WHITE, 0.18),
    "SINE/Alu": _mix(SINE_BASE, WHITE, 0.32),
    "SINE/7SL": _mix(SINE_BASE, WHITE, 0.40),
    "SINE/5S": _mix(SINE_BASE, BLACK, 0.08),
    "SINE/RTE": _mix(SINE_BASE, BLACK, 0.16),
    "SINE/Deu": _mix(SINE_BASE, WHITE, 0.10),
    "SINE/MIR": _mix(SINE_BASE, WHITE, 0.26),
    "SINE/U": _mix(SINE_BASE, WHITE, 0.46),
    "SINE/tRNA-Alu": _mix(SINE_BASE, WHITE, 0.28),
    "SINE/tRNA-RTE": _mix(SINE_BASE, BLACK, 0.06),
    "SINE/tRNA-V": _mix(SINE_BASE, WHITE, 0.14),
    "SINE/tRNA-7SL": _mix(SINE_BASE, WHITE, 0.36),
    "SINE/tRNA-CR1": _mix(SINE_BASE, WHITE, 0.22),
}

# --- RC/Helitron: brown (distinct from DNA warm tones) ---
RC_SHADES: Dict[str, Tuple[float, float, float, float]] = {
    "RC/Helitron": _hex_rgba("#8C6D31"),
}

# --- PLE / Retroposon: distinct accents, still CVD-friendly ---
OTHER_TE_SHADES: Dict[str, Tuple[float, float, float, float]] = {
    "PLE": _mix(OI_SKYBLUE, OI_PURPLE, 0.35),
    "Retroposon/SVA": _mix(OI_PURPLE, OI_VERMIL, 0.20),
}


def class_color(cls: str) -> Tuple[float, float, float, float]:
    """
    Deterministic color lookup for each class.
    Only "Unknown" is gray.
    """
    if cls == "Unknown":
        return UNKNOWN_COLOR
    if cls == "Other":
        return OTHER_COLOR

    if cls in DNA_SHADES:
        return DNA_SHADES[cls]
    if cls in LTR_SHADES:
        return LTR_SHADES[cls]
    if cls in LINE_SHADES:
        return LINE_SHADES[cls]
    if cls in SINE_SHADES:
        return SINE_SHADES[cls]
    if cls in RC_SHADES:
        return RC_SHADES[cls]
    if cls in OTHER_TE_SHADES:
        return OTHER_TE_SHADES[cls]

    # Broad-category deterministic fallbacks
    if cls.startswith("DNA/") or cls == "DNA":
        return DNA_SHADES["DNA"]
    if cls.startswith("LTR/") or cls == "LTR":
        return LTR_SHADES["LTR"]
    if cls.startswith("LINE/") or cls == "LINE":
        return LINE_SHADES["LINE"]
    if cls.startswith("SINE/") or cls == "SINE":
        return SINE_SHADES["SINE"]
    if cls.startswith("RC/") or cls == "RC/Helitron":
        return RC_SHADES["RC/Helitron"]

    return FALLBACK_SLATE


@dataclass
class DivsumData:
    class_total_bp: Dict[str, int]
    div_percent: Dict[str, Dict[int, float]]  # class -> divbin -> percent
    div_bins: List[int]
    seen_classes: List[str]


def get_genome_size_from_twobit(two_bit: Path) -> int:
    cmd = ["twoBitInfo", "-noNs", str(two_bit), "stdout"]
    try:
        out = subprocess.check_output(cmd, text=True)
    except FileNotFoundError as e:
        raise RuntimeError("twoBitInfo not found in PATH. Install UCSC tools or use --genome-size.") from e
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"twoBitInfo failed: {' '.join(cmd)}") from e

    gsize = 0
    for line in out.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 2 and parts[1].isdigit():
            gsize += int(parts[1])
    if gsize <= 0:
        raise RuntimeError("Genome size computed from twoBitInfo is 0; check input and -noNs behavior.")
    return gsize


def parse_divsum(divsum_path: Path, genome_size: int, max_div: int = 50) -> DivsumData:
    class_total_bp: Dict[str, int] = defaultdict(int)
    div_percent: Dict[str, Dict[int, float]] = defaultdict(lambda: defaultdict(float))
    div_bins_set: Set[int] = set()
    seen_classes: List[str] = []

    in_coverage_section = False
    in_class_section = False
    rec = False

    with divsum_path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line:
                continue

            if line.startswith("Coverage for each repeat"):
                in_class_section = False
                in_coverage_section = True
                continue

            if line.startswith("Class\tRepeat"):
                in_class_section = True
                continue

            if in_class_section and ("\t" in line) and not line.strip().startswith("-"):
                fields = line.split("\t")
                if len(fields) == 5:
                    cls = fix_name(fields[0])
                    try:
                        abs_len = int(fields[2])
                    except ValueError:
                        continue
                    class_total_bp[cls] += abs_len
                else:
                    raise RuntimeError(
                        f"{divsum_path} looks like an old calcDivergenceFromAlign.pl format "
                        f"(expected 5 tab-separated fields in Class section)."
                    )
                continue

            if in_coverage_section and line.startswith("Div "):
                fields = line.split()
                seen_classes = [fix_name(x) for x in fields[1:]]
                rec = True
                continue

            if rec:
                fields = line.split()
                if not fields:
                    continue
                try:
                    div_bin = int(fields[0])
                except ValueError:
                    continue
                if div_bin > max_div:
                    continue

                vals = fields[1:]
                for i, v in enumerate(vals):
                    if i >= len(seen_classes):
                        break
                    cls = seen_classes[i]
                    try:
                        bp = float(v)
                    except ValueError:
                        bp = 0.0
                    div_percent[cls][div_bin] += 100.0 * bp / float(genome_size)

                div_bins_set.add(div_bin)

    div_bins = sorted(div_bins_set)
    return DivsumData(
        class_total_bp=dict(class_total_bp),
        div_percent=div_percent,
        div_bins=div_bins,
        seen_classes=seen_classes,
    )


def read_classes_file(path: Path) -> List[str]:
    """
    Read class labels from a text file.
    - One label per line
    - Blank lines ignored
    - Lines starting with # ignored
    """
    out: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            out.append(line)
    return out


def select_landscape_classes(
    data: DivsumData,
    min_count_bp: int,
    features: Optional[List[str]],
    classes_file: Optional[Path],
) -> List[str]:
    present = set(data.class_total_bp.keys())

    ordered = [c for c in GRAPH_LABELS_ORDER if c in present and data.class_total_bp.get(c, 0) > min_count_bp]
    if not ordered:
        raise RuntimeError("No classes passed filters (min-count / presence). Try lowering --min-count or check divsum.")

    requested: List[str] = []
    if classes_file is not None:
        if not classes_file.exists():
            raise RuntimeError(f"--classes-file not found: {classes_file}")
        requested.extend(read_classes_file(classes_file))
    if features:
        requested.extend(features)

    if requested:
        want = [fix_name(x) for x in requested]
        want_set = set(want)
        filtered = [c for c in ordered if c in want_set]

        missing = [w for w in want if w not in present]
        if missing:
            sys.stderr.write("Warning: requested classes not present in divsum totals:\n")
            for w in sorted(set(missing)):
                sys.stderr.write(f"  - {w}\n")

        if not filtered:
            raise RuntimeError("After applying --features/--classes-file, no classes remain to plot.")
        return filtered

    return ordered


def build_figure_and_axes(figsize: Tuple[float, float], legend_mode: str, panel_ratio: Tuple[float, float]):
    if legend_mode == "panel":
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(1, 2, width_ratios=[panel_ratio[0], panel_ratio[1]], wspace=0.05)
        ax = fig.add_subplot(gs[0, 0])
        ax_leg = fig.add_subplot(gs[0, 1])
        ax_leg.axis("off")
        return fig, ax, ax_leg
    else:
        fig, ax = plt.subplots(figsize=figsize)
        return fig, ax, None


def plot_landscape(
    data: DivsumData,
    classes: List[str],
    out_path: Path,
    title: str,
    figsize: Tuple[float, float],
    dpi: int,
    show_grid: bool,
    ylim: Optional[Tuple[float, float]],
    legend_mode: str,
    legend_right: float,
    legend_ncol: int,
    legend_fontsize: int,
    panel_ratio: Tuple[float, float],
):
    x = data.div_bins
    if not x:
        raise RuntimeError("No divergence bins found to plot (check --max-div).")

    classes = [fix_name(c) for c in classes]

    fig, ax, ax_leg = build_figure_and_axes(figsize, legend_mode, panel_ratio)
    bottoms = [0.0] * len(x)

    for cls in classes:
        color = class_color(cls)
        y = [float(data.div_percent.get(cls, {}).get(b, 0.0)) for b in x]
        ax.bar(
            x,
            y,
            bottom=bottoms,
            width=0.8,
            align="center",
            edgecolor="white",
            linewidth=0.2,
            label=cls,
            color=color,
        )
        bottoms = [b0 + y0 for b0, y0 in zip(bottoms, y)]

    ax.set_title(title)
    ax.set_xlabel("Kimura substitution level (CpG adjusted)")
    ax.set_ylabel("Percent of genome")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)

    if show_grid:
        ax.yaxis.grid(True, linestyle="-", linewidth=0.5, alpha=0.25)
    else:
        ax.grid(False)

    if max(x) - min(x) >= 10:
        ax.set_xticks(list(range(min(x), max(x) + 1, 5)))

    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])

    # Legend: alphabetical (does NOT change stacking order), no title
    handles, labels = ax.get_legend_handles_labels()
    pairs = sorted(zip(labels, handles), key=lambda z: z[0].lower())
    labels_sorted = [p[0] for p in pairs]
    handles_sorted = [p[1] for p in pairs]

    if legend_mode == "panel":
        ax_leg.legend(
            handles_sorted,
            labels_sorted,
            loc="upper left",
            frameon=False,
            ncol=legend_ncol,
            fontsize=legend_fontsize,
        )
    else:
        fig.subplots_adjust(right=legend_right)
        ax.legend(
            handles_sorted,
            labels_sorted,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            frameon=False,
            ncol=legend_ncol,
            fontsize=legend_fontsize,
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format=out_path.suffix.lstrip("."), dpi=dpi)
    plt.close(fig)


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Create an Interspersed Repeat Landscape plot from RepeatMasker *.divsum (SVG/PDF/PNG)."
    )
    p.add_argument("--div", required=True, type=Path, help="Divergence summary file (*.divsum)")

    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--genome-size", type=int, help="Genome size used for percentage calculations")
    g.add_argument("--twoBit", type=Path, help="Genome .2bit file; genome size computed via twoBitInfo -noNs")

    p.add_argument("--title", default="Interspersed Repeat Landscape", help="Plot title")
    p.add_argument("--out", required=True, type=Path, help="Output figure path (.svg, .pdf, .png)")

    p.add_argument("--min-count", type=int, default=0, help="Minimum class bp in totals to include (like Perl -minCount)")

    p.add_argument(
        "--features",
        nargs="+",
        default=None,
        help="Optional list of class labels to plot. Names are normalized (e.g., DNA/En-Spm -> DNA/CMC).",
    )
    p.add_argument(
        "--classes-file",
        type=Path,
        default=None,
        help="TXT file with class labels to plot (one per line; # comments allowed). Names are normalized.",
    )

    p.add_argument("--max-div", type=int, default=50, help="Maximum divergence bin to include (default: 50)")
    p.add_argument("--figsize", type=float, nargs=2, default=(10.0, 6.0), help="Figure size inches: W H")
    p.add_argument("--dpi", type=int, default=300, help="DPI for raster outputs (ignored for SVG/PDF vector)")
    p.add_argument("--grid", action="store_true", help="Show a light horizontal grid on Y")

    p.add_argument("--ylim", type=float, nargs=2, default=None, help="Y-axis limits: ymin ymax (e.g., 0 3.5)")

    p.add_argument(
        "--legend-mode",
        default="margin",
        choices=["margin", "panel"],
        help="Legend layout. 'margin' reserves a right margin; 'panel' uses a separate legend column.",
    )
    p.add_argument(
        "--legend-right",
        type=float,
        default=0.78,
        help="For --legend-mode margin: right boundary for axes (smaller -> more room for legend).",
    )
    p.add_argument("--legend-ncol", type=int, default=1, help="Legend columns.")
    p.add_argument("--legend-fontsize", type=int, default=8, help="Legend font size.")
    p.add_argument(
        "--legend-panel-ratio",
        type=float,
        nargs=2,
        default=(4.2, 1.3),
        help="For --legend-mode panel: width ratios plot legend (e.g., 4.2 1.3)",
    )

    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_argparser().parse_args(argv)

    if not args.div.exists() or args.div.stat().st_size == 0:
        raise SystemExit(f"Error: --div file missing or empty: {args.div}")

    if args.twoBit:
        genome_size = get_genome_size_from_twobit(args.twoBit)
    else:
        genome_size = int(args.genome_size)
        if genome_size <= 0:
            raise SystemExit("Error: --genome-size must be > 0")

    data = parse_divsum(args.div, genome_size=genome_size, max_div=args.max_div)

    classes = select_landscape_classes(
        data=data,
        min_count_bp=args.min_count,
        features=args.features,
        classes_file=args.classes_file,
    )

    present_in_div = set(data.div_percent.keys())
    for cls in sorted(present_in_div):
        if cls not in GRAPH_LABELS_ORDER and cls not in {"Satellite", "Segmental", "Structural_RNA", "Simple_repeat"}:
            sys.stderr.write(f"Warning: class '{cls}' not in canonical label list; it will not be graphed.\n")

    ylim = tuple(args.ylim) if args.ylim is not None else None
    panel_ratio = (float(args.legend_panel_ratio[0]), float(args.legend_panel_ratio[1]))

    plot_landscape(
        data=data,
        classes=classes,
        out_path=args.out,
        title=args.title.replace("_", " "),
        figsize=(float(args.figsize[0]), float(args.figsize[1])),
        dpi=int(args.dpi),
        show_grid=bool(args.grid),
        ylim=ylim,
        legend_mode=args.legend_mode,
        legend_right=float(args.legend_right),
        legend_ncol=int(args.legend_ncol),
        legend_fontsize=int(args.legend_fontsize),
        panel_ratio=panel_ratio,
    )

    sys.stderr.write(f"Wrote: {args.out}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
