#!/usr/bin/env python3
"""
secretion_report.py - turn a combined SecretionGenie summary into a single
self-contained HTML report describing which secretion systems were predicted
in which taxonomic groups.

Designed for the output of a whole-RefSeq screen: the summary is streamed line
by line, so a multi-gigabyte combined.csv with tens of millions of rows is
processed in bounded memory. Nothing but pandas / numpy / matplotlib is needed
and the resulting HTML embeds every figure, so it can be emailed or dropped on
a web server as one file.

Inputs
------
  --summary       combined.csv from the screen (the '<out>/secretiongenie-summary.csv'
                  files concatenated). '####' separator lines are tolerated, as is
                  the presence or absence of the trailing 'seq' column.
  --assembly-info ncbi_assembly_info.tsv (NCBI assembly summary). Used for the
                  accession -> organism / taxid / group / proteome size mapping.
                  Both flavours are handled: RefSeq (GCF_ in column 1) and
                  GenBank (GCA_ in column 1 with the paired GCF_ in column 18).
  --genome-stats  optional genome_stats.tsv[.gz] written by the batch controller
                  (accession, n_proteins, n_systems). Supplying it makes every
                  prevalence a fraction of the genomes actually screened rather
                  than of the genomes that happen to carry a system.
  --taxdump       optional directory holding NCBI's nodes.dmp and names.dmp. With
                  it, genomes are grouped at a real taxonomic rank (phylum by
                  default). Without it, the report falls back to the genus parsed
                  from organism_name plus the coarse 'group' column.
  --lineage       optional TSV of 'accession <tab> lineage' (semicolon-separated,
                  e.g. taxonkit output) as an alternative to --taxdump.

Usage
-----
  python3 secretion_report.py \
      --summary   results/summaries/combined.csv \
      --assembly-info ~/databases/ncbi_assembly_info.tsv \
      --genome-stats  results/genome_stats.tsv.gz \
      --taxdump   ~/databases/taxdump \
      --rank phylum \
      --outdir    results/report

Outputs
-------
  <outdir>/secretion_report.html      the report (self-contained)
  <outdir>/tables/*.tsv               every aggregate behind the figures
"""

from collections import defaultdict, Counter
import argparse
import base64
import gzip
import html
import io
import os
import random
import sys
import textwrap
import time

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap

# matplotlib 3.11 deprecated boxplot(vert=False) in favour of orientation=
_MPL_VER = tuple(int(x) for x in matplotlib.__version__.split(".")[:2])
BOXPLOT_HORIZONTAL = ({"orientation": "horizontal"} if _MPL_VER >= (3, 11)
                      else {"vert": False})

# ---------------------------------------------------------------------------
# look and feel
#
# Figures are built to publication standard: vector PDF/SVG plus a
# high-resolution PNG, with text kept as text (fonttype 42 / svg no-convert)
# so labels stay editable in Illustrator, Inkscape or Affinity.
# ---------------------------------------------------------------------------

INK = "#28251D"
MUTED = "#7A7974"
FAINT = "#BAB9B4"
ACCENT = "#20808D"
ACCENT_DARK = "#12525C"
SECOND = "#A84B2F"
GRID = "#D4D1CA"

# single-hue sequential ramp built from the accent (safe for greyscale printing)
SEQUENTIAL = LinearSegmentedColormap.from_list(
    "mab_teal", ["#FBFBF9", "#DCEAEC", "#AFD2D7", "#74B0B9", "#3D909B", "#20808D", "#12525C"])


def _mix(hex_colour, target, fraction):
    hex_colour = hex_colour.lstrip("#")
    r, g, b = (int(hex_colour[i:i + 2], 16) for i in (0, 2, 4))
    tr, tg, tb = target
    return "#%02x%02x%02x" % (int(r + (tr - r) * fraction),
                              int(g + (tg - g) * fraction),
                              int(b + (tb - b) * fraction))


_BASE8 = ["#20808D", "#A84B2F", "#1B474D", "#8FC7CE",
          "#944454", "#D19900", "#848456", "#6E522B"]
# 24 well-separated categorical colours: the base sequence, then lighter and
# darker derivatives of it, so a 19-system legend stays inside the palette
CATEGORICAL = (_BASE8
               + [_mix(c, (255, 255, 255), 0.45) for c in _BASE8]
               + [_mix(c, (0, 0, 0), 0.35) for c in _BASE8])
# hatch marks back up colour for the derived groups (colour-vision safety)
HATCHES = [None] * 8 + ["///"] * 8 + ["..."] * 8


def cat_colour(i):
    return CATEGORICAL[i % len(CATEGORICAL)]


def cat_hatch(i):
    return HATCHES[i % len(HATCHES)]


plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
    # keep text as text in vector output so it stays editable downstream
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial", "Liberation Sans"],
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.labelsize": 10.5,
    "axes.linewidth": 0.8,
    "axes.edgecolor": GRID,
    "axes.labelcolor": INK,
    "axes.titlecolor": INK,
    "axes.labelpad": 5,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK,
    "ytick.labelcolor": INK,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "grid.alpha": 0.9,
    "legend.frameon": False,
    "legend.fontsize": 9,
    "legend.handlelength": 1.4,
    "legend.handleheight": 0.9,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "figure.constrained_layout.use": False,
})


class FigureSink(object):
    """writes each figure once per requested format and returns a base64 PNG
       for embedding in the HTML report"""

    def __init__(self, outdir=None, formats=("pdf", "svg", "png"), dpi=600, embed_dpi=200):
        self.outdir = outdir
        self.formats = [f.strip().lower() for f in formats if f.strip()]
        self.dpi = dpi
        self.embed_dpi = embed_dpi
        self.written = []
        self._n = 0
        if outdir:
            os.makedirs(outdir, exist_ok=True)

    def emit(self, fig, name):
        self._n += 1
        if self.outdir:
            stem = "fig%02d_%s" % (self._n, name)
            for ext in self.formats:
                path = os.path.join(self.outdir, "%s.%s" % (stem, ext))
                fig.savefig(path, format=ext, dpi=self.dpi, facecolor="white")
                self.written.append(path)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=self.embed_dpi, facecolor="white")
        plt.close(fig)
        return base64.b64encode(buf.getvalue()).decode("ascii")


SINK = FigureSink(outdir=None)


def panel_label(ax, letter, dx=-34, dy=20):
    """the A / B / C panel letters journals ask for, positioned in absolute
       points above and to the left of the axes so they never sit on the title"""
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction",
                xytext=(dx, dy), textcoords="offset points",
                fontsize=14, fontweight="bold", color=INK, va="center", ha="left",
                annotation_clip=False)


def style(ax, xlabel=None, ylabel=None, title=None, grid_axis="y"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(False)
    if grid_axis in ("x", "both"):
        ax.xaxis.grid(True)
    if grid_axis in ("y", "both"):
        ax.yaxis.grid(True)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, loc="left", pad=11)
    return ax


def fig_to_b64(fig, name="figure"):
    return SINK.emit(fig, name)


def smart_open(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


def tick_compact(v, _pos=None):
    """short axis tick labels: 250k rather than 250,000"""
    v = float(v)
    if abs(v) >= 1e6:
        return ("%.1fM" % (v / 1e6)).replace(".0M", "M")
    if abs(v) >= 1e3:
        return ("%.0fk" % (v / 1e3))
    return "%g" % v


def human(n):
    try:
        n = float(n)
    except (TypeError, ValueError):
        return str(n)
    if n >= 1e9:
        return "%.2f B" % (n / 1e9)
    if n >= 1e6:
        return "%.2f M" % (n / 1e6)
    if n >= 1e3:
        return "%s" % format(int(n), ",")
    return "%g" % n


# ===========================================================================
# System key: what the pipeline actually detects, and what it might mean
#
# The distinction matters and is kept explicit everywhere in the report:
#   detected  - the machine components found and the quorum they satisfied.
#               This is a structural claim about gene content and gene order.
#   roles     - what systems of this class are known to do in organisms where
#               they have been studied. This is NOT a claim about the genome in
#               which the system was found: no effector, substrate, expression
#               or phenotype is measured here.
#
# Definitions of the models follow the TXSScan package
# (https://github.com/macsy-models/TXSScan) and its source publications
# (Abby et al. 2016 Sci Rep 6:23080; Denise et al. 2019 PLoS Biol 17:e3000390;
# Bongiovanni et al. 2024 Nat Commun 15:429).
# ===========================================================================

SYSTEM_KEY = {
    "T1SS": {
        "full": "Type I secretion system",
        "family": "ABC-transporter export",
        "detected": "ABC transporter + membrane-fusion protein + outer-membrane factor "
                    "encoded together (abc, mfp, omf)",
        "roles": "One-step export of large exoproteins across both membranes: haemolysins, "
                 "proteases, lipases, adhesins, S-layer and biofilm proteins",
        "markers": "T1SS_abc, T1SS_mfp, T1SS_omf",
        "notes": "omf is treated as a loner and multi-system component, so it may be shared "
                 "between systems or sit away from the rest of the locus",
    },
    "T2SS": {
        "full": "Type II secretion system",
        "family": "Type IV filament superfamily",
        "detected": "Gsp/Xcp secreton: secretin GspD plus the pseudopilus and assembly "
                    "platform (gspC-gspM)",
        "roles": "Two-step secretion of folded periplasmic proteins: toxins, cellulases, "
                 "chitinases, lipases, phosphatases",
        "markers": "T2SS_gspD (secretin), T2SS_gspE (ATPase), T2SS_gspF",
        "notes": "Shares ancestry with type IV pili; discrimination from T4aP/MSH rests on the "
                 "full component set rather than on any single gene",
    },
    "T3SS": {
        "full": "Type III secretion system (injectisome)",
        "family": "Flagellum-related export",
        "detected": "Sct core: export apparatus (sctR, sctS, sctT, sctU, sctV), ATPase sctN, "
                    "C-ring sctQ, secretin sctC",
        "roles": "Contact-dependent injection of effectors directly into eukaryotic host "
                 "cells; central to many host-associated and symbiotic lifestyles",
        "markers": "T3SS_sctC, T3SS_sctN, T3SS_sctV",
        "notes": "Homologous to the flagellar export apparatus; the Flagellum model is "
                 "searched in parallel and flagellar genes are forbidden in this model",
    },
    "T4aP": {
        "full": "Type IVa pilus",
        "family": "Type IV filament superfamily",
        "detected": "PilA fibre with PilB/PilT ATPases, PilC platform, PilQ secretin and the "
                    "PilMNOP alignment complex",
        "roles": "Twitching motility, surface adhesion, natural transformation, "
                 "microcolony and biofilm formation, phage receptor",
        "markers": "T4aP_pilB, T4aP_pilT (retraction), T4aP_pilQ",
        "notes": "Presence of the retraction ATPase pilT is the usual separator from T4bP "
                 "and from Tad",
    },
    "T4bP": {
        "full": "Type IVb pilus",
        "family": "Type IV filament superfamily",
        "detected": "T4bP pilin plus dedicated assembly ATPase and platform",
        "roles": "Adhesion and colonisation; the class includes the toxin-coregulated pilus, "
                 "bundle-forming pilus, longus, Cof and the R64 thin pilus",
        "markers": "T4bP_pilA, T4bP_pilB, T4bP_pilC",
        "notes": "Often plasmid or island encoded, so it is a frequent multi-locus call",
    },
    "Tad": {
        "full": "Tad (tight-adherence) pilus",
        "family": "Type IV filament superfamily",
        "detected": "Flp pre-pilin with TadA ATPase, TadB/TadC platform, TadZ and the "
                    "RcpA/RcpC secretin module",
        "roles": "Tight adherence to surfaces and host tissue, biofilm initiation, "
                 "autoaggregation",
        "markers": "Tad_flp, Tad_tadA, Tad_tadZ",
        "notes": "Widespread and frequently mobile; a common source of Tad calls in taxa "
                 "with no described pilus",
    },
    "MSH": {
        "full": "Mannose-sensitive haemagglutinin pilus",
        "family": "Type IV filament superfamily",
        "detected": "MshA-type pilus with the mshG platform and mshJKL assembly components",
        "roles": "Attachment to abiotic and chitinous surfaces, early biofilm formation, "
                 "phage receptor",
        "markers": "MSH_mshG, MSH_mshJ, MSH_mshL",
        "notes": "Closely related to T4aP; described mainly from Vibrio and relatives",
    },
    "ComM": {
        "full": "Competence machinery of monoderms",
        "family": "Type IV filament superfamily",
        "detected": "ComG pseudopilus with comEA/comEB and the ComEC membrane channel, "
                    "and the ComK regulator where present",
        "roles": "Natural transformation: binding, uptake and translocation of "
                 "environmental DNA into the cytoplasm",
        "markers": "ComM_comEC (channel), ComM_comGB, ComM_comK",
        "notes": "Modelled for monoderms (the classical Bacillus machinery); calls in "
                 "diderm lineages deserve manual inspection",
    },
    "Flagellum": {
        "full": "Bacterial flagellum",
        "family": "Flagellum-related export",
        "detected": "Flagellar export apparatus and basal body core (flgB, flgC, fliE and the "
                    "Sct-homologous export components)",
        "roles": "Swimming and swarming motility; in several taxa also export of "
                 "non-flagellar substrates",
        "markers": "Flg_flgB, Flg_fliE, Flg_sctN_FLG",
        "notes": "The model is designed to separate flagella from T3SS injectisomes, not to "
                 "annotate the complete flagellum",
    },
    "pT4SSt": {
        "full": "Protein-secreting type IV secretion system, MPF_T class",
        "family": "Type IV secretion",
        "detected": "VirB/VirD4-like machinery (virB1-virB10 homologues) in a protein-"
                    "secretion configuration",
        "roles": "Translocation of effector proteins into eukaryotic or bacterial cells; "
                 "the class contains the T4ASS machineries typified by Agrobacterium VirB "
                 "and the Helicobacter Cag apparatus",
        "markers": "T4SS_T_virB4, T4SS_T_virB10, T4SS_T_virB2",
        "notes": "Conjugative T4SS models were moved to CONJscan; only protein-secreting "
                 "variants are searched here, so a call is not evidence of conjugation",
    },
    "pT4SSi": {
        "full": "Protein-secreting type IV secretion system, MPF_I class",
        "family": "Type IV secretion",
        "detected": "IncI-like tra machinery (traE, traP, traY homologues) in a protein-"
                    "secretion configuration",
        "roles": "Effector translocation by the type IVB branch, which is the class "
                 "containing the Dot/Icm machineries of Legionella and Coxiella",
        "markers": "T4SS_I_traE, T4SS_I_traP, T4SS_I_traY",
        "notes": "Type IVB systems are closely related to IncI plasmid conjugation systems, "
                 "hence the MPF_I label",
    },
    "T5aSS": {
        "full": "Type Va secretion: classical autotransporter",
        "family": "Autotransporter",
        "detected": "A single polypeptide carrying the C-terminal autotransporter "
                    "beta-barrel translocator domain (PF03797)",
        "roles": "Surface display or release of adhesins, proteases, esterases and "
                 "cytotoxins; adhesion, aggregation, serum resistance",
        "markers": "T5aSS_PF03797",
        "notes": "A one-gene model, so counts scale with the number of autotransporter genes "
                 "in a genome rather than with the number of machines",
    },
    "T5bSS": {
        "full": "Type Vb secretion: two-partner secretion",
        "family": "Autotransporter",
        "detected": "TpsB transporter of the ShlB/FhaC/HecB family, with its TpsA partner "
                    "where co-localised",
        "roles": "Secretion of large exoproteins such as filamentous haemagglutinins, "
                 "haemolysins and contact-dependent growth inhibition toxins",
        "markers": "T5bSS_PF03865",
        "notes": "Frequently present in several copies per genome",
    },
    "T5cSS": {
        "full": "Type Vc secretion: trimeric autotransporter",
        "family": "Autotransporter",
        "detected": "YadA-like trimeric autotransporter anchor domain (PF03895)",
        "roles": "Trimeric surface adhesins mediating binding to collagen, fibronectin and "
                 "epithelial surfaces; autoaggregation and serum resistance",
        "markers": "T5cSS_PF03895",
        "notes": "Also a one-gene model; short anchor domains give low sequence coverage, so "
                 "profile coverage is computed on the profile, not the protein",
    },
    "T6SSi": {
        "full": "Type VI secretion system, subtype i",
        "family": "Contractile injection system",
        "detected": "Full phage-tail-like machine: TssBC sheath, Hcp/TssD tube, VgrG/TssI "
                    "spike, TssJLM membrane complex, baseplate and ClpV/TssH ATPase",
        "roles": "Contact-dependent injection of toxins into neighbouring bacteria or "
                 "eukaryotic cells; interbacterial antagonism, niche competition, and in "
                 "some taxa metal-ion acquisition",
        "markers": "T6SSi_tssM, T6SSi_tssH, T6SSi_tssD",
        "notes": "The broadly distributed subtype; multiple copies per genome are common and "
                 "usually correspond to genuinely distinct machines",
    },
    "T6SSii": {
        "full": "Type VI secretion system, subtype ii",
        "family": "Contractile injection system",
        "detected": "The Francisella pathogenicity island machinery (iglA/iglB/iglC with "
                    "pdpA-pdpD)",
        "roles": "Intracellular survival: phagosome escape and cytosolic replication in "
                 "the lineages where it has been characterised",
        "markers": "T6SSii_iglA, T6SSii_iglC, T6SSii_pdpA",
        "notes": "Very narrow known distribution, so calls outside Francisella and relatives "
                 "are worth inspecting individually",
    },
    "T6SSiii": {
        "full": "Type VI secretion system, subtype iii",
        "family": "Contractile injection system",
        "detected": "The Bacteroidota-type machine, including its distinct membrane complex "
                    "as redefined in 2024",
        "roles": "Antagonism between gut and environmental Bacteroidota; strain competition "
                 "and community structuring",
        "markers": "T6SSiii_tssD, T6SSiii_tssH, T6SSiii_tssP",
        "notes": "Model updated per Bongiovanni et al. 2024; largely restricted to "
                 "Bacteroidota",
    },
    "T9SS": {
        "full": "Type IX secretion system (Por secretion system)",
        "family": "Bacteroidota-specific export",
        "detected": "GldK/GldL/GldM/GldN core with PorV and the SprA translocon, plus the "
                    "SprE/SprT loners",
        "roles": "Gliding motility and surface delivery of adhesins and peptidases, "
                 "including the gingipains of periodontal pathogens",
        "markers": "T9SS_gldK_TIGR03525, T9SS_sprA_PF14349, T9SS_porV",
        "notes": "sprA, sprE, sprT and porQ are loners in the model, so they frequently sit "
                 "far from the gld core",
    },
    "Archaeal-T4P": {
        "full": "Archaeal type IV pili superfamily",
        "family": "Type IV filament superfamily",
        "detected": "arCOG-defined archaeal pilus machinery: the assembly ATPase, membrane "
                    "platform, prepilin peptidase and pilin set",
        "roles": "Motility via the archaellum, UV-inducible DNA exchange through Ups pili, "
                 "substrate binding by the bindosome, surface adhesion and biofilm",
        "markers": "Archaeal-T4P_arCOG01817_ATP, Archaeal-T4P_arCOG01809_IM",
        "notes": "The model spans several distinct archaeal machines, so a call identifies "
                 "the superfamily rather than one specific appendage",
    },
}


# ---------------------------------------------------------------------------
# Functional grouping. Calling all nineteen models "secretion systems" hides the
# fact that they do very different things: only some of them inject effectors
# into another cell, and one of them imports DNA rather than exporting anything.
# ---------------------------------------------------------------------------

FUNCTIONAL_CATEGORY = {
    # deliver proteins/effectors directly into a host or competing cell
    "T3SS": "Direct effector delivery",
    "pT4SSt": "Direct effector delivery",
    "pT4SSi": "Direct effector delivery",
    "T6SSi": "Direct effector delivery",
    "T6SSii": "Direct effector delivery",
    "T6SSiii": "Direct effector delivery",
    # release toxins, enzymes and adhesins into the surroundings or the surface
    "T1SS": "Extracellular protein secretion",
    "T2SS": "Extracellular protein secretion",
    "T5aSS": "Extracellular protein secretion",
    "T5bSS": "Extracellular protein secretion",
    "T5cSS": "Extracellular protein secretion",
    "T9SS": "Extracellular protein secretion",
    # surface attachment, motility, biofilm, colonisation
    "T4aP": "Adhesion, motility, colonisation",
    "T4bP": "Adhesion, motility, colonisation",
    "Tad": "Adhesion, motility, colonisation",
    "MSH": "Adhesion, motility, colonisation",
    "Flagellum": "Adhesion, motility, colonisation",
    "Archaeal-T4P": "Adhesion, motility, colonisation",
    # import rather than export
    "ComM": "DNA uptake and competence",
}

CATEGORY_ORDER = ["Direct effector delivery", "Extracellular protein secretion",
                  "Adhesion, motility, colonisation", "DNA uptake and competence"]

CATEGORY_BLURB = {
    "Direct effector delivery": "inject proteins into host or competing cells",
    "Extracellular protein secretion": "release toxins, enzymes and adhesins outside the cell",
    "Adhesion, motility, colonisation": "attachment, twitching or swimming motility, biofilm",
    "DNA uptake and competence": "natural transformation, DNA acquisition",
}

CATEGORY_COLOUR = {
    "Direct effector delivery": "#A84B2F",
    "Extracellular protein secretion": "#20808D",
    "Adhesion, motility, colonisation": "#848456",
    "DNA uptake and competence": "#7A39BB",
}


def category_of(system):
    return FUNCTIONAL_CATEGORY.get(system, "Unclassified model")


# the classes whose defining activity is delivering effectors into another cell:
# an unexpected occurrence of one of these is the kind of result worth chasing
EFFECTOR_DELIVERY = tuple(s for s, c in FUNCTIONAL_CATEGORY.items()
                          if c == "Direct effector delivery")

# models that are expected in one domain only
EXPECTED_DOMAIN = {"Archaeal-T4P": "archaea"}


def key_for(system):
    return SYSTEM_KEY.get(system, {
        "full": system, "family": "not in key", "detected": "-", "roles": "-",
        "markers": "-", "notes": "no key entry for this model name",
    })


# ===========================================================================
# 1. taxonomy
# ===========================================================================

RANKS = ["superkingdom", "domain", "phylum", "class", "order", "family", "genus", "species"]


def load_taxdump(taxdump_dir, wanted_rank):
    """parse nodes.dmp / names.dmp and return taxid -> (rank_name, lineage_dict)"""
    nodes_path = os.path.join(taxdump_dir, "nodes.dmp")
    names_path = os.path.join(taxdump_dir, "names.dmp")
    if not (os.path.isfile(nodes_path) and os.path.isfile(names_path)):
        sys.stderr.write("  taxdump: nodes.dmp / names.dmp not found in %s, skipping\n" % taxdump_dir)
        return None, None

    parent = {}
    rank = {}
    with open(nodes_path, errors="replace") as handle:
        for line in handle:
            fields = line.split("\t|\t")
            if len(fields) < 3:
                continue
            try:
                tid = int(fields[0])
                parent[tid] = int(fields[1])
            except ValueError:
                continue
            rank[tid] = fields[2].replace("\t|", "").strip()

    name = {}
    with open(names_path, errors="replace") as handle:
        for line in handle:
            fields = line.split("\t|\t")
            if len(fields) < 4:
                continue
            if fields[3].replace("\t|", "").strip() != "scientific name":
                continue
            try:
                name[int(fields[0])] = fields[1].replace("\t|", "").strip()
            except ValueError:
                continue

    sys.stderr.write("  taxdump: %s nodes, %s names\n" % (human(len(parent)), human(len(name))))

    cache = {}

    def lineage_of(taxid):
        """walk to the root, collecting the wanted ranks"""
        if taxid in cache:
            return cache[taxid]
        chain = []
        tid = taxid
        seen = set()
        while tid and tid not in seen and tid != 1:
            seen.add(tid)
            chain.append(tid)
            tid = parent.get(tid, 0)
        out = {}
        for tid in chain:
            r = rank.get(tid, "no rank")
            if r in RANKS and r not in out:
                out[r] = name.get(tid, "unclassified")
        cache[taxid] = out
        return out

    return lineage_of, wanted_rank


def load_lineage_tsv(path, wanted_rank):
    """accession <tab> semicolon-separated lineage (e.g. taxonkit lineage output)"""
    mapping = {}
    with smart_open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            acc = fields[0].strip()
            parts = [p.strip() for p in fields[1].split(";") if p.strip()]
            mapping[acc] = parts
    sys.stderr.write("  lineage file: %s accessions\n" % human(len(mapping)))
    # position in a standard 7-rank lineage
    idx = {"domain": 0, "superkingdom": 0, "phylum": 1, "class": 2,
           "order": 3, "family": 4, "genus": 5, "species": 6}.get(wanted_rank, 1)
    return {acc: (parts[idx] if len(parts) > idx and parts[idx] else "unclassified")
            for acc, parts in mapping.items()}


def load_assembly_info(path, lineage_of=None, rank="phylum", lineage_map=None):
    """accession -> dict(organism, genus, group, taxid, n_proteins, taxon)"""
    info = {}
    n_rows = 0
    with smart_open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 25:
                continue
            n_rows += 1
            acc = f[0] if f[0].startswith("GCF_") else (f[17] if len(f) > 17 and f[17].startswith("GCF_") else "")
            if not acc:
                continue
            organism = f[7] if len(f) > 7 else "na"
            genus = organism.split()[0] if organism and organism != "na" else "unclassified"
            if genus.startswith("[") or genus.startswith("'"):
                genus = genus.strip("[]'")
            group = f[24] if len(f) > 24 and f[24] not in ("", "na") else "unclassified"
            try:
                taxid = int(f[5])
            except (ValueError, IndexError):
                taxid = 0
            try:
                n_prot = int(f[35]) if len(f) > 35 and f[35].isdigit() else 0
            except (ValueError, IndexError):
                n_prot = 0

            taxon = None
            if lineage_map is not None:
                taxon = lineage_map.get(acc)
            if taxon is None and lineage_of is not None and taxid:
                taxon = lineage_of(taxid).get(rank)
            if not taxon:
                # fall back: genus is always available, group gives the domain
                taxon = genus if rank == "genus" else ("%s (no %s)" % (group, rank) if rank != "group" else group)

            info[acc] = {"organism": organism, "genus": genus, "group": group,
                         "taxid": taxid, "n_proteins": n_prot, "taxon": taxon}
    sys.stderr.write("  assembly info: %s rows, %s GCF accessions mapped\n"
                     % (human(n_rows), human(len(info))))
    return info


# ===========================================================================
# 2. streaming pass over the combined summary
# ===========================================================================

class Aggregates(object):
    """everything the report needs, accumulated in bounded memory"""

    RESERVOIR = 40000

    def __init__(self):
        self.n_rows = 0
        self.n_separators = 0
        self.n_systems = 0
        self.systems_per_type = Counter()             # system -> number of systems
        self.genomes_per_type = defaultdict(set)      # system -> set of accessions (bounded below)
        self.per_genome = defaultdict(Counter)        # accession -> {system: n}
        self.hit_types = Counter()                    # cluster / loner / multi_system ...
        self.gene_status = Counter()
        self.component_counts = defaultdict(Counter)  # system -> {function: n systems containing it}
        self.hits_per_system = defaultdict(list)      # system -> reservoir of hit counts
        self.wholeness = defaultdict(list)            # system -> reservoir of wholeness
        self.score = defaultdict(list)                # system -> reservoir of score
        self.loci = defaultdict(Counter)              # system -> {nb_loci: n}
        self.replicons_seen = 0
        self.bad_rows = 0
        self._rng = random.Random(1)
        # notable findings: the best few instances of each system per genus, so
        # rare-but-convincing occurrences can be surfaced after prevalence is known
        self.best_per_group = defaultdict(list)   # (group, system) -> [record, ...]
        self.cross_domain = defaultdict(list)     # system -> [record, ...] wrong domain
        self.group_of = {}                        # accession -> grouping label (genus)
        self.keep_per_group = 3

    def _reservoir(self, store, key, value):
        buf = store[key]
        if len(buf) < self.RESERVOIR:
            buf.append(value)
        else:
            j = self._rng.randrange(len(buf) + 1)
            if j < self.RESERVOIR:
                buf[j] = value

    def close_system(self, system, accession, functions, n_hits, wholeness, score, nb_loci,
                     group=None, domain=None):
        if system is None:
            return
        self.n_systems += 1
        self.systems_per_type[system] += 1
        self.per_genome[accession][system] += 1
        for func in functions:
            self.component_counts[system][func] += 1
        self._reservoir(self.hits_per_system, system, n_hits)
        if wholeness is not None:
            self._reservoir(self.wholeness, system, wholeness)
        if score is not None:
            self._reservoir(self.score, system, score)
        self.loci[system][nb_loci] += 1

        # keep the strongest few instances per (grouping taxon, system)
        record = (round(wholeness, 4) if wholeness is not None else 0.0,
                  round(score, 2) if score is not None else 0.0,
                  n_hits, nb_loci, accession)
        if group is not None:
            self.group_of[accession] = group
            bucket = self.best_per_group[(group, system)]
            if len(bucket) < self.keep_per_group:
                bucket.append(record)
                bucket.sort(reverse=True)
            elif record > bucket[-1]:
                bucket[-1] = record
                bucket.sort(reverse=True)

        # a model expected in one domain, found in the other
        expected = EXPECTED_DOMAIN.get(system)
        if domain and expected and domain != expected and domain != "unclassified":
            if len(self.cross_domain[system]) < 400:
                self.cross_domain[system].append(record + (domain,))
        elif domain == "archaea" and expected is None and system != "ComM":
            if len(self.cross_domain[system]) < 400:
                self.cross_domain[system].append(record + (domain,))


def stream_summary(path, agg, max_rows=None, progress_every=2_000_000, meta=None):
    """single pass over the combined summary.

    Rows belonging to one predicted system are contiguous (SecretionGenie writes
    one block per system, closed by a '####' line), so systems are delimited by a
    change of system_id. The '####' lines are counted independently and used as a
    consistency check on the number of systems found.
    """
    t0 = time.time()
    meta = meta or {}
    with smart_open(path) as handle:
        header_line = handle.readline()
        if not header_line:
            raise SystemExit("empty summary file: %s" % path)
        header = header_line.rstrip("\n").split(",")
        col = {name: i for i, name in enumerate(header)}
        required = ["file", "system", "system_id", "function", "gene_status",
                    "hit_type", "system_wholeness", "system_score", "nb_loci"]
        missing = [c for c in required if c not in col]
        if missing:
            raise SystemExit("summary is missing expected column(s): %s\nheader was: %s"
                             % (", ".join(missing), ",".join(header)))
        n_fields = len(header)

        cur_id = None
        cur_sys = None
        cur_acc = None
        cur_group = None
        cur_domain = None
        cur_funcs = set()
        cur_hits = 0
        cur_whole = None
        cur_score = None
        cur_loci = 1

        for line in handle:
            if line[0] == "#":
                agg.n_separators += 1
                continue
            f = line.rstrip("\n").split(",")
            if len(f) < n_fields - 1:
                agg.bad_rows += 1
                continue
            agg.n_rows += 1
            if max_rows and agg.n_rows > max_rows:
                break
            if progress_every and agg.n_rows % progress_every == 0:
                sys.stderr.write("    %s rows in %.0fs\n" % (human(agg.n_rows), time.time() - t0))
                sys.stderr.flush()

            sid = f[col["system_id"]]
            if sid != cur_id:
                agg.close_system(cur_sys, cur_acc, cur_funcs, cur_hits,
                                 cur_whole, cur_score, cur_loci,
                                 group=cur_group, domain=cur_domain)
                cur_id = sid
                cur_sys = f[col["system"]]
                cur_acc = f[col["file"]]
                rec = meta.get(cur_acc)
                cur_group = rec["genus"] if rec else "unclassified"
                cur_domain = rec["group"] if rec else "unclassified"
                cur_funcs = set()
                cur_hits = 0
                try:
                    cur_whole = float(f[col["system_wholeness"]])
                except ValueError:
                    cur_whole = None
                try:
                    cur_score = float(f[col["system_score"]])
                except ValueError:
                    cur_score = None
                try:
                    cur_loci = int(f[col["nb_loci"]])
                except ValueError:
                    cur_loci = 1

            cur_hits += 1
            cur_funcs.add(f[col["function"]])
            agg.hit_types[f[col["hit_type"]]] += 1
            agg.gene_status[f[col["gene_status"]]] += 1

        agg.close_system(cur_sys, cur_acc, cur_funcs, cur_hits,
                         cur_whole, cur_score, cur_loci,
                         group=cur_group, domain=cur_domain)

    sys.stderr.write("  streamed %s rows / %s systems / %s separators in %.0fs\n"
                     % (human(agg.n_rows), human(agg.n_systems),
                        human(agg.n_separators), time.time() - t0))
    if agg.bad_rows:
        sys.stderr.write("  %s malformed rows skipped\n" % human(agg.bad_rows))
    if agg.n_separators and abs(agg.n_separators - agg.n_systems) > max(10, 0.001 * agg.n_systems):
        sys.stderr.write("  WARNING: %s '####' separators but %s systems counted - the file may "
                         "interleave systems; counts are based on system_id changes\n"
                         % (human(agg.n_separators), human(agg.n_systems)))
    return agg


# ===========================================================================
# 3. figures
# ===========================================================================

def fig_systems_overview(order, systems_per_type, genomes_carrying, n_screened):
    fig, axes = plt.subplots(1, 2, figsize=(11, max(3.2, 0.32 * len(order))))
    y = np.arange(len(order))

    counts = [systems_per_type.get(s, 0) for s in order]
    axes[0].barh(y, counts, color=ACCENT, height=0.68)
    axes[0].set_yticks(y, order)
    axes[0].invert_yaxis()
    style(axes[0], xlabel="systems predicted", title="Predicted systems by type", grid_axis="x")
    span = max(counts) if counts else 1
    for yi, v in zip(y, counts):
        axes[0].text(v + span * 0.012, yi, human(v), va="center", fontsize=7.5, color=MUTED)
    axes[0].set_xlim(0, span * 1.16)

    pct = [100.0 * genomes_carrying.get(s, 0) / n_screened if n_screened else 0 for s in order]
    axes[1].barh(y, pct, color=SECOND, height=0.68)
    axes[1].set_yticks(y, [])
    axes[1].tick_params(axis="y", length=0)
    axes[1].invert_yaxis()
    style(axes[1], xlabel="% of screened genomes carrying \u22651", grid_axis="x",
          title="Prevalence across genomes")
    span2 = max(pct) if pct else 1
    for yi, v in zip(y, pct):
        axes[1].text(v + span2 * 0.012, yi, "%.1f%%" % v, va="center", fontsize=7.5, color=MUTED)
    axes[1].set_xlim(0, max(1.0, span2 * 1.18))
    panel_label(axes[0], "A")
    panel_label(axes[1], "B", dx=-22)
    fig.tight_layout()
    return fig_to_b64(fig, "systems_by_type")


def fig_per_genome_distribution(per_genome_counters, n_zero):
    """per_genome_counters: list of Counter(system -> n) for genomes with >=1 system"""
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.4))
    totals = [sum(c.values()) for c in per_genome_counters]
    distinct_counts = [len(c) for c in per_genome_counters]
    values = np.array(totals) if totals else np.array([0])
    top = int(np.percentile(values, 99.5)) if values.size else 1
    top = max(top, 5)
    bins = np.arange(-0.5, top + 1.5, 1)
    allv = np.concatenate([np.zeros(n_zero, dtype=int), values]) if n_zero else values
    axes[0].hist(np.clip(allv, 0, top), bins=bins, color=ACCENT, edgecolor="white", linewidth=0.4)
    style(axes[0], xlabel="systems per genome (clipped at %d)" % top, ylabel="genomes",
          title="Systems per genome")

    distinct = np.array(distinct_counts) if distinct_counts else np.array([0])
    alld = np.concatenate([np.zeros(n_zero, dtype=int), distinct]) if n_zero else distinct
    dmax = int(alld.max()) if alld.size else 1
    axes[1].hist(alld, bins=np.arange(-0.5, dmax + 1.5, 1), color=SECOND,
                 edgecolor="white", linewidth=0.4)
    style(axes[1], xlabel="distinct system types per genome", ylabel="genomes",
          title="System-type richness per genome")
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")
    fig.tight_layout()
    return fig_to_b64(fig, "systems_per_genome")


def fig_taxon_heatmap(matrix, taxa, systems, denom, title, value_label):
    fig, ax = plt.subplots(figsize=(min(1.05 + 0.52 * len(systems), 13),
                                    max(3.0, 0.34 * len(taxa) + 1.4)))
    data = np.array(matrix, dtype=float)
    im = ax.imshow(data, aspect="auto", cmap=SEQUENTIAL, vmin=0,
                   vmax=max(1.0, np.nanpercentile(data, 99) if data.size else 1))
    ax.set_xticks(np.arange(len(systems)), systems, rotation=45, ha="right")
    labels = ["%s  (n=%s)" % (t, human(denom.get(t, 0))) for t in taxa]
    ax.set_yticks(np.arange(len(taxa)), labels)
    ax.grid(False)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.set_title(title, loc="left", pad=10)
    cbar = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.015)
    cbar.set_label(value_label, fontsize=9)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=9, width=0.8, length=3)
    if data.size and len(taxa) * len(systems) <= 400:
        hi = im.norm(np.nanmax(data)) if np.nanmax(data) else 1
        for i in range(len(taxa)):
            for j in range(len(systems)):
                v = data[i, j]
                if v >= 0.5:
                    ax.text(j, i, "%.0f" % v if v >= 10 else "%.1f" % v,
                            ha="center", va="center", fontsize=6.5,
                            color="white" if im.norm(v) > 0.6 * hi else INK)
    fig.tight_layout()
    return fig_to_b64(fig, "prevalence_by_taxon")


def fig_composition(taxa, systems, matrix, denom):
    """stacked share of system types within each taxon"""
    data = np.array(matrix, dtype=float)
    totals = data.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    frac = 100.0 * data / totals
    fig, ax = plt.subplots(figsize=(11, max(3.0, 0.34 * len(taxa) + 1.2)))
    y = np.arange(len(taxa))
    left = np.zeros(len(taxa))
    for j, s in enumerate(systems):
        ax.barh(y, frac[:, j], left=left, height=0.7, label=s,
                color=cat_colour(j), hatch=cat_hatch(j),
                edgecolor="white", linewidth=0.5)
        left += frac[:, j]
    ax.set_yticks(y, ["%s  (%s systems)" % (t, human(denom.get(t, 0))) for t in taxa])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    style(ax, xlabel="% of all systems in that taxon", grid_axis="x",
          title="System-type composition by taxon")
    ax.legend(ncol=min(6, len(systems)), loc="upper center",
              bbox_to_anchor=(0.5, -0.16 if len(taxa) > 6 else -0.28), fontsize=9)
    fig.tight_layout()
    return fig_to_b64(fig, "composition_by_taxon")


def fig_cooccurrence(per_genome_sets, systems):
    n = len(systems)
    idx = {s: i for i, s in enumerate(systems)}
    both = np.zeros((n, n))
    alone = np.zeros(n)
    for sset in per_genome_sets:
        present = [idx[s] for s in sset if s in idx]
        for i in present:
            alone[i] += 1
        for a in range(len(present)):
            for b in range(a + 1, len(present)):
                both[present[a], present[b]] += 1
                both[present[b], present[a]] += 1
    with np.errstate(divide="ignore", invalid="ignore"):
        cond = 100.0 * both / alone[:, None]
    cond = np.nan_to_num(cond)
    np.fill_diagonal(cond, np.nan)

    fig, ax = plt.subplots(figsize=(min(1.6 + 0.5 * n, 11), min(1.4 + 0.46 * n, 9.5)))
    im = ax.imshow(cond, cmap=SEQUENTIAL, vmin=0, vmax=100)
    ax.set_xticks(np.arange(n), systems, rotation=45, ha="right")
    ax.set_yticks(np.arange(n), systems)
    ax.grid(False)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.set_title("Co-occurrence: % of genomes with the row system that also carry the column system",
                 loc="left", pad=10, fontsize=10)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("% of row-system genomes", fontsize=9)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=9, width=0.8, length=3)
    if n <= 22:
        for i in range(n):
            for j in range(n):
                if i == j or np.isnan(cond[i, j]) or cond[i, j] < 1:
                    continue
                ax.text(j, i, "%.0f" % cond[i, j], ha="center", va="center", fontsize=7,
                        color="white" if cond[i, j] > 60 else INK)
    fig.tight_layout()
    return fig_to_b64(fig, "cooccurrence")


def fig_quality(order, wholeness, score, hits_per_system):
    fig, axes = plt.subplots(1, 3, figsize=(12.5, max(3.0, 0.28 * len(order) + 1.2)))
    for ax, store, label, colour in (
            (axes[0], wholeness, "system wholeness", ACCENT),
            (axes[1], score, "system score", SECOND),
            (axes[2], hits_per_system, "components per system", "#D19900")):
        data = [store.get(s, []) or [0] for s in order]
        bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False,
                        medianprops=dict(color="white", linewidth=1.2),
                        **BOXPLOT_HORIZONTAL)
        for patch in bp["boxes"]:
            patch.set_facecolor(colour)
            patch.set_edgecolor(colour)
            patch.set_alpha(0.85)
        for whisker in bp["whiskers"] + bp["caps"]:
            whisker.set_color(MUTED)
            whisker.set_linewidth(0.8)
        ax.set_yticks(np.arange(1, len(order) + 1),
                      order if ax is axes[0] else [""] * len(order))
        if ax is not axes[0]:
            ax.tick_params(axis="y", length=0)
        ax.invert_yaxis()
        style(ax, xlabel=label, grid_axis="x")
    for letter, ax in zip("ABC", axes):
        panel_label(ax, letter, dx=-34 if ax is axes[0] else -22)
    fig.tight_layout()
    # figure-level title: a three-panel title is too long to sit on one axes
    # without running into the panel letters
    fig.suptitle("Model completeness, score and size of the predicted systems",
                 x=0.0, y=1.08, ha="left", fontsize=12.5, fontweight="bold", color=INK)
    return fig_to_b64(fig, "system_quality")


def fig_components(system, comp_counts, n_systems_of_type, max_genes=26):
    items = comp_counts.most_common(max_genes)
    if not items:
        return None
    genes = [g for g, _ in items][::-1]
    pct = [100.0 * c / n_systems_of_type for _, c in items][::-1]
    fig, ax = plt.subplots(figsize=(7.2, max(2.4, 0.26 * len(genes) + 0.9)))
    y = np.arange(len(genes))
    ax.barh(y, pct, color=ACCENT, height=0.68)
    ax.set_yticks(y, genes)
    ax.set_xlim(0, 105)
    style(ax, xlabel="%% of %s systems containing the component" % human(n_systems_of_type),
          title="%s component prevalence" % system, grid_axis="x")
    for yi, v in zip(y, pct):
        ax.text(v + 1.2, yi, "%.0f" % v, va="center", fontsize=8.5, color=MUTED)
    fig.tight_layout()
    return fig_to_b64(fig, "components_%s" % system.replace("/", "-"))


def fig_loci(order, loci):
    fig, ax = plt.subplots(figsize=(8.6, max(2.6, 0.3 * len(order) + 1.0)))
    y = np.arange(len(order))
    single, multi = [], []
    for s in order:
        c = loci.get(s, Counter())
        tot = sum(c.values()) or 1
        single.append(100.0 * c.get(1, 0) / tot)
        multi.append(100.0 * sum(v for k, v in c.items() if k and k > 1) / tot)
    ax.barh(y, single, color="#C9D9DB", height=0.68, label="single locus")
    ax.barh(y, multi, left=single, color=ACCENT, height=0.68, label="two or more loci")
    ax.set_yticks(y, order)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    style(ax, xlabel="% of systems", grid_axis="x",
          title="Scattered (multi-locus) versus single-locus systems")
    ax.legend(loc="lower right")
    fig.tight_layout()
    return fig_to_b64(fig, "single_vs_multi_locus")


def fig_proteome_vs_systems(pairs):
    if len(pairs) < 50:
        return None
    prot = np.array([p for p, _ in pairs], dtype=float)
    nsys = np.array([s for _, s in pairs], dtype=float)
    keep = (prot > 200) & (prot < np.percentile(prot, 99.5))
    prot, nsys = prot[keep], nsys[keep]
    if prot.size < 50:
        return None
    bins = np.linspace(prot.min(), prot.max(), 26)
    which = np.digitize(prot, bins) - 1
    xs, ys, lo, hi = [], [], [], []
    for b in range(len(bins) - 1):
        sel = nsys[which == b]
        if sel.size < 20:
            continue
        xs.append(0.5 * (bins[b] + bins[b + 1]))
        ys.append(sel.mean())
        lo.append(np.percentile(sel, 25))
        hi.append(np.percentile(sel, 75))
    if not xs:
        return None
    fig, ax = plt.subplots(figsize=(7.6, 3.4))
    ax.fill_between(xs, lo, hi, color=ACCENT, alpha=0.18, linewidth=0)
    ax.plot(xs, ys, color=ACCENT, linewidth=2)
    style(ax, xlabel="predicted proteins in the genome", ylabel="systems per genome",
          title="Secretion-system load versus proteome size")
    ax.text(0.99, 0.04, "line = mean, band = interquartile range", transform=ax.transAxes,
            ha="right", fontsize=9, color=MUTED)
    fig.tight_layout()
    return fig_to_b64(fig, "systems_vs_proteome")


def fig_system_key(order):
    """the legend key: every model, its swatch, and the functional tier it sits in"""
    rows = []   # (kind, label, colour_index, category)
    for cat in list(CATEGORY_ORDER) + ["Unclassified model"]:
        if cat == "Unclassified model":
            members = [sy for sy in order if category_of(sy) not in CATEGORY_ORDER]
        else:
            members = [sy for sy in order if category_of(sy) == cat]
        if not members:
            continue
        rows.append(("HEADER", cat, None, cat))
        rows.append(("BLURB", CATEGORY_BLURB.get(cat, ""), None, cat))
        for sy in members:
            rows.append(("ITEM", sy, order.index(sy), cat))

    # split into two columns; if a tier is cut in half, repeat its header
    split = (len(rows) + 2) // 2
    while split < len(rows) and rows[split][0] == "BLURB":
        split += 1
    left, right = rows[:split], rows[split:]
    if right and right[0][0] == "ITEM":
        right = [("HEADER", right[0][3] + " (continued)", None, right[0][3])] + right
    per_col = max(len(left), len(right))

    row_h = 0.30
    fig, ax = plt.subplots(figsize=(13.2, per_col * row_h + 0.75))
    ax.set_xlim(0, 2)
    ax.set_ylim(0, per_col)
    ax.axis("off")

    for col, column in enumerate((left, right)):
        x0 = col * 1.0
        for row, (kind, label, idx, cat) in enumerate(column):
            y = per_col - row - 0.5
            if kind == "HEADER":
                ax.text(x0 + 0.005, y, label.upper(), fontsize=10, fontweight="bold",
                        color=CATEGORY_COLOUR.get(cat, MUTED), va="center", ha="left")
            elif kind == "BLURB":
                ax.text(x0 + 0.005, y + 0.10, label, fontsize=8.8, color=MUTED,
                        va="center", ha="left", style="italic")
            else:
                ax.add_patch(plt.Rectangle((x0 + 0.014, y - 0.10), 0.048, 0.20,
                                           facecolor=cat_colour(idx), hatch=cat_hatch(idx),
                                           edgecolor="white", linewidth=0.6, clip_on=False))
                ax.text(x0 + 0.080, y, label, fontsize=10, fontweight="bold", color=INK,
                        va="center", ha="left")
                ax.text(x0 + 0.285, y, key_for(label)["full"], fontsize=9.5, color=INK,
                        va="center", ha="left")
    ax.set_title("Key to the nineteen models, grouped by what the machine does",
                 loc="left", pad=14)
    fig.tight_layout()
    return fig_to_b64(fig, "system_key")


def fig_by_category(order, systems_per_type, genomes_carrying, n_screened, cat_carriers):
    """the functional tiers: how much of the catalogue is effector delivery, and
       how many genomes carry at least one system of each tier"""
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.3),
                             gridspec_kw={"width_ratios": [1.25, 1.0]})

    cats = [c for c in CATEGORY_ORDER if any(category_of(sy) == c for sy in order)]
    # panel A: stacked contribution of each system to its tier
    y = np.arange(len(cats))
    left = np.zeros(len(cats))
    for sy in order:
        cat = category_of(sy)
        if cat not in cats:
            continue
        row = cats.index(cat)
        vals = np.zeros(len(cats))
        vals[row] = systems_per_type.get(sy, 0)
        axes[0].barh(y, vals, left=left, height=0.62, color=cat_colour(order.index(sy)),
                     hatch=cat_hatch(order.index(sy)), edgecolor="white", linewidth=0.5)
        left += vals
    axes[0].set_yticks(y, ["%s\n%s" % (c, CATEGORY_BLURB.get(c, "")) for c in cats])
    axes[0].invert_yaxis()
    style(axes[0], xlabel="systems predicted (stacked by model)", grid_axis="x",
          title="What the catalogue is made of")
    for row, c in enumerate(cats):
        total = sum(systems_per_type.get(sy, 0) for sy in order if category_of(sy) == c)
        axes[0].text(total + left.max() * 0.012, row, human(total), va="center",
                     fontsize=9, color=MUTED)
    axes[0].set_xlim(0, left.max() * 1.16)
    axes[0].xaxis.set_major_formatter(mticker.FuncFormatter(tick_compact))
    axes[0].xaxis.set_major_locator(mticker.MaxNLocator(6))

    # panel B: genomes carrying at least one system of the tier
    pct = [100.0 * cat_carriers.get(c, 0) / max(1, n_screened) for c in cats]
    axes[1].barh(y, pct, height=0.62,
                 color=[CATEGORY_COLOUR.get(c, ACCENT) for c in cats])
    axes[1].set_yticks(y, [])
    axes[1].tick_params(axis="y", length=0)
    axes[1].invert_yaxis()
    style(axes[1], xlabel="% of screened genomes with \u22651", grid_axis="x",
          title="Genomes carrying each tier")
    for row, v in enumerate(pct):
        axes[1].text(v + max(pct) * 0.015, row, "%.1f%%" % v, va="center",
                     fontsize=9.5, color=MUTED)
    axes[1].set_xlim(0, max(1.0, max(pct) * 1.2))
    panel_label(axes[0], "A", dx=-150)
    panel_label(axes[1], "B", dx=-22)
    fig.tight_layout()
    return fig_to_b64(fig, "functional_categories")


def fig_category_by_taxon(taxa, cats, matrix, denom, rank):
    """tier-level prevalence per taxon: the compact version of the big heatmap"""
    fig, ax = plt.subplots(figsize=(min(3.4 + 1.7 * len(cats), 12),
                                   max(3.0, 0.36 * len(taxa) + 1.5)))
    data = np.array(matrix, dtype=float)
    im = ax.imshow(data, aspect="auto", cmap=SEQUENTIAL, vmin=0, vmax=100)
    ax.set_xticks(np.arange(len(cats)),
                  [textwrap.fill(c, 15) for c in cats],
                  rotation=0, ha="center", fontsize=9)
    ax.set_yticks(np.arange(len(taxa)),
                  ["%s  (n=%s)" % (t, human(denom.get(t, 0))) for t in taxa])
    ax.grid(False)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    for i in range(len(taxa)):
        for j in range(len(cats)):
            v = data[i, j]
            ax.text(j, i, "%.0f" % v if v >= 1 else "", ha="center", va="center",
                    fontsize=8.5, color="white" if v > 55 else INK)
    ax.set_title("Prevalence of each functional tier by %s" % rank, loc="left", pad=11)
    cbar = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.015)
    cbar.set_label("% of genomes with \u22651 system", fontsize=9)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=9, width=0.8, length=3)
    fig.tight_layout()
    return fig_to_b64(fig, "category_by_taxon")


def fig_notable(records, top_n=18):
    """rare-but-convincing occurrences: the leads worth checking by hand.

    The informative axis is rarity, not completeness: every candidate has already
    passed the completeness filter, so completeness is shown as an annotation and
    the horizontal position says how unusual the occurrence is inside its genus."""
    recs = records[:top_n]
    if not recs:
        return None
    fig, ax = plt.subplots(figsize=(12.6, max(3.2, 0.44 * len(recs) + 1.4)))
    y = np.arange(len(recs))[::-1]
    prev = [max(r["prevalence"], 0.004) for r in recs]
    colours = [CATEGORY_COLOUR.get(category_of(r["system"]), ACCENT) for r in recs]
    lo = min(prev) / 2.2
    ax.set_xscale("log")
    ax.hlines(y, lo, prev, color=GRID, linewidth=1.6, zorder=1)
    ax.scatter(prev, y, s=95, color=colours, zorder=3, edgecolor="white", linewidth=0.9)
    ax.set_yticks(y, ["%s  \u2014  %s" % (r["organism"], r["system"]) for r in recs])
    hi = max(prev)
    ax.set_xlim(lo, hi * 18)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, _: ("%g%%" % v) if 0.01 <= v <= hi * 2.5 else ""))
    style(ax, xlabel="carriers as % of that genus (log scale) \u2014 further left is rarer",
          grid_axis="x",
          title="Rare, high-confidence occurrences: candidate leads for manual review")
    for yi, r, x in zip(y, recs, prev):
        ax.text(x * 1.5, yi, "%s of %s %s  \u00b7  %.0f%% complete, %d components"
                % (human(r["carriers"]), human(r["group_genomes"]), r["group"],
                   100 * r["wholeness"], r["components"]),
                va="center", fontsize=8.6, color=MUTED)
    handles = [plt.Line2D([0], [0], marker="o", linestyle="", markersize=9,
                          markerfacecolor=CATEGORY_COLOUR[c], markeredgecolor="white",
                          label=c)
               for c in CATEGORY_ORDER if any(category_of(r["system"]) == c for r in recs)]
    fig.tight_layout(rect=(0, 0.055 if len(recs) > 6 else 0.11, 1, 1))
    if handles:
        fig.legend(handles=handles, loc="lower center", ncol=min(4, len(handles)),
                   fontsize=9, frameon=False)
    return fig_to_b64(fig, "notable_occurrences")


def fig_repertoire(top_genomes):
    """the genomes carrying the widest range of machines, split by functional tier"""
    if not top_genomes:
        return None
    fig, ax = plt.subplots(figsize=(11.8, max(3.0, 0.42 * len(top_genomes) + 1.3)))
    y = np.arange(len(top_genomes))[::-1]
    left = np.zeros(len(top_genomes))
    for cat in CATEGORY_ORDER:
        vals = np.array([g["types_by_tier"].get(cat, 0) for g in top_genomes], dtype=float)
        ax.barh(y, vals, left=left, height=0.62, color=CATEGORY_COLOUR[cat],
                edgecolor="white", linewidth=0.6, label=cat)
        left += vals
    ax.set_yticks(y, [g["organism"] for g in top_genomes])
    style(ax, xlabel="distinct system types in one genome (stacked by tier)", grid_axis="x",
          title="Widest repertoires in the screen")
    for yi, g, tot in zip(y, top_genomes, left):
        ax.text(tot + left.max() * 0.015, yi, "%d types  \u00b7  %d systems"
                % (g["n_types"], g["n_systems"]), va="center", fontsize=8.8, color=MUTED)
    ax.set_xlim(0, left.max() * 1.34)
    handles, labels = ax.get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0.055 if len(top_genomes) > 6 else 0.10, 1, 1))
    fig.legend(handles, labels, loc="lower center", ncol=min(4, len(labels)),
               fontsize=9, frameon=False)
    return fig_to_b64(fig, "widest_repertoires")


# ===========================================================================
# 4. HTML assembly
# ===========================================================================

CSS = """
:root {
  --bg: #F7F6F2; --surface: #FFFFFF; --border: #D4D1CA;
  --ink: #28251D; --muted: #7A7974; --accent: #20808D; --accent-dark: #12525C;
}
* { box-sizing: border-box; }
body {
  margin: 0; background: var(--bg); color: var(--ink);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
  font-size: 16px; line-height: 1.6;
}
header.hero {
  background: var(--accent-dark); color: #F4F7F7; padding: 40px 44px 34px;
}
header.hero h1 { margin: 0 0 6px; font-size: 30px; line-height: 1.2; letter-spacing: -0.01em; }
header.hero p { margin: 0; color: #BCE2E7; font-size: 15px; }
nav.toc {
  position: sticky; top: 0; z-index: 5; background: var(--surface);
  border-bottom: 1px solid var(--border); padding: 12px 44px; font-size: 14px;
}
nav.toc a { color: var(--accent); text-decoration: none; margin-right: 22px; }
nav.toc a:hover { text-decoration: underline; }
main { max-width: 1180px; margin: 0 auto; padding: 8px 44px 80px; }
section { margin: 44px 0 0; }
h2 {
  font-size: 22px; margin: 0 0 8px; padding-bottom: 8px;
  border-bottom: 2px solid var(--accent); display: inline-block;
}
h3 { font-size: 17px; margin: 34px 0 6px; }
p.sub { color: var(--muted); font-size: 14.5px; max-width: 68em; margin: 6px 0 18px; }
p.sub strong { color: var(--ink); }
figure { margin: 22px 0; background: var(--surface); border: 1px solid var(--border);
         border-radius: 8px; padding: 16px; }
figure img { width: 100%; height: auto; display: block; }
figcaption { color: var(--muted); font-size: 13.5px; margin-top: 10px; line-height: 1.5; }
.cards { display: flex; flex-wrap: wrap; gap: 14px; margin: 20px 0 8px; }
.card {
  flex: 1 1 170px; background: var(--surface); border: 1px solid var(--border);
  border-radius: 8px; padding: 16px 18px;
}
.card .v { font-size: 27px; font-weight: 700; font-variant-numeric: tabular-nums lining-nums; }
.card .l { color: var(--muted); font-size: 12.5px; text-transform: uppercase;
           letter-spacing: 0.06em; margin-top: 4px; }
.tablewrap { overflow-x: auto; background: var(--surface); border: 1px solid var(--border);
             border-radius: 8px; margin: 18px 0; }
table { border-collapse: collapse; width: 100%; font-size: 13.5px; }
th, td { text-align: left; padding: 9px 14px; border-bottom: 1px solid var(--border);
         vertical-align: top; }
th { background: #F1EFE9; font-weight: 600; white-space: nowrap; position: sticky; top: 0; }
td.num { text-align: right; font-variant-numeric: tabular-nums lining-nums; white-space: nowrap; }
tr:last-child td { border-bottom: none; }
tr:hover td { background: #FBFBF9; }
code { background: #F1EFE9; padding: 1px 5px; border-radius: 4px; font-size: 13px; }
pre { background: var(--surface); border: 1px solid var(--border); border-radius: 8px;
      padding: 14px 16px; overflow-x: auto; font-size: 13px; }
a { color: var(--accent); }
footer { color: var(--muted); font-size: 13px; border-top: 1px solid var(--border);
         margin-top: 50px; padding-top: 16px; }
ul.refs { font-size: 13.5px; color: var(--muted); }
ul.refs li { margin-bottom: 5px; }
"""


def figure_html(b64, caption):
    if not b64:
        return ""
    return ('<figure><img src="data:image/png;base64,%s" alt="%s"/>'
            '<figcaption>%s</figcaption></figure>' % (b64, html.escape(caption), caption))


def table_html(headers, rows, numeric_cols=()):
    out = ['<div class="tablewrap"><table><thead><tr>']
    for h in headers:
        out.append("<th>%s</th>" % html.escape(str(h)))
    out.append("</tr></thead><tbody>")
    for row in rows:
        out.append("<tr>")
        for i, cell in enumerate(row):
            cls = ' class="num"' if i in numeric_cols else ""
            out.append("<td%s>%s</td>" % (cls, html.escape(str(cell))))
        out.append("</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


def card(value, label):
    return '<div class="card"><div class="v">%s</div><div class="l">%s</div></div>' % (
        html.escape(str(value)), html.escape(str(label)))


def write_tsv(path, headers, rows):
    with open(path, "w") as out:
        out.write("\t".join(str(h) for h in headers) + "\n")
        for row in rows:
            out.write("\t".join("" if c is None else str(c) for c in row) + "\n")


# ===========================================================================
# 5. main
# ===========================================================================

def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="build an HTML report from a combined SecretionGenie summary",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--summary", required=True,
                    help="combined secretiongenie-summary.csv (optionally .gz)")
    ap.add_argument("--assembly-info", required=True,
                    help="NCBI assembly summary TSV used for the screen")
    ap.add_argument("--genome-stats",
                    help="genome_stats.tsv[.gz] from the batch controller: accession, "
                         "n_proteins, n_systems. Gives correct prevalence denominators")
    ap.add_argument("--taxdump", help="directory with NCBI nodes.dmp and names.dmp")
    ap.add_argument("--lineage", help="TSV of accession <tab> semicolon-separated lineage")
    ap.add_argument("--rank", default="phylum",
                    help="taxonomic rank used for grouping in the taxonomy figures")
    ap.add_argument("--outdir", default="secretion_report",
                    help="output directory for the HTML, figures/ and tables/")
    ap.add_argument("--max-rows", type=int, help="stop after this many data rows (testing)")
    ap.add_argument("--top-taxa", type=int, default=14,
                    help="how many taxa to show in the taxonomy figures")
    ap.add_argument("--component-figures", type=int, default=6,
                    help="component-prevalence figures for the N most abundant systems")
    ap.add_argument("--title", default="Secretion systems across the RefSeq prokaryotes",
                    help="report title")
    ap.add_argument("--rare-max-carriers", type=float, default=0.02,
                    help="a system counts as rare in a genus when at most this fraction of "
                         "that genus carries it")
    ap.add_argument("--rare-min-genomes", type=int, default=25,
                    help="minimum genomes sequenced for a genus before rarity means anything")
    ap.add_argument("--rare-min-wholeness", type=float, default=0.8,
                    help="minimum model completeness for a rare occurrence to be reported")
    ap.add_argument("--notable-top", type=int, default=60,
                    help="how many candidate occurrences to keep in the report table")
    ap.add_argument("--figure-formats", default="pdf,svg,png",
                    help="publication figure formats written to <outdir>/figures "
                         "(comma separated; 'none' to skip)")
    ap.add_argument("--dpi", type=int, default=600, help="resolution of the written figures")
    ap.add_argument("--embed-dpi", type=int, default=200,
                    help="resolution of the PNGs embedded in the HTML")
    return ap.parse_args(argv)


def main(argv=None):
    global SINK
    args = parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    tables_dir = os.path.join(args.outdir, "tables")
    os.makedirs(tables_dir, exist_ok=True)
    formats = [] if args.figure_formats.strip().lower() in ("none", "") else \
        args.figure_formats.split(",")
    figures_dir = os.path.join(args.outdir, "figures") if formats else None
    SINK = FigureSink(outdir=figures_dir, formats=formats, dpi=args.dpi,
                      embed_dpi=args.embed_dpi)

    # ---------------- taxonomy ----------------------------------------------
    sys.stderr.write("[1/6] taxonomy\n")
    lineage_of = None
    lineage_map = None
    if args.lineage:
        lineage_map = load_lineage_tsv(args.lineage, args.rank)
    elif args.taxdump:
        lineage_of, _rank_used = load_taxdump(args.taxdump, args.rank)
    else:
        sys.stderr.write("  no --taxdump / --lineage: grouping falls back to genus and domain\n")

    # ---------------- assembly info and screened set ------------------------
    sys.stderr.write("[2/6] assembly info\n")
    info = load_assembly_info(args.assembly_info, lineage_of=lineage_of,
                              rank=args.rank, lineage_map=lineage_map)

    screened = None
    proteome_pairs = []
    stats_systems = {}
    stats_proteins = {}
    if args.genome_stats:
        screened = set()
        with smart_open(args.genome_stats) as handle:
            for line in handle:
                if line.startswith("#") or line.lower().startswith("accession"):
                    continue
                f = line.rstrip("\n").split("\t")
                if not f or not f[0]:
                    continue
                screened.add(f[0])
                try:
                    n_prot = int(f[1])
                except (IndexError, ValueError):
                    n_prot = 0
                try:
                    n_sys = int(f[2])
                except (IndexError, ValueError):
                    n_sys = 0
                stats_systems[f[0]] = n_sys
                stats_proteins[f[0]] = n_prot
                if n_prot:
                    proteome_pairs.append((n_prot, n_sys))
        sys.stderr.write("  genome_stats: %s genomes screened\n" % human(len(screened)))

    # ---------------- streaming pass ----------------------------------------
    sys.stderr.write("[3/6] streaming %s\n" % args.summary)
    agg = stream_summary(args.summary, Aggregates(), max_rows=args.max_rows, meta=info)

    if screened is None:
        screened = set(agg.per_genome)
        sys.stderr.write("  no --genome-stats: denominators use the %s genomes with at least "
                         "one system, so prevalences are upper bounds\n" % human(len(screened)))
    n_screened = max(1, len(screened))
    n_with = len(agg.per_genome)
    n_zero = max(0, len(screened) - n_with)

    if not proteome_pairs:
        for acc, counter in agg.per_genome.items():
            n_prot = info.get(acc, {}).get("n_proteins", 0)
            if n_prot:
                proteome_pairs.append((n_prot, sum(counter.values())))

    # ---------------- aggregate by taxon ------------------------------------
    sys.stderr.write("[4/6] aggregating by %s\n" % args.rank)
    order = [s for s, _ in agg.systems_per_type.most_common()]
    genomes_carrying = Counter()
    for acc, counter in agg.per_genome.items():
        for s in counter:
            genomes_carrying[s] += 1

    taxon_of = {}
    for acc in screened:
        rec = info.get(acc)
        taxon_of[acc] = rec["taxon"] if rec else "not in assembly info"
    taxon_genomes = Counter(taxon_of.values())
    top_taxa = [t for t, _ in taxon_genomes.most_common(args.top_taxa)]

    # taxon x system: genomes carrying, and systems counted
    carriers_matrix = defaultdict(Counter)
    counts_matrix = defaultdict(Counter)
    for acc, counter in agg.per_genome.items():
        t = taxon_of.get(acc, "not in assembly info")
        for s, n in counter.items():
            carriers_matrix[t][s] += 1
            counts_matrix[t][s] += n
    prevalence = [[100.0 * carriers_matrix[t].get(s, 0) / max(1, taxon_genomes[t])
                   for s in order] for t in top_taxa]
    composition = [[counts_matrix[t].get(s, 0) for s in order] for t in top_taxa]
    systems_per_taxon = {t: sum(counts_matrix[t].values()) for t in top_taxa}

    # ---------------- functional tiers --------------------------------------
    cat_carriers = Counter()          # tier -> genomes carrying >=1
    cat_systems = Counter()           # tier -> systems
    for acc, counter in agg.per_genome.items():
        for t in {category_of(sy) for sy in counter}:
            cat_carriers[t] += 1
        for sy, n in counter.items():
            cat_systems[category_of(sy)] += n

    cats_present = [c for c in CATEGORY_ORDER if any(category_of(sy) == c for sy in order)]
    cat_carriers_by_taxon = defaultdict(Counter)
    for acc, counter in agg.per_genome.items():
        t = taxon_of.get(acc, "not in assembly info")
        for c in {category_of(sy) for sy in counter}:
            cat_carriers_by_taxon[t][c] += 1
    cat_prevalence = [[100.0 * cat_carriers_by_taxon[t].get(c, 0) / max(1, taxon_genomes[t])
                       for c in cats_present] for t in top_taxa]

    # ---------------- notable findings --------------------------------------
    sys.stderr.write("[5/6] finding rare, high-confidence occurrences\n")
    genus_genomes = Counter()
    for acc in screened:
        rec = info.get(acc)
        genus_genomes[rec["genus"] if rec else "unclassified"] += 1
    genus_carriers = defaultdict(Counter)
    for acc, counter in agg.per_genome.items():
        rec = info.get(acc)
        g = rec["genus"] if rec else "unclassified"
        for sy in counter:
            genus_carriers[g][sy] += 1

    notable = []
    for (group, system), bucket in agg.best_per_group.items():
        n_group = genus_genomes.get(group, 0)
        if n_group < args.rare_min_genomes or group == "unclassified":
            continue
        carriers = genus_carriers[group].get(system, 0)
        if carriers > max(1, args.rare_max_carriers * n_group):
            continue
        wholeness, score, n_comp, nb_loci, acc = bucket[0]
        if wholeness < args.rare_min_wholeness:
            continue
        rec = info.get(acc, {})
        notable.append({
            "accession": acc, "organism": rec.get("organism", acc), "group": group,
            "group_genomes": n_group, "carriers": carriers, "system": system,
            "category": category_of(system), "wholeness": wholeness, "score": score,
            "components": n_comp, "loci": nb_loci,
            "spotlight": system in EFFECTOR_DELIVERY,
            "prevalence": 100.0 * carriers / max(1, n_group),
        })
    notable.sort(key=lambda r: (not r["spotlight"], -r["wholeness"], -r["score"],
                                r["prevalence"], -r["group_genomes"]))
    notable = notable[:args.notable_top]
    sys.stderr.write("  %d candidate occurrences (rare in genus, completeness \u2265 %.2f)\n"
                     % (len(notable), args.rare_min_wholeness))

    repertoires = []
    for acc, counter in agg.per_genome.items():
        types_by_tier = Counter()
        for sy in counter:
            types_by_tier[category_of(sy)] += 1
        repertoires.append({
            "types_by_tier": types_by_tier,
            "accession": acc,
            "organism": info.get(acc, {}).get("organism", acc),
            "n_types": len(counter),
            "n_systems": sum(counter.values()),
            "tiers": "%d/4 tiers" % len(types_by_tier),
            "systems": ", ".join("%s x%d" % (sy, n) for sy, n in counter.most_common()),
        })
    repertoires.sort(key=lambda r: (-r["n_types"], -r["n_systems"]))
    top_repertoires = repertoires[:15]

    crossover = []
    for system, recs in agg.cross_domain.items():
        for wholeness, score, n_comp, nb_loci, acc, domain in recs:
            if wholeness < args.rare_min_wholeness:
                continue
            crossover.append({
                "system": system, "accession": acc,
                "organism": info.get(acc, {}).get("organism", acc),
                "domain": domain, "wholeness": wholeness, "score": score,
                "components": n_comp,
            })
    crossover.sort(key=lambda r: (-r["wholeness"], -r["score"]))

    # ---------------- figures ------------------------------------------------
    sys.stderr.write("[6/6] figures and HTML\n")
    figs = {}
    figs["overview"] = fig_systems_overview(order, agg.systems_per_type,
                                            genomes_carrying, n_screened)
    figs["per_genome"] = fig_per_genome_distribution(list(agg.per_genome.values()), n_zero)
    figs["key"] = fig_system_key(order)
    figs["categories"] = fig_by_category(order, agg.systems_per_type, genomes_carrying,
                                         n_screened, cat_carriers)
    figs["category_taxon"] = fig_category_by_taxon(top_taxa, cats_present, cat_prevalence,
                                                   taxon_genomes, args.rank)
    figs["notable"] = fig_notable(notable)
    figs["repertoire"] = fig_repertoire(top_repertoires)
    figs["heatmap"] = fig_taxon_heatmap(
        prevalence, top_taxa, order, taxon_genomes,
        "Prevalence of each system by %s" % args.rank,
        "% of genomes with \u22651 system")
    figs["composition"] = fig_composition(top_taxa, order, composition, systems_per_taxon)
    figs["cooccurrence"] = fig_cooccurrence(
        [set(c) for c in agg.per_genome.values()], order[:min(len(order), 20)])
    figs["quality"] = fig_quality(order, agg.wholeness, agg.score, agg.hits_per_system)
    figs["loci"] = fig_loci(order, agg.loci)
    figs["proteome"] = fig_proteome_vs_systems(proteome_pairs)
    component_figs = []
    for system in order[:args.component_figures]:
        b64 = fig_components(system, agg.component_counts[system],
                             agg.systems_per_type[system])
        if b64:
            component_figs.append((system, b64))

    # ---------------- tables -------------------------------------------------
    key_rows = []
    for cat in CATEGORY_ORDER:
        for sy in [x for x in order if category_of(x) == cat]:
            k = key_for(sy)
            key_rows.append([sy, k["full"], cat, k["detected"], k["roles"], k["markers"],
                             format(agg.systems_per_type[sy], ","),
                             "%.2f" % (100.0 * genomes_carrying[sy] / n_screened)])
    # stable column order: grouped by tier, alphabetical inside a tier, so the
    # layout of this file does not shift between runs or datasets
    system_cols = []
    for cat in CATEGORY_ORDER:
        system_cols.extend(sorted(sy for sy in order if category_of(sy) == cat))
    system_cols.extend(sorted(sy for sy in order if category_of(sy) not in CATEGORY_ORDER))

    genome_rows = []
    for acc in sorted(screened):
        rec = info.get(acc, {})
        counter = agg.per_genome.get(acc)
        n_sys = sum(counter.values()) if counter else stats_systems.get(acc, 0)
        n_prot = stats_proteins.get(acc) or rec.get("n_proteins", 0)
        ranks = {}
        taxid = rec.get("taxid", 0)
        if lineage_of is not None and taxid:
            ranks = lineage_of(taxid)
        row = [acc, n_sys, n_prot, len(counter) if counter else 0]
        # one column per functional tier, then one per system type
        for cat in CATEGORY_ORDER:
            row.append(sum(n for sy, n in (counter or {}).items()
                           if category_of(sy) == cat))
        for sy in system_cols:
            row.append((counter or {}).get(sy, 0))
        for r in RANKS:
            value = ranks.get(r, "")
            if not value:
                # fall back to what the assembly summary itself provides
                if r in ("superkingdom", "domain"):
                    value = rec.get("group", "")
                elif r == "genus":
                    value = rec.get("genus", "")
            row.append(value or "NA")
        row.append(taxid or "NA")
        row.append(rec.get("organism", "NA"))
        genome_rows.append(row)
    tier_cols = ["n_" + c.split(",")[0].split(" and ")[0].strip().lower().replace(" ", "_")
                 for c in CATEGORY_ORDER]
    write_tsv(os.path.join(tables_dir, "per_genome_summary.tsv"),
              ["accession", "n_systems", "n_proteins", "n_system_types"] + tier_cols +
              system_cols + RANKS + ["taxid", "organism_name"],
              genome_rows)
    sys.stderr.write("  per-genome table: %s rows in %s\n"
                     % (human(len(genome_rows)),
                        os.path.join(tables_dir, "per_genome_summary.tsv")))

    write_tsv(os.path.join(tables_dir, "system_key.tsv"),
              ["system", "full_name", "functional_category", "what_was_detected",
               "potential_roles", "marker_components", "n_systems", "pct_genomes"], key_rows)

    write_tsv(os.path.join(tables_dir, "systems_by_type.tsv"),
              ["system", "functional_category", "n_systems", "genomes_carrying",
               "pct_of_screened_genomes"],
              [[s, category_of(s), agg.systems_per_type[s], genomes_carrying[s],
                "%.3f" % (100.0 * genomes_carrying[s] / n_screened)] for s in order])

    write_tsv(os.path.join(tables_dir, "prevalence_by_taxon.tsv"),
              [args.rank, "genomes"] + order,
              [[t, taxon_genomes[t]] + ["%.2f" % v for v in row]
               for t, row in zip(top_taxa, prevalence)])

    write_tsv(os.path.join(tables_dir, "systems_by_taxon.tsv"),
              [args.rank, "genomes"] + order,
              [[t, taxon_genomes[t]] + [counts_matrix[t].get(s, 0) for s in order]
               for t in top_taxa])

    write_tsv(os.path.join(tables_dir, "notable_occurrences.tsv"),
              ["accession", "organism", "genus", "genus_genomes", "genus_carriers",
               "pct_of_genus", "system", "functional_category", "wholeness", "score",
               "components", "loci", "effector_delivery"],
              [[r["accession"], r["organism"], r["group"], r["group_genomes"], r["carriers"],
                "%.2f" % r["prevalence"], r["system"], r["category"],
                "%.3f" % r["wholeness"], "%.2f" % r["score"], r["components"], r["loci"],
                "yes" if r["spotlight"] else "no"] for r in notable])

    write_tsv(os.path.join(tables_dir, "widest_repertoires.tsv"),
              ["accession", "organism", "n_system_types", "n_systems", "tiers", "systems"],
              [[r["accession"], r["organism"], r["n_types"], r["n_systems"], r["tiers"],
                r["systems"]] for r in repertoires[:500]])

    if crossover:
        write_tsv(os.path.join(tables_dir, "cross_domain_occurrences.tsv"),
                  ["system", "accession", "organism", "genome_domain", "wholeness",
                   "score", "components"],
                  [[r["system"], r["accession"], r["organism"], r["domain"],
                    "%.3f" % r["wholeness"], "%.2f" % r["score"], r["components"]]
                   for r in crossover[:500]])

    comp_rows = []
    for system in order:
        total = agg.systems_per_type[system]
        for func, n in agg.component_counts[system].most_common():
            comp_rows.append([system, func, n, "%.2f" % (100.0 * n / max(1, total))])
    write_tsv(os.path.join(tables_dir, "component_prevalence.tsv"),
              ["system", "component", "n_systems_with_component", "pct_of_systems"], comp_rows)

    # ---------------- HTML ---------------------------------------------------
    sections = []
    mean_per_carrier = agg.n_systems / max(1, n_with)
    cards = "".join([
        card(human(n_screened), "genomes screened"),
        card(human(agg.n_systems), "systems predicted"),
        card(len(order), "system types"),
        card("%.1f%%" % (100.0 * n_with / n_screened), "genomes with \u22651 system"),
        card("%.1f" % mean_per_carrier, "systems per carrier"),
    ])

    overview_rows = [[s, category_of(s), format(agg.systems_per_type[s], ","),
                      format(genomes_carrying[s], ","),
                      "%.2f" % (100.0 * genomes_carrying[s] / n_screened)] for s in order]

    sections.append("""
<section id="overview">
  <h2>Overview</h2>
  <p class="sub">%s rows of the combined summary were streamed, delimiting %s predicted
  systems across %s screened genomes. A predicted system is a set of components that satisfied
  its model's quorum and co-localisation rules; component hits that failed those rules are in the
  per-batch rejected-candidate files, not here.</p>
  <div class="cards">%s</div>
  %s
  %s
  %s
</section>""" % (
        human(agg.n_rows), human(agg.n_systems), human(n_screened), cards,
        figure_html(figs["overview"],
                    "A: predicted systems by type. B: the share of screened genomes carrying "
                    "at least one system of that type."),
        figure_html(figs["per_genome"],
                    "A: how many systems a single genome carries. B: how many distinct system "
                    "types it carries. Genomes with no system are included in both."),
        table_html(["system", "tier", "systems", "genomes", "% of genomes"],
                   overview_rows, numeric_cols=(2, 3, 4))))

    key_table_rows = []
    for cat in CATEGORY_ORDER:
        for sy in [x for x in order if category_of(x) == cat]:
            k = key_for(sy)
            key_table_rows.append([sy, k["full"], cat, k["detected"], k["roles"],
                                   format(agg.systems_per_type[sy], ","),
                                   "%.2f" % (100.0 * genomes_carrying[sy] / n_screened)])

    sections.append("""
<section id="key">
  <h2>System key: what was detected, and what it might mean</h2>
  <p class="sub">Two different kinds of statement are kept apart throughout this report.
  <strong>Detected</strong> is a structural result: these components were found, they
  co-localised, and they satisfied the model quorum. <strong>Potential role</strong> is
  background biology from organisms where the machine has been studied &mdash; no effector,
  substrate, expression level or phenotype is measured here, so a detected system is a
  capability, not a demonstrated activity.</p>
  %s
  %s
  %s
  %s
</section>""" % (
        figure_html(figs["key"],
                    "Key to every model in the screen, grouped by what the machine does. "
                    "Swatches match the composition figure, so this doubles as an expanded legend."),
        figure_html(figs["categories"],
                    "Left: how the catalogue divides between the four functional tiers. "
                    "Right: the share of screened genomes carrying at least one system of each tier."),
        figure_html(figs["category_taxon"],
                    "Tier-level prevalence per taxon &mdash; the compact companion to the "
                    "per-model heatmap below."),
        table_html(["system", "full name", "tier", "what was detected", "potential roles",
                    "systems", "% genomes"], key_table_rows, numeric_cols=(5, 6))))

    spotlight_rows = [
        [r["organism"], r["system"], r["category"], r["group"],
         "%d / %s" % (r["carriers"], human(r["group_genomes"])),
         "%.2f" % r["prevalence"], "%.3f" % r["wholeness"], "%.2f" % r["score"],
         r["components"], r["loci"], r["accession"]] for r in notable]
    cross_rows = [[r["system"], r["organism"], r["domain"], "%.3f" % r["wholeness"],
                   "%.2f" % r["score"], r["components"], r["accession"]]
                  for r in crossover[:60]]
    rep_rows = [[r["organism"], r["n_types"], r["n_systems"], r["tiers"], r["systems"],
                 r["accession"]] for r in top_repertoires]
    n_spot = sum(1 for r in notable if r["spotlight"])

    cross_block = ("<h3>Systems found outside their expected domain</h3>" +
                   table_html(["system", "organism", "genome domain", "wholeness", "score",
                               "components", "accession"], cross_rows, numeric_cols=(3, 4, 5))
                   if cross_rows else
                   '<h3>Systems found outside their expected domain</h3>'
                   '<p class="sub">None passed the completeness threshold, which is the '
                   'expected result for a clean screen.</p>')

    sections.append("""
<section id="notable">
  <h2>Notable findings</h2>
  <p class="sub">Occurrences that are <strong>rare inside their own genus</strong> yet
  structurally convincing: at most %.0f%% of the genus carries the system, the genus has at
  least %d screened genomes, and the detected system is at least %.0f%% complete. %d of the %d
  candidates are effector-delivery systems (T3SS, T4SS, T6SS), the tier where an unexpected
  occurrence changes how a strain is interpreted. Treat this as a shortlist for manual and
  literature review: rarity in this dataset is not the same as absence from the literature, and
  a rare call can equally be a contaminated assembly, a misassigned taxon, or a mobile element
  captured in one strain.</p>
  %s
  %s
  <h3>Widest repertoires</h3>
  <p class="sub">Genomes carrying the largest number of distinct system types. These are the
  strains where the full functional range is available in one organism.</p>
  %s
  %s
  %s
</section>""" % (
        args.rare_max_carriers * 100, args.rare_min_genomes, args.rare_min_wholeness * 100,
        n_spot, len(notable),
        (figure_html(figs["notable"],
                     "Rare but complete occurrences, effector-delivery systems first. The "
                     "annotation gives how many genomes of that genus carry the system.") +
         table_html(["organism", "system", "tier", "genus", "carriers / genus", "% of genus",
                     "wholeness", "score", "components", "loci", "accession"],
                    spotlight_rows, numeric_cols=(5, 6, 7, 8, 9))
         if notable else
         '<p class="sub">No occurrence passed all three filters. Either the screened genera are '
         'too uniform for rarity to mean anything at these settings, or the thresholds are too '
         'strict &mdash; raise <code>--rare-max-carriers</code>, lower '
         '<code>--rare-min-genomes</code>, or relax <code>--rare-min-wholeness</code> and '
         'rebuild.</p>'),
        "",
        figure_html(figs["repertoire"],
                    "Genomes carrying the widest range of distinct system types, stacked by tier."),
        table_html(["organism", "system types", "systems", "tiers", "systems detected",
                    "accession"], rep_rows, numeric_cols=(1, 2)),
        cross_block))

    sections.append("""
<section id="taxonomy">
  <h2>Taxonomic distribution</h2>
  <p class="sub">The %d most heavily sequenced %s-level groups. Prevalence is the percentage of
  screened genomes in that group carrying at least one system of the type; composition is the
  share of all systems found in the group. Sequencing effort is deeply uneven across taxa, so
  read prevalence as a property of the sequenced sample rather than of the clade.</p>
  %s
  %s
</section>""" % (
        len(top_taxa), args.rank,
        figure_html(figs["heatmap"],
                    "Percentage of screened genomes in each group carrying at least one "
                    "system of the type. Group sizes are given on the axis."),
        figure_html(figs["composition"],
                    "Share of all systems detected in each group, by system type.")))

    sections.append("""
<section id="cooccurrence">
  <h2>Co-occurrence</h2>
  <p class="sub">Conditional prevalence: of the genomes carrying the row system, what percentage
  also carry the column system. The matrix is asymmetric on purpose &mdash; a rare system may sit
  almost always beside a common one without the reverse being true.</p>
  %s
</section>""" % figure_html(
        figs["cooccurrence"],
        "Conditional co-occurrence between the most abundant system types."))

    comp_blocks = "".join(
        figure_html(b64, "Component prevalence across the %s %s systems detected. A component "
                         "below 100%% is either accessory in the model or absent from that "
                         "particular locus." % (human(agg.systems_per_type[s]), s))
        for s, b64 in component_figs)

    hit_rows = [[k, format(v, ","), "%.2f" % (100.0 * v / max(1, sum(agg.hit_types.values())))]
                for k, v in agg.hit_types.most_common()]
    status_rows = [[k, format(v, ","),
                    "%.2f" % (100.0 * v / max(1, sum(agg.gene_status.values())))]
                   for k, v in agg.gene_status.most_common()]

    sections.append("""
<section id="quality">
  <h2>Quality and architecture</h2>
  <p class="sub">Model completeness (wholeness) is the fraction of the model's mandatory and
  accessory components that were found; the score is MacSyFinder's, rewarding mandatory
  components inside clusters and penalising redundancy. Together they separate confident calls
  from marginal ones.</p>
  %s
  %s
  %s
  <h3>Component prevalence</h3>
  %s
  <h3>Hit classification</h3>
  %s
  %s
</section>""" % (
        figure_html(figs["quality"],
                    "Distribution of wholeness, score and component count per system type. "
                    "Boxes are quartiles, whiskers the 1.5 IQR range, outliers omitted."),
        figure_html(figs["loci"],
                    "Systems detected in one locus versus spread over two or more, which the "
                    "models allow for the multi-locus system types."),
        figure_html(figs["proteome"],
                    "Systems per genome against proteome size. Larger genomes carry more "
                    "systems, so proteome size is a confounder in any cross-taxon comparison."),
        comp_blocks,
        table_html(["hit type", "hits", "% of hits"], hit_rows, numeric_cols=(1, 2)),
        table_html(["gene status", "hits", "% of hits"], status_rows, numeric_cols=(1, 2))))

    figure_note = ("Publication figures were written to <code>figures/</code> as %s at %d dpi; "
                   "text in the PDF and SVG stays editable text rather than outlines."
                   % (", ".join(f.strip().upper() for f in formats), args.dpi)) if formats else \
        "Figure files were not written (<code>--figure-formats none</code>)."

    sections.append("""
<section id="methods">
  <h2>Methods and caveats</h2>
  <p class="sub">Every protein set was searched with the TXSScan profile library using
  SecretionGenie, which reimplements MacSyFinder v2 selection: HMMER gathering cut-offs where a
  profile defines them, otherwise an E-value ceiling; an i-Evalue ceiling and a minimum profile
  coverage per hit; one best hit per protein; co-localisation within each model's
  <code>inter_gene_max_space</code>; then the model quorum over mandatory, accessory and
  forbidden components, with exchangeable components and multi-locus combinations resolved into
  a best solution per genome.</p>
  <p class="sub">%s Aggregate tables are in <code>tables/</code>, one TSV per figure, so every
  number here can be recomputed or replotted independently.</p>
  <h3>What this analysis does not show</h3>
  <p class="sub">Detection is structural. A predicted system says that the genes are present and
  arranged as the model requires; it does not show that the machine is expressed, assembled or
  functional, and it identifies no effectors or substrates. Prevalence reflects the sequenced
  sample, which is dominated by clinical and model organisms. Draft assemblies fragment loci and
  lose systems whose components fall on different contigs, so absence in a fragmented assembly is
  weak evidence. Labels inherit the blind spots of profile HMMs: components too divergent for the
  profiles are missed, and a quorum that just fails leaves no system at all.</p>
  <h3>References</h3>
  <ul class="refs">
    <li>Abby SS et al. (2016) Identification of protein secretion systems in bacterial genomes.
      <em>Sci Rep</em> 6:23080. https://doi.org/10.1038/srep23080</li>
    <li>Denise R, Abby SS, Rocha EPC (2019) Diversification of the type IV filament superfamily
      into machines for adhesion, protein secretion, DNA uptake and motility.
      <em>PLoS Biol</em> 17:e3000390. https://doi.org/10.1371/journal.pbio.3000390</li>
    <li>N&eacute;ron B et al. (2023) MacSyFinder v2: improved modelling and search engine to
      identify molecular systems in genomes. <em>Peer Community J</em> 3:e28.
      https://doi.org/10.24072/pcjournal.250</li>
    <li>Bongiovanni TR et al. (2024) Assembly of a unique membrane complex in type VI secretion
      systems of Bacteroidota. <em>Nat Commun</em> 15:429.
      https://doi.org/10.1038/s41467-023-44426-1</li>
    <li>TXSScan models: https://github.com/macsy-models/TXSScan</li>
  </ul>
</section>""" % figure_note)

    toc = ('<nav class="toc">'
           '<a href="#overview">Overview</a>'
           '<a href="#key">System key</a>'
           '<a href="#notable">Notable findings</a>'
           '<a href="#taxonomy">Taxonomic distribution</a>'
           '<a href="#cooccurrence">Co-occurrence</a>'
           '<a href="#quality">Quality and architecture</a>'
           '<a href="#methods">Methods</a>'
           '</nav>')

    subtitle = ("%s genomes &middot; %s predicted systems &middot; grouped by %s &middot; built %s"
                % (human(n_screened), human(agg.n_systems), html.escape(args.rank),
                   time.strftime("%Y-%m-%d %H:%M")))

    doc = ("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\"/>\n"
           "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"/>\n"
           "<title>" + html.escape(args.title) + "</title>\n<style>" + CSS + "</style>\n"
           "</head>\n<body>\n"
           "<header class=\"hero\"><h1>" + html.escape(args.title) + "</h1>"
           "<p>" + subtitle + "</p></header>\n" + toc + "\n<main>\n" +
           "\n".join(sections) +
           "\n<footer>Built by secretion_report.py from a SecretionGenie screen. "
           "Detected systems are structural predictions, not demonstrated activities."
           "</footer>\n</main>\n</body>\n</html>\n")

    if SINK.written:
        write_tsv(os.path.join(tables_dir, "figure_manifest.tsv"),
                  ["figure_name", "file"],
                  [[os.path.basename(pth).split("_", 1)[1].rsplit(".", 1)[0],
                    os.path.relpath(pth, args.outdir)] for pth in SINK.written])

    out_path = os.path.join(args.outdir, "secretion_report.html")
    with open(out_path, "w") as handle:
        handle.write(doc)

    sys.stderr.write("\nwrote %s (%.1f MB)\n" % (out_path, len(doc) / 1e6))
    sys.stderr.write("aggregate tables in %s\n" % tables_dir)
    if SINK.written:
        sys.stderr.write("%d publication figure files in %s (%s at %d dpi; vector text stays "
                         "editable)\n" % (len(SINK.written), figures_dir,
                                          "/".join(f.strip() for f in formats), args.dpi))
    return 0


if __name__ == "__main__":
    sys.exit(main())
