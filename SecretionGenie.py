#!/usr/bin/env python3
"""
SecretionGenie - part of the MagicLamp suite.

HMM-based identification and categorization of bacterial protein secretion
systems and evolutionarily related surface appendages (T1SS-T9SS, T4aP/T4bP,
Tad, MSH, ComM, Flagellum, Archaeal-T4P), using the TXSScan HMM profile
library and MacSyFinder model definitions:

    https://github.com/macsy-models/TXSScan

SecretionGenie re-implements, in a single self-contained script, the search,
filtering, co-localization ("gene neighborhood / island / cluster") and
quorum logic of MacSyFinder v2 / MacSyLib (Neron et al. 2023,
Peer Community Journal; Abby et al. 2016, Sci Rep; Denise et al. 2019,
PLoS Biol), so that gene clusters are called with the same rules TXSScan
itself uses:

  * hmmsearch is run per profile with --cut_ga when the profile carries a GA
    (gathering) bit-score threshold, otherwise with -E <e_value_search>
    (MacSyFinder default behaviour; --no_cut_ga forces -E everywhere).
  * domain hits are filtered on independent E-value (i_evalue_sel, def 0.001)
    and on profile coverage (coverage_profile, def 0.5), coverage being
    (hmm_to - hmm_from + 1) / profile_length.
  * at most one hit is kept per protein (best bit score across all profiles).
  * hits are clustered per replicon/contig: two hits co-localize when the
    number of intervening non-matched genes is <= inter_gene_max_space
    (the gene-level value when defined, else the model-level value; when both
    genes define one, the minimum is taken). Circular topology is supported.
  * a group of co-localized hits is only a cluster if it encodes more than one
    gene, is not exclusively 'neutral', or is a single hit of a 'loner' gene
    (true loner) or of a model with min_genes_required == 1.
  * clusters (and, for multi_loci models, every combination of clusters) are
    challenged against the model quorum: min_mandatory_genes_required,
    min_genes_required, and absence of 'forbidden' genes. 'exchangeables'
    count for the gene they replace.
  * candidate systems are scored exactly like MacSyFinder: mandatory 1.0,
    accessory 0.5, neutral 0.0, x0.8 for exchangeables, x0.7 out of cluster
    (true loner / multi-system hit), -1.5 per redundant function across loci.
  * 'multi_system' genes can rescue rejected candidates, and the reported set
    of systems is the compatible combination (no shared hits, except
    multi_model / multi_system genes) with the maximum total score.

Two CSV outputs are produced, in the same style as ATPGenie:

  * <out>/secretiongenie-summary.csv   : one block of rows per detected system
                                         (gene-level detail), blocks separated
                                         by a '####' line. Every row carries the
                                         gene's strand, nucleotide start and end,
                                         gene length in nt and protein length in
                                         aa, taken from the Prodigal FASTA headers
                                         or, with --gbk, from the CDS features of
                                         the GenBank file. Where the protein file
                                         has no coordinates, the four positional
                                         fields are NA and protein length is still
                                         reported.
  * <out>/secretiongenie.heatmap.csv   : heatmap-compatible matrix, rows are
                                         system types (and, optionally with
                                         --genes, individual components),
                                         columns are genomes/bins

Developed by Arkadiy Garber. Please send comments and inquiries to
ark@midauthorbio.com
"""

from collections import defaultdict
import xml.etree.ElementTree as ET
import itertools
import argparse
import textwrap
import shutil
import glob
import sys
import os
import re


# =====================================================================
# Helper functions (mirrored from ATPGenie.py / Lucifer.py for
# stylistic consistency with the rest of the MagicLamp suite)
# =====================================================================

def lastItem(ls):
    x = ''
    for i in ls:
        x = i
    return x


def RemoveDuplicates(ls):
    empLS = []
    for i in ls:
        if i not in empLS:
            empLS.append(i)
    return empLS


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length - 1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x) - 1]


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
    outString = "".join(emptyList)
    return outString


def fasta(fasta_file):
    seq = ''
    header = ''
    Dict = defaultdict(lambda: 'EMPTY')
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:].split(" ")[0]
                seq = ''
            else:
                header = i[1:].split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    return Dict


def prodigal_coords(faa_path):
    """ORF -> (start, end, strand) read from Prodigal-style FASTA headers.

    Prodigal writes '>orf_1 # 190 # 255 # 1 # ID=1_1;partial=00;...', so the
    nucleotide coordinates and the strand of every called gene are already in
    the protein file. Headers without those fields (a plain '>orf_1', as written
    by most GenBank-to-FASTA converters) simply yield no entry, and the caller
    reports NA for that ORF rather than guessing.
    """
    coords = {}
    try:
        handle = open(faa_path)
    except (IOError, OSError):
        return coords
    with handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            fields = [f.strip() for f in line[1:].rstrip("\n").split(" # ")]
            if len(fields) < 4:
                continue
            orf = fields[0].split(" ")[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
                strand = int(fields[3])
            except ValueError:
                continue
            coords[orf] = (min(start, end), max(start, end), "-" if strand < 0 else "+")
    return coords


def genbank_cds_coords(gbk_path):
    """locus_tag (and protein_id) -> (start, end, strand) for every CDS feature.

    Locations are read from the feature table, so complement() gives the strand
    and join()/order() collapse to the outermost bounds of the coding region.
    Partial markers ('<1..>500') are tolerated.
    """
    coords = {}
    in_cds = False
    location = ""
    current = None
    try:
        handle = open(gbk_path, errors="replace")
    except (IOError, OSError):
        return coords

    def finish(loc):
        numbers = re.findall(r"\d+", loc)
        if len(numbers) < 2:
            return None
        values = [int(n) for n in numbers]
        return (min(values), max(values), "-" if "complement" in loc else "+")

    with handle:
        for line in handle:
            raw = line.rstrip("\n")
            feature = re.match(r"^ {5}(\S+)\s+(.*)$", raw)
            if feature:
                in_cds = feature.group(1) == "CDS"
                location = feature.group(2).strip() if in_cds else ""
                current = None
                continue
            if not in_cds:
                continue
            stripped = raw.strip()
            if stripped.startswith("/"):
                if current is None:
                    current = finish(location)
                tag = re.match(r'^/(locus_tag|protein_id|old_locus_tag)="([^"]+)"', stripped)
                if tag and current:
                    coords.setdefault(tag.group(2), current)
            elif stripped and current is None:
                # continuation of a multi-line location
                location += stripped
    return coords


def coord_fields(coords, orf, seq):
    """the five extra output columns: strand, start, end, gene length, protein length"""
    protein = str(seq).replace("*", "") if seq and seq != "EMPTY" else ""
    prot_len = str(len(protein)) if protein else "NA"
    entry = coords.get(orf) if coords else None
    if entry:
        start, end, strand = entry
        return [strand, str(start), str(end), str(end - start + 1), prot_len]
    return ["NA", "NA", "NA", "NA", prot_len]


def filt(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls.append(string)
    ls = filt(ls, ["", "\n"])
    return ls


def csvSafe(value):
    """the summary file is a plain CSV: strip anything that would break it"""
    s = str(value)
    s = s.replace(",", ";").replace("\n", " ").replace("\r", " ")
    return s


# =====================================================================
# Model objects (MacSyFinder grammar, version 2.x)
# =====================================================================

MANDATORY = "mandatory"
ACCESSORY = "accessory"
NEUTRAL = "neutral"
FORBIDDEN = "forbidden"


def as_bool(value, default=False):
    if value is None:
        return default
    return str(value).strip().lower() in ("1", "true", "yes", "on")


class Gene(object):
    """a gene of a model; 'ref' is set for exchangeable (analog/homolog) genes"""

    def __init__(self, name, status, loner=False, multi_system=False,
                 multi_model=False, inter_gene_max_space=None, ref=None):
        self.name = name
        self.status = status
        self.loner = loner
        self.multi_system = multi_system
        self.multi_model = multi_model
        self.inter_gene_max_space = inter_gene_max_space
        self.ref = ref                 # the gene this one can replace
        self.exchangeables = []

    @property
    def is_exchangeable(self):
        return self.ref is not None

    def function(self):
        """the biological function: own name, or the name of the gene replaced"""
        return self.ref.name if self.ref is not None else self.name


class Model(object):

    def __init__(self, name, fqn, path):
        self.name = name
        self.fqn = fqn
        self.path = path
        self.inter_gene_max_space = None
        self.min_mandatory_genes_required = None
        self.min_genes_required = None
        self.max_nb_genes = None
        self.multi_loci = False
        self.genes = []                # top-level genes only

    # ---------------------------------------------------------------
    def all_genes(self):
        """top level genes + their exchangeables"""
        out = []
        for gene in self.genes:
            out.append(gene)
            out.extend(gene.exchangeables)
        return out

    def gene_by_name(self):
        return {g.name: g for g in self.all_genes()}

    def genes_of_status(self, status):
        return [g for g in self.genes if g.status == status]

    @property
    def mandatory_genes(self):
        return self.genes_of_status(MANDATORY)

    @property
    def accessory_genes(self):
        return self.genes_of_status(ACCESSORY)

    @property
    def neutral_genes(self):
        return self.genes_of_status(NEUTRAL)

    @property
    def forbidden_genes(self):
        return self.genes_of_status(FORBIDDEN)


def parse_model(xml_path, definitions_root):
    """parse one MacSyFinder model definition file"""
    tree = ET.parse(xml_path)
    root = tree.getroot()
    if root.tag != "model":
        raise ValueError("%s is not a MacSyFinder model definition" % xml_path)

    name = os.path.basename(xml_path)[:-len(".xml")]
    rel = os.path.relpath(xml_path, definitions_root)
    fqn = os.path.splitext(rel)[0].replace(os.sep, "/")
    model = Model(name, fqn, xml_path)

    model.inter_gene_max_space = int(root.get("inter_gene_max_space"))
    model.min_mandatory_genes_required = root.get("min_mandatory_genes_required")
    model.min_genes_required = root.get("min_genes_required")
    model.max_nb_genes = root.get("max_nb_genes")
    model.multi_loci = as_bool(root.get("multi_loci"), False)

    for gene_node in root.findall("gene"):
        igms = gene_node.get("inter_gene_max_space")
        gene = Gene(gene_node.get("name"),
                    gene_node.get("presence"),
                    loner=as_bool(gene_node.get("loner")),
                    multi_system=as_bool(gene_node.get("multi_system")),
                    multi_model=as_bool(gene_node.get("multi_model")),
                    inter_gene_max_space=int(igms) if igms is not None else None)
        # exchangeables (MacSyFinder >= 2 grammar), also accept the v1 tags
        for holder_tag in ("exchangeables", "homologs", "analogs"):
            for holder in gene_node.findall(holder_tag):
                for ex_node in holder.findall("gene"):
                    ex_igms = ex_node.get("inter_gene_max_space")
                    ex = Gene(ex_node.get("name"),
                              gene.status,
                              loner=as_bool(ex_node.get("loner"), gene.loner),
                              multi_system=as_bool(ex_node.get("multi_system"), gene.multi_system),
                              multi_model=as_bool(ex_node.get("multi_model"), gene.multi_model),
                              inter_gene_max_space=int(ex_igms) if ex_igms is not None
                              else gene.inter_gene_max_space,
                              ref=gene)
                    gene.exchangeables.append(ex)
        model.genes.append(gene)

    # defaults, as in MacSyFinder
    nb_mandatory = len(model.mandatory_genes)
    nb_full = len(model.mandatory_genes) + len(model.accessory_genes)
    model.min_mandatory_genes_required = int(model.min_mandatory_genes_required) \
        if model.min_mandatory_genes_required is not None else nb_mandatory
    model.min_genes_required = int(model.min_genes_required) \
        if model.min_genes_required is not None else nb_mandatory
    model.max_nb_genes = int(model.max_nb_genes) \
        if model.max_nb_genes is not None else nb_full
    if model.max_nb_genes < 1:
        model.max_nb_genes = 1
    return model


# =====================================================================
# Hits
# =====================================================================

class CoreHit(object):
    """one protein / one profile match that passed the HMMER filters"""

    __slots__ = ("orf", "gene_name", "replicon", "position", "seq_len",
                 "i_eval", "score", "cov_profile", "cov_seq", "begin", "end")

    def __init__(self, orf, gene_name, replicon, position, seq_len,
                 i_eval, score, cov_profile, cov_seq, begin, end):
        self.orf = orf
        self.gene_name = gene_name
        self.replicon = replicon
        self.position = position
        self.seq_len = seq_len
        self.i_eval = i_eval
        self.score = score
        self.cov_profile = cov_profile
        self.cov_seq = cov_seq
        self.begin = begin
        self.end = end


class ModelHit(object):
    """a CoreHit bound to a gene of one particular model"""

    __slots__ = ("core", "gene", "status", "loner", "multi_system",
                 "multi_model", "counterpart")

    def __init__(self, core, gene):
        self.core = core
        self.gene = gene
        self.status = gene.status
        # 'loner'/'multi_system' below are the *hit* properties: they become
        # True only once the hit is recognized as being out of a regular
        # cluster (true loner) or imported from another cluster/system.
        self.loner = False
        self.multi_system = False
        self.multi_model = gene.multi_model
        self.counterpart = []

    @property
    def position(self):
        return self.core.position

    @property
    def score(self):
        return self.core.score

    @property
    def orf(self):
        return self.core.orf

    def function(self):
        return self.gene.function()

    def hit_type(self):
        if self.loner and self.multi_system:
            return "loner_multi_system"
        if self.loner:
            return "loner"
        if self.multi_system:
            return "multi_system"
        return "cluster"


class HitWeight(object):
    """MacSyFinder scoring weights"""

    def __init__(self, mandatory=1.0, accessory=0.5, neutral=0.0, itself=1.0,
                 exchangeable=0.8, out_of_cluster=0.7, redundancy_penalty=1.5):
        self.mandatory = mandatory
        self.accessory = accessory
        self.neutral = neutral
        self.itself = itself
        self.exchangeable = exchangeable
        self.out_of_cluster = out_of_cluster
        self.redundancy_penalty = redundancy_penalty


def get_best_hits(hits):
    """keep, for each protein, only the best-scoring hit (MacSyFinder's
       get_best_hits with key='score')"""
    register = defaultdict(list)
    for hit in hits:
        register[(hit.replicon, hit.position)].append(hit)
    best = []
    for group in register.values():
        group.sort(key=lambda h: (-h.score, h.i_eval, h.gene_name))
        best.append(group[0])
    return best


# =====================================================================
# Clusters (gene neighborhoods / genomic islands / loci)
# =====================================================================

class Cluster(object):

    _counter = itertools.count(1)

    def __init__(self, hits, model, weights):
        self.hits = sorted(hits, key=lambda h: h.position)
        self.model = model
        self.weights = weights
        self.id = "c%d" % next(Cluster._counter)
        self._score = None

    def __len__(self):
        return len(self.hits)

    @property
    def replicon(self):
        return self.hits[0].core.replicon

    @property
    def loner(self):
        """True when made of hits of one single gene and that gene is a loner"""
        return len({h.gene.name for h in self.hits}) == 1 and self.hits[0].gene.loner

    @property
    def multi_system(self):
        return len(self.hits) == 1 and self.hits[0].gene.multi_system

    @property
    def is_out_of_cluster(self):
        """the cluster is a lone hit imported in the system (loner / multi system)"""
        return any(h.loner or h.multi_system for h in self.hits) and len(self.hits) == 1

    def functions(self):
        return frozenset(h.function() for h in self.hits)

    def fulfilled_function(self, *names):
        wanted = set()
        for n in names:
            wanted.add(n.name if isinstance(n, Gene) else n)
        return self.functions() & wanted

    @property
    def score(self):
        if self._score is not None:
            return self._score
        w = self.weights
        out_of_cluster = self.is_out_of_cluster or self.loner or self.multi_system
        seen = {}
        for hit in self.hits:
            if hit.status == MANDATORY:
                hit_score = w.mandatory
            elif hit.status == ACCESSORY:
                hit_score = w.accessory
            elif hit.status == NEUTRAL:
                hit_score = w.neutral
            else:
                # forbidden genes never contribute to the score; in ordered mode
                # their presence rejects the candidate upstream, in unordered mode
                # they are only reported
                hit_score = 0.0
            hit_score *= w.exchangeable if hit.gene.is_exchangeable else w.itself
            if out_of_cluster:
                hit_score *= w.out_of_cluster
            func = hit.function()
            # only one occurrence of each function counts, the best one
            if func not in seen or hit_score > seen[func]:
                seen[func] = hit_score
        self._score = sum(seen.values())
        return self._score

    def positions(self):
        return [h.position for h in self.hits]


def colocates(h1, h2, model, replicon_len, topology, forced_igms=None):
    """MacSyFinder's _colocates: number of intervening genes <= inter_gene_max_space"""
    dist = h2.position - h1.position - 1
    if forced_igms is not None:
        igms = forced_igms
    else:
        d1 = h1.gene.inter_gene_max_space
        d2 = h2.gene.inter_gene_max_space
        if d1 is None and d2 is None:
            igms = model.inter_gene_max_space
        elif d1 is None:
            igms = d2
        elif d2 is None:
            igms = d1
        else:
            igms = min(d1, d2)
    if 0 <= dist <= igms:
        return True
    if dist <= 0 and topology == "circular":
        # h1 and h2 overlap the origin of replication
        dist = replicon_len - h1.position + h2.position - 1
        return dist <= igms
    return False


def scaffold_to_cluster(scaffold, model, weights):
    """a co-localized group of hits becomes a Cluster only if it encodes more
       than one gene and is not exclusively neutral; a single-gene group is
       kept when the gene is a loner, or when min_genes_required == 1"""
    if not scaffold:
        return None
    gene_types = {h.gene.name for h in scaffold}
    if len(gene_types) > 1:
        if all(h.status == NEUTRAL for h in scaffold):
            return None
        return Cluster(scaffold, model, weights)
    cluster = Cluster(scaffold, model, weights)
    if cluster.loner:
        return cluster           # squashed later into a true loner
    if model.min_genes_required == 1:
        if scaffold[0].status == NEUTRAL:
            return None
        return cluster
    return None


def clusterize_hits(hits, model, weights, replicon_len, topology, forced_igms=None):
    """build the co-localized groups of hits of one replicon for one model"""
    clusters = []
    if not hits:
        return clusters
    hits = sorted(hits, key=lambda h: (h.position, -h.score))
    # one hit per protein (the best one)
    hits = [next(grp) for _, grp in itertools.groupby(hits, key=lambda h: h.position)]

    scaffold = [hits[0]]
    previous = hits[0]
    for hit in hits[1:]:
        if colocates(previous, hit, model, replicon_len, topology, forced_igms):
            scaffold.append(hit)
        else:
            cluster = scaffold_to_cluster(scaffold, model, weights)
            if cluster is not None:
                clusters.append(cluster)
            scaffold = [hit]
        previous = hit

    cluster = scaffold_to_cluster(scaffold, model, weights)
    if cluster is not None:
        clusters.append(cluster)
    elif topology == "circular":
        # the trailing scaffold may join the first cluster across the origin
        if clusters and colocates(scaffold[-1], clusters[0].hits[0], model,
                                  replicon_len, topology, forced_igms):
            merged = Cluster(scaffold + clusters[0].hits, model, weights)
            clusters[0] = merged
        elif colocates(scaffold[-1], hits[0], model, replicon_len, topology, forced_igms):
            cluster = scaffold_to_cluster(scaffold + [hits[0]], model, weights)
            if cluster is not None:
                clusters.append(cluster)
    if topology == "circular" and len(clusters) > 1:
        if colocates(clusters[-1].hits[-1], clusters[0].hits[0], model,
                     replicon_len, topology, forced_igms):
            clusters[0] = Cluster(clusters[-1].hits + clusters[0].hits, model, weights)
            clusters = clusters[:-1]
    return clusters


def get_true_loners(clusters, model, weights):
    """split the clusters into regular clusters and 'true loners' (a cluster made
       of hits of a single loner gene: the gene is allowed to sit alone)"""
    true_clusters = []
    loner_hits = defaultdict(list)
    for cluster in clusters:
        if cluster.loner:
            loner_hits[cluster.hits[0].function()].extend(cluster.hits)
        else:
            true_clusters.append(cluster)

    true_loners = {}
    for func, hits in loner_hits.items():
        # tag the hits as loners (and loner+multi_system when relevant),
        # then keep only the best-scoring one for that function
        for hit in hits:
            hit.loner = True
            if hit.gene.multi_system:
                hit.multi_system = True
        best = sorted(hits, key=lambda h: (-h.score, h.core.i_eval, h.position))[0]
        best.counterpart = [h for h in hits if h is not best]
        true_loners[func] = Cluster([best], model, weights)
    return true_clusters, true_loners


def combine_clusters(clusters, true_loners, multi_loci=False, max_combinations=200000,
                     max_loci=None):
    """generate the combinations of clusters (+ loners) to be challenged against
       the model quorum, as MacSyFinder's combine_clusters does"""
    if not clusters:
        combinations = []
    elif multi_loci:
        top = len(clusters)
        if max_loci is not None:
            top = min(top, max_loci)
        combinations = []
        for i in range(1, top + 1):
            for comb in itertools.combinations(clusters, i):
                combinations.append(comb)
                if len(combinations) >= max_combinations:
                    break
            if len(combinations) >= max_combinations:
                sys.stderr.write("  [warning] too many loci for %s: cluster combinations "
                                 "truncated at %d\n" % (clusters[0].model.name, max_combinations))
                break
    else:
        combinations = [(clst,) for clst in clusters]

    loner_items = list(true_loners.items())
    loner_combinations = []
    for i in range(1, len(loner_items) + 1):
        loner_combinations.extend(itertools.combinations(loner_items, i))

    with_loners = []
    for loner_comb in loner_combinations:
        loner_functions = [item[0] for item in loner_comb]
        loners = [item[1] for item in loner_comb]
        if combinations:
            for one in combinations:
                to_add = True
                for clst in one:
                    if clst.fulfilled_function(*loner_functions):
                        to_add = False
                        break
                if to_add:
                    with_loners.append(tuple(list(one) + loners))
        # a definition may be fulfilled by loners only
        with_loners.append(tuple(loners))
    return list(combinations) + with_loners


# =====================================================================
# Systems and quorum
# =====================================================================

class System(object):

    _counter = itertools.count(1)

    def __init__(self, model, clusters, weights, genome):
        self.model = model
        self.clusters = list(clusters)
        self.weights = weights
        self.genome = genome
        self.id = "%s_%s_%d" % (os.path.basename(genome), model.name, next(System._counter))
        self._score = None

    @property
    def replicon(self):
        return self.clusters[0].replicon

    @property
    def hits(self):
        out = []
        for cluster in self.clusters:
            out.extend(cluster.hits)
        return sorted(out, key=lambda h: h.position)

    def core_hits(self):
        return {h.core for h in self.hits}

    def functions(self):
        out = set()
        for cluster in self.clusters:
            out |= cluster.functions()
        return out

    def fulfilled_function(self, *names):
        out = set()
        for cluster in self.clusters:
            out |= cluster.fulfilled_function(*names)
        return out

    @property
    def position(self):
        pos = [h.position for h in self.hits]
        return min(pos), max(pos)

    @property
    def loci_nb(self):
        return len([c for c in self.clusters if not (c.loner or c.multi_system)])

    @property
    def wholeness(self):
        present = self.functions()
        wanted = [g.name for g in self.model.mandatory_genes + self.model.accessory_genes]
        return sum(1 for n in wanted if n in present) / float(self.model.max_nb_genes)

    @property
    def score(self):
        """MacSyFinder's system score: sum of the regular cluster scores, minus
           the redundancy penalty, plus the out-of-cluster loner/multi-system
           contributions whose function is not already fulfilled"""
        if self._score is not None:
            return self._score
        regular = []
        out_of_cluster = []
        for cluster in self.clusters:
            if cluster.loner or cluster.multi_system:
                out_of_cluster.append(cluster)
            else:
                regular.append(cluster)

        score = sum(c.score for c in regular)
        for gene in self.model.mandatory_genes + self.model.accessory_genes:
            nb = sum(1 for c in regular if c.fulfilled_function(gene.name))
            if nb:
                score -= (nb - 1) * self.weights.redundancy_penalty

        regular_functions = set()
        for cluster in regular:
            regular_functions |= cluster.functions()
        for cluster in out_of_cluster:
            func = cluster.hits[0].function()
            if func not in regular_functions:
                score += cluster.score
        self._score = score
        return score

    def is_compatible(self, other):
        """two systems are compatible when they do not share hits; hits of
           multi_system genes may be shared between occurrences of the same
           model, hits of multi_model genes between different models"""
        if self.model is other.model:
            mine = {h.core for h in self.hits if not (h.loner or h.multi_system)}
            theirs = {h.core for h in other.hits if not (h.loner or h.multi_system)}
            return not (mine & theirs)
        mine = {h.core: h for h in self.hits}
        theirs = {h.core: h for h in other.hits}
        for core in set(mine) & set(theirs):
            if not (mine[core].multi_model and theirs[core].multi_model):
                return False
        return True

    def multisystem_hits(self):
        return {h for h in self.hits if h.gene.multi_system}


class RejectedCandidate(object):

    def __init__(self, model, clusters, reasons):
        self.model = model
        self.clusters = list(clusters)
        self.reasons = reasons

    @property
    def hits(self):
        out = []
        for cluster in self.clusters:
            out.extend(cluster.hits)
        return sorted(out, key=lambda h: h.position)

    def fulfilled_function(self, *names):
        out = set()
        for cluster in self.clusters:
            out |= cluster.fulfilled_function(*names)
        return out


def match_quorum(model, clusters, weights, genome, ignore_forbidden=False):
    """MacSyFinder's OrderedMatchMaker.match: check min_mandatory_genes_required,
       min_genes_required and the absence of forbidden genes.
       In 'unordered' mode (ignore_forbidden=True) the forbidden genes are only
       reported, not used to reject the candidate, exactly as MacSyFinder's
       UnorderedMatchMaker does."""
    mandatory = set()
    accessory = set()
    neutral = set()
    forbidden = set()
    forbidden_hits = []
    for cluster in clusters:
        for hit in cluster.hits:
            func = hit.function()
            if hit.status == MANDATORY:
                mandatory.add(func)
            elif hit.status == ACCESSORY:
                accessory.add(func)
            elif hit.status == NEUTRAL:
                neutral.add(func)
            elif hit.status == FORBIDDEN:
                forbidden.add(func)
                forbidden_hits.append(hit)

    reasons = []
    if forbidden and not ignore_forbidden:
        reasons.append("%d forbidden gene occurrence(s): %s"
                       % (len(forbidden_hits),
                          "; ".join(h.gene.name for h in forbidden_hits)))
    if len(mandatory) < model.min_mandatory_genes_required:
        reasons.append("the quorum of mandatory genes required (%d) is not reached: %d"
                       % (model.min_mandatory_genes_required, len(mandatory)))
    if len(accessory) + len(mandatory) < model.min_genes_required:
        reasons.append("the quorum of genes required (%d) is not reached: %d"
                       % (model.min_genes_required, len(accessory) + len(mandatory)))

    if reasons:
        return RejectedCandidate(model, clusters, reasons)
    return System(model, clusters, weights, genome)


def combine_multisystems(rejected, multi_system_clusters):
    """try to rescue rejected candidates with multi_system hits found elsewhere"""
    new_combinations = []
    ms_combinations = []
    for i in range(1, len(multi_system_clusters) + 1):
        ms_combinations.extend(itertools.combinations(multi_system_clusters, i))
    for candidate in rejected:
        for ms_comb in ms_combinations:
            functions = [c.hits[0].function() for c in ms_comb]
            if not candidate.fulfilled_function(*functions):
                new_combinations.append(tuple(candidate.clusters + list(ms_comb)))
    return new_combinations


# =====================================================================
# Best solution: maximum-score set of mutually compatible systems
# =====================================================================

def find_best_solution(systems, max_nodes=60):
    """among candidate systems, keep the combination that shares no hit and
       maximizes the sum of the scores (MacSyFinder's find_best_solutions).
       Exact search (Bron-Kerbosch on the compatibility graph) for small
       instances, greedy otherwise."""
    if not systems:
        return []
    systems = sorted(systems, key=lambda s: (-s.score, s.position[0], s.model.name))
    if len(systems) > max_nodes:
        return _greedy_solution(systems)

    index = {s: i for i, s in enumerate(systems)}
    neighbours = {i: set() for i in range(len(systems))}
    for si, sj in itertools.combinations(systems, 2):
        if si.is_compatible(sj):
            neighbours[index[si]].add(index[sj])
            neighbours[index[sj]].add(index[si])

    best = {"score": None, "clique": []}
    node_budget = [400000]

    def expand(r, p, x):
        if node_budget[0] <= 0:
            return
        node_budget[0] -= 1
        if not p and not x:
            score = sum(systems[i].score for i in r)
            if best["score"] is None or score > best["score"]:
                best["score"] = score
                best["clique"] = list(r)
            return
        pivot = max(p | x, key=lambda v: len(neighbours[v]))
        for v in list(p - neighbours[pivot]):
            expand(r | {v}, p & neighbours[v], x & neighbours[v])
            p = p - {v}
            x = x | {v}

    expand(set(), set(range(len(systems))), set())
    if node_budget[0] <= 0 and not best["clique"]:
        return _greedy_solution(systems)
    solution = [systems[i] for i in best["clique"]]
    solution.sort(key=lambda s: (s.replicon, s.position[0], s.model.name))
    return solution


def _greedy_solution(systems):
    kept = []
    for candidate in sorted(systems, key=lambda s: (-s.score, s.position[0], s.model.name)):
        if all(candidate.is_compatible(k) for k in kept):
            kept.append(candidate)
    kept.sort(key=lambda s: (s.replicon, s.position[0], s.model.name))
    return kept


# =====================================================================
# HMM library handling
# =====================================================================

def profile_info(hmm_path):
    """return (profile length, has GA threshold)"""
    length = None
    has_ga = False
    with open(hmm_path) as handle:
        for line in handle:
            if line.startswith("LENG"):
                length = int(line.split()[1])
            elif line.startswith("GA "):
                has_ga = True
            elif line.startswith("HMM "):
                break
    return length, has_ga


def parse_domtblout(path, gene_name, profile_len, i_evalue_sel, coverage_profile,
                    position_of, replicon_of, length_of):
    """parse hmmsearch --domtblout and apply MacSyFinder's hit filters"""
    hits = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 22:
                continue
            orf = fields[0]
            try:
                seq_len = int(fields[2])
                qlen = int(fields[5])
                i_eval = float(fields[12])
                score = float(fields[13])
                hmm_from = int(fields[15])
                hmm_to = int(fields[16])
                ali_from = int(fields[17])
                ali_to = int(fields[18])
            except ValueError:
                continue
            if orf not in position_of:
                continue
            if i_eval > i_evalue_sel:
                continue
            plen = profile_len or qlen
            cov_profile = (hmm_to - hmm_from + 1) / float(plen)
            if cov_profile < coverage_profile:
                continue
            cov_seq = (ali_to - ali_from + 1) / float(length_of.get(orf, seq_len) or seq_len)
            hits.append(CoreHit(orf, gene_name, replicon_of[orf], position_of[orf],
                                length_of.get(orf, seq_len), i_eval, score,
                                cov_profile, cov_seq, ali_from, ali_to))
    return hits


def hmm_format_tag(hmm_path):
    """the HMMER save-format version of a profile ('HMMER3/b', 'HMMER3/f', ...).
       hmmsearch cannot read a file that mixes formats, so profiles are grouped
       by format when the concatenated library is built."""
    try:
        with open(hmm_path) as handle:
            first = handle.readline().strip()
    except OSError:
        return "unknown"
    tag = first.split("[")[0].strip() if first else "unknown"
    return tag if tag else "unknown"


def build_hmm_db(needed_genes, profile_meta, db_dir):
    """Concatenate the profiles into a handful of HMM libraries so that a couple
       of hmmsearch calls per genome replace one call per profile (~4x faster on
       a bacterial genome). Results are unchanged because:

         * the NAME line of every model is rewritten to the profile file name,
           so the 'query name' column of the domtblout is the TXSScan gene name;
         * profiles are grouped by threshold class, so those carrying a GA
           bit-score cut-off are searched with --cut_ga and the others with -E,
           exactly as MacSyFinder does per profile;
         * profiles are additionally grouped by HMMER save format, because
           hmmsearch refuses a library that mixes HMMER3/b and HMMER3/f models.

       :return: list of (path, threshold_kind) with threshold_kind in
                ('cut_ga', 'evalue')
    """
    os.makedirs(db_dir, exist_ok=True)
    manifest_path = os.path.join(db_dir, "txss_hmm_db.manifest")

    groups = defaultdict(list)
    for gene_name in sorted(needed_genes):
        kind = "cut_ga" if profile_meta[gene_name][1] else "evalue"
        fmt = hmm_format_tag(needed_genes[gene_name]).replace("/", "_").replace(" ", "")
        groups[(kind, fmt)].append(gene_name)

    plan = []
    manifest_lines = []
    for (kind, fmt), genes in sorted(groups.items()):
        path = os.path.join(db_dir, "txss_%s_%s.hmm" % (kind, fmt))
        plan.append((path, kind, genes))
        manifest_lines.extend("%s\t%s\t%s" % (kind, fmt, g) for g in genes)
    manifest = "\n".join(manifest_lines) + "\n"

    reusable = os.path.isfile(manifest_path)
    if reusable:
        try:
            reusable = open(manifest_path).read() == manifest and \
                all(os.path.isfile(path) for path, _, _ in plan)
        except OSError:
            reusable = False
    if reusable:
        return [(path, kind) for path, kind, _ in plan]

    for path, _, genes in plan:
        tmp = path + ".tmp"
        with open(tmp, "w") as out:
            for gene_name in genes:
                renamed = False
                last = "\n"
                with open(needed_genes[gene_name]) as handle:
                    for line in handle:
                        if not renamed and line.startswith("NAME"):
                            out.write("NAME  %s\n" % gene_name)
                            renamed = True
                        else:
                            out.write(line)
                        last = line
                if not last.endswith("\n"):
                    # several TXSScan profiles have no trailing newline: without
                    # this, concatenation glues '//' onto the next file's header
                    # and hmmsearch aborts with 'bad file format'
                    out.write("\n")
        os.replace(tmp, path)

    tmp = manifest_path + ".tmp"
    with open(tmp, "w") as out:
        out.write(manifest)
    os.replace(tmp, manifest_path)
    return [(path, kind) for path, kind, _ in plan]


def parse_domtblout_multi(path, needed_genes, i_evalue_sel, coverage_profile,
                          position_of, replicon_of, length_of):
    """parse a domtblout produced by searching a concatenated HMM library:
       the query name column carries the TXSScan gene name and the qlen column
       the profile length. Same filters as parse_domtblout()."""
    hits = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 22:
                continue
            orf = fields[0]
            gene_name = fields[3]
            if orf not in position_of or gene_name not in needed_genes:
                continue
            try:
                seq_len = int(fields[2])
                qlen = int(fields[5])
                i_eval = float(fields[12])
                score = float(fields[13])
                hmm_from = int(fields[15])
                hmm_to = int(fields[16])
                ali_from = int(fields[17])
                ali_to = int(fields[18])
            except ValueError:
                continue
            if i_eval > i_evalue_sel:
                continue
            cov_profile = (hmm_to - hmm_from + 1) / float(qlen)
            if cov_profile < coverage_profile:
                continue
            cov_seq = (ali_to - ali_from + 1) / float(length_of.get(orf, seq_len) or seq_len)
            hits.append(CoreHit(orf, gene_name, replicon_of[orf], position_of[orf],
                                length_of.get(orf, seq_len), i_eval, score,
                                cov_profile, cov_seq, ali_from, ali_to))
    return hits


# =====================================================================
# main
# =====================================================================

def main():

    parser = argparse.ArgumentParser(
        prog="MagicLamp.py SecretionGenie",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(r'''
        *******************************************************

        SecretionGenie - part of the MagicLamp suite.
        Developed by Arkadiy Garber.
        Please send comments and inquiries to ark@midauthorbio.com

          ____                    _   _              ____            _
         / ___|  ___  ___ _ __ ___| |_(_) ___  _ __  / ___| ___ _ __ (_) ___
         \___ \ / _ \/ __| '__/ _ \ __| |/ _ \| '_ \| |  _ / _ \ '_ \| |/ _ \
          ___) |  __/ (__| | |  __/ |_| | (_) | | | | |_| |  __/ | | | |  __/
         |____/ \___|\___|_|  \___|\__|_|\___/|_| |_|\____|\___|_| |_|_|\___|

                 protein secretion systems & related appendages
                    (TXSScan profiles + MacSyFinder logic)

        Two-stage analysis:
          1. Profile genomes/proteomes with the TXSScan HMM library
             (hmmsearch --cut_ga when the profile carries a GA threshold,
             otherwise -E; hits filtered on i-evalue and profile coverage).
          2. Reconstruct gene neighborhoods/islands (co-localization with
             inter_gene_max_space, loners, multi-loci) and validate them
             against the TXSScan model quorum (min_mandatory_genes_required,
             min_genes_required, forbidden genes), exactly as MacSyFinder
             version 2 does, then keep the best-scoring compatible solution.

        *******************************************************
        '''))

    parser.add_argument('-bin_dir', type=str, help="directory of bins/genomes/assemblies", default="NA")

    parser.add_argument('-bin_ext', type=str, help="extension for bins (do not include the period)", default="NA")

    parser.add_argument('-txss_dir', type=str, default="NA",
                        help="path to the TXSScan package (the directory containing 'definitions/' and "
                             "'profiles/'). If not provided, SecretionGenie looks at the $txss_models and "
                             "$txss_hmms environment variables, then next to MagicLamp.py "
                             "(hmms/txss), then in the current directory (TXSScan/).")

    parser.add_argument('-models', type=str, default="all",
                        help="comma-separated list of systems to search for, e.g. "
                             "T2SS,T3SS,T6SSi (default=all). Names may also be given fully "
                             "qualified, e.g. bacteria/diderm/T6SSi")

    parser.add_argument('-d', type=int, default=0,
                        help="maximum distance (number of intervening genes) between two genes to be "
                             "considered part of the same genomic cluster. By default (0) the "
                             "inter_gene_max_space defined by each TXSScan model is used, which is the "
                             "behaviour of MacSyFinder/TXSScan. Provide an integer to override it for "
                             "every model.")

    parser.add_argument('-out', type=str, help="name output directory (default=secretiongenie_out)",
                        default="secretiongenie_out")

    parser.add_argument('-t', type=int, help="number of threads to use for HMMSEARCH (default=1)", default=1)

    parser.add_argument('-e_value_search', type=float, default=0.1,
                        help="maximum e-value for hmmsearch reporting, used for profiles that do not "
                             "carry a GA bit-score threshold (default=0.1, as in MacSyFinder)")

    parser.add_argument('-i_evalue_sel', type=float, default=0.001,
                        help="maximum independent e-value for a domain hit to be selected "
                             "(default=0.001, as in MacSyFinder)")

    parser.add_argument('-coverage_profile', type=float, default=0.5,
                        help="minimum fraction of the HMM profile that must be covered by the "
                             "alignment (default=0.5, as in MacSyFinder)")

    parser.add_argument('-topology', type=str, default="linear",
                        help="'linear' or 'circular'. Contigs of draft genomes/bins are best treated as "
                             "linear (default); use 'circular' for closed replicons.")

    parser.add_argument('-redundancy_penalty', type=float, default=1.5,
                        help="penalty applied when a function is encoded by several loci of the same "
                             "system (default=1.5)")

    parser.add_argument('-max_loci', type=int, default=6,
                        help="for multi_loci models, maximum number of loci combined into one candidate "
                             "system (default=6). Guards against combinatorial explosion on very "
                             "fragmented assemblies.")

    parser.add_argument('--no_cut_ga', type=str, const=True, nargs="?",
                        help="include this flag to ignore the GA bit-score thresholds carried by the "
                             "profiles and use -e_value_search for every profile")

    parser.add_argument('--unordered', type=str, const=True, nargs="?",
                        help="include this flag to skip the co-localization step and call systems on "
                             "the gene quorum only (MacSyFinder's 'unordered' db-type). Useful for "
                             "highly fragmented metagenomic assemblies.")

    parser.add_argument('--gbk', type=str, help="include this flag if your bins are in Genbank format",
                        const=True, nargs="?")

    parser.add_argument('--meta', type=str,
                        help="include this flag if the provided contigs are from metagenomic/"
                             "metatranscriptomic assemblies", const=True, nargs="?")

    parser.add_argument('--genes', type=str, const=True, nargs="?",
                        help="include this flag to add one heatmap row per system component "
                             "(system|gene) in addition to the per-system rows")

    parser.add_argument('--norm', type=str, const=True, nargs="?",
                        help="include this flag if you would like the counts for each secretion system "
                             "category to be normalized to the number of predicted ORFs in each genome "
                             "or metagenome. Without normalization, SecretionGenie will create a "
                             "heatmap-compatible CSV output with raw counts. With normalization, "
                             "SecretionGenie will create a heatmap-compatible CSV with 'normalized "
                             "abundances'.")

    parser.add_argument('--keep_hmm', type=str, const=True, nargs="?",
                        help="include this flag to keep the raw hmmsearch output files")

    parser.add_argument('--fast_hmm', type=str, const=True, nargs="?",
                        help="include this flag to search a single concatenated HMM library per "
                             "genome instead of calling hmmsearch once per profile. Results are "
                             "identical (each profile keeps its own GA threshold) but it is about "
                             "4x faster, which matters when screening thousands of genomes.")

    parser.add_argument('-hmm_db_dir', type=str, default="NA",
                        help="where to build/find the concatenated HMM library used by --fast_hmm "
                             "(default: <out>/hmm_db). Point several parallel SecretionGenie "
                             "processes at one shared directory to build it only once.")

    parser.add_argument('--build_hmm_db', type=str, const=True, nargs="?",
                        help="build the concatenated HMM library for --fast_hmm and exit. Run this "
                             "once before launching parallel jobs that share -hmm_db_dir.")

    args = parser.parse_known_args()[0]

    # -----------------------------------------------------------------
    # Locate the TXSScan library (mirrors ATPGenie's discovery logic)
    # -----------------------------------------------------------------

    def have_lib(root):
        if not root:
            return False
        return os.path.isdir(os.path.join(root, "definitions")) and \
            os.path.isdir(os.path.join(root, "profiles"))

    candidates = []
    if args.txss_dir != "NA":
        candidates.append(args.txss_dir)
    for var in ("txss_models", "txss_hmms", "TXSSCAN"):
        if os.environ.get(var):
            candidates.append(os.environ[var])
            candidates.append(allButTheLast(os.environ[var], "/"))
    magiclamp = shutil.which("MagicLamp.py")
    if magiclamp:
        candidates.append(os.path.join(os.path.dirname(magiclamp), "hmms", "txss"))
    candidates.append(os.path.join(os.getcwd(), "TXSScan"))
    candidates.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "TXSScan"))

    TXSSdir = ""
    for candidate in candidates:
        if have_lib(candidate):
            TXSSdir = candidate
            break

    if not TXSSdir:
        print("SecretionGenie could not locate the TXSScan library.")
        print("Please clone it (git clone https://github.com/macsy-models/TXSScan.git) and provide the")
        print("path with -txss_dir, or export it as $txss_models. The directory must contain the two")
        print("subdirectories: definitions/ and profiles/")
        raise SystemExit

    definitionsDir = os.path.join(TXSSdir, "definitions")
    profilesDir = os.path.join(TXSSdir, "profiles")

    # -----------------------------------------------------------------
    # Argument validation
    # -----------------------------------------------------------------

    print("checking arguments")

    if args.bin_dir != "NA":
        binDir = args.bin_dir + "/"
        binDirLS = os.listdir(args.bin_dir)
        print(".")
    else:
        print("Looks like you did not provide a directory of genomes/bins or assemblies.")
        print("Exiting")
        raise SystemExit

    if args.bin_ext != "NA":
        print(".")
    else:
        print("Looks like you did not provide an extension for your genomes/bins or assemblies, so "
              "SecretionGenie does not know which files in the provided directory are FASTA files that "
              "you would like analyzed.")
        print("Exiting")
        raise SystemExit

    if lastItem(args.out) == "/":
        outDirectory = "%s" % args.out[0:len(args.out) - 1]
    else:
        outDirectory = "%s" % args.out

    try:
        os.listdir(outDirectory)
        print("Looks like you already have a directory with the name: " + outDirectory)
        print("Ok, proceeding with analysis!")
    except FileNotFoundError:
        print(".")
        os.system("mkdir -p %s" % outDirectory)
    os.system("mkdir -p %s/ORF_calls" % outDirectory)
    os.system("mkdir -p %s/HMM_results" % outDirectory)

    if not shutil.which("hmmsearch"):
        print("hmmsearch was not found in your $PATH. Please install HMMER3.")
        raise SystemExit
    if not args.gbk and not shutil.which("prodigal"):
        print("prodigal was not found in your $PATH. Please install Prodigal (or use --gbk).")
        raise SystemExit

    topology = args.topology.strip().lower()
    if topology not in ("linear", "circular"):
        print("-topology must be either 'linear' or 'circular'")
        raise SystemExit

    weights = HitWeight(redundancy_penalty=args.redundancy_penalty)
    forced_igms = args.d if args.d and args.d > 0 else None

    # -----------------------------------------------------------------
    # Load the models
    # -----------------------------------------------------------------

    allModels = []
    for xml in sorted(glob.glob(os.path.join(definitionsDir, "**", "*.xml"), recursive=True)):
        try:
            allModels.append(parse_model(xml, definitionsDir))
        except Exception as err:
            print("could not parse model %s (%s)" % (xml, err))

    if args.models.strip().lower() in ("all", "na"):
        models = allModels
    else:
        wanted = [m.strip() for m in args.models.split(",") if m.strip()]
        models = []
        for name in wanted:
            hit = [m for m in allModels if m.name == name or m.fqn == name]
            if not hit:
                print("model '%s' is not part of the TXSScan package; available models are:" % name)
                print("  " + ", ".join(sorted(m.name for m in allModels)))
                raise SystemExit
            models.extend(hit)
    models = sorted(RemoveDuplicates(models), key=lambda m: m.name)

    print("%d TXSScan model(s) will be searched: %s" % (len(models), ", ".join(m.name for m in models)))

    # the non-redundant set of profiles needed by these models
    neededGenes = {}
    for model in models:
        for gene in model.all_genes():
            neededGenes[gene.name] = os.path.join(profilesDir, gene.name + ".hmm")
    missing = [g for g, p in neededGenes.items() if not os.path.isfile(p)]
    if missing:
        print("the following HMM profiles are missing from %s:" % profilesDir)
        print("  " + ", ".join(sorted(missing)))
        raise SystemExit

    profileMeta = {}
    for gene_name, path in neededGenes.items():
        profileMeta[gene_name] = profile_info(path)

    print("%d HMM profiles will be used" % len(neededGenes))

    hmmDbDir = args.hmm_db_dir if args.hmm_db_dir != "NA" else os.path.join(outDirectory, "hmm_db")
    hmmDbs = []
    if args.build_hmm_db or args.fast_hmm:
        hmmDbs = build_hmm_db(neededGenes, profileMeta, hmmDbDir)
    if args.build_hmm_db:
        print("concatenated HMM library built in %s" % hmmDbDir)
        for path, kind in hmmDbs:
            print("  %s (%s)" % (path, kind))
        print("re-run with --fast_hmm -hmm_db_dir %s to use it" % hmmDbDir)
        raise SystemExit
    if args.fast_hmm:
        print("fast mode: %d hmmsearch call(s) per genome against %s"
              % (len(hmmDbs), hmmDbDir))

    print("All required arguments provided!")
    print("")

    # -----------------------------------------------------------------
    # ORF calling (Prodigal) or GenBank protein extraction
    # -----------------------------------------------------------------

    BinDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    coordsByGenome = {}     # genome -> {ORF: (start, end, strand)}
    genomes = []
    for i in sorted(binDirLS):
        if lastItem(i.split(".")) != args.bin_ext:
            continue
        genomes.append(i)
        cell = i
        if not args.gbk:
            try:
                testFile = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
                print("ORFs for %s found. Skipping Prodigal, and going with %s-proteins.faa" % (i, i))
                for line in testFile:
                    if re.match(r'>', line):
                        if re.findall(r'\|]', line):
                            print("Looks like one of your fasta files has a header containing the "
                                  "character: \\|")
                            print("Unfortunately, this is a problem for SecretionGenie because it uses "
                                  "that character as delimiter to store important information.")
                            print("Please rename your FASTA file headers")
                            raise SystemExit
                testFile.close()
            except FileNotFoundError:
                binFile = open("%s/%s" % (binDir, i), "r")
                for line in binFile:
                    if re.match(r'>', line):
                        if re.findall(r'\|]', line):
                            print("Looks like one of your fasta files has a header containing the "
                                  "character: \\|")
                            print("Unfortunately, this is a problem for SecretionGenie because it uses "
                                  "that character as delimiter to store important information.")
                            print("Please rename your FASTA file headers")
                            raise SystemExit
                binFile.close()

                print("Finding ORFs for " + cell)
                if args.meta:
                    os.system("prodigal -i %s/%s -a %s/ORF_calls/%s-proteins.faa "
                              "-o %s/ORF_calls/%s-prodigal.out -p meta -q"
                              % (binDir, i, outDirectory, i, outDirectory, i))
                else:
                    os.system("prodigal -i %s/%s -a %s/ORF_calls/%s-proteins.faa "
                              "-o %s/ORF_calls/%s-prodigal.out -q"
                              % (binDir, i, outDirectory, i, outDirectory, i))
        else:
            os.system('gtt-genbank-to-AA-seqs -i %s/%s -o %s/%s.faa' % (binDir, i, outDirectory, i))
            faa = fasta(open("%s/%s.faa" % (outDirectory, i)))

            # nucleotide coordinates and strand come from the CDS features of the
            # source GenBank file, keyed by locus_tag / protein_id
            gbkCoords = genbank_cds_coords("%s/%s" % (binDir, i))
            altCoords = {}
            renamedCoords = {}

            gbkDict = defaultdict(list)
            counter = 0
            count = 0
            for gbkline in open("%s/%s" % (binDir, i)):
                if re.findall(r'/locus_tag', gbkline.rstrip()):
                    count += 1

            if count > 0:
                locus = "unknown"
                seenTags = defaultdict(set)
                for gbkline in open("%s/%s" % (binDir, i)):
                    ls = gbkline.rstrip()
                    if re.findall(r'LOCUS', ls):
                        locus = ls.split("       ")[1].split(" ")[0]
                    if re.findall(r'/locus_tag', ls):
                        locusTag = remove(ls.split("=")[1], ["\""])
                        counter += 1
                    if counter > 0:
                        # a locus_tag appears on both the 'gene' and the 'CDS'
                        # feature of the same gene, so keep only its first
                        # occurrence: duplicated ORFs would double every rank and
                        # halve the apparent distance between neighbouring genes
                        if locusTag not in seenTags[locus]:
                            seenTags[locus].add(locusTag)
                            gbkDict[locus].append(locusTag)
                        counter = 0
            else:
                locus = "unknown"
                seenAlt = defaultdict(set)
                for gbkline in open("%s/%s" % (binDir, i)):
                    ls = gbkline.rstrip()
                    if re.findall(r'LOCUS', ls):
                        locus = ls.split("       ")[1].split(" ")[0]
                    if re.findall(r'gene   ', ls):
                        gene = ls.split("            ")[1]
                        start = remove(gene.split("..")[0], ["c", "o", "m", "p", "l", "e", "m",
                                                             "e", "n", "t", "(", ")"])
                        end = remove(gene.split("..")[1], ["c", "o", "m", "p", "l", "e", "m",
                                                           "e", "n", "t", "(", ")"])
                        altContigName = (locus + "_" + start + "_" + end)
                        try:
                            altCoords[altContigName] = (
                                min(int(start), int(end)), max(int(start), int(end)),
                                "-" if re.findall(r'complement', ls) else "+")
                        except ValueError:
                            pass
                        counter += 1
                    if counter > 0:
                        if altContigName not in seenAlt[locus]:
                            seenAlt[locus].add(altContigName)
                            gbkDict[locus].append(altContigName)
                        counter = 0

            idxOut = open("%s/ORF_calls/%s-proteins.idx" % (outDirectory, i), "w")
            faaOut = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "w")
            for gbkkey1 in gbkDict.keys():
                counter = 0
                for gbkey2 in gbkDict[gbkkey1]:
                    counter += 1
                    if len(faa[gbkey2]) > 0 and faa[gbkey2] != "EMPTY":
                        newOrf = gbkkey1 + "_" + str(counter)
                        idxOut.write(gbkey2 + "," + newOrf + "\n")
                        entry = gbkCoords.get(gbkey2) or altCoords.get(gbkey2)
                        if entry:
                            renamedCoords[newOrf] = entry
                        faaOut.write(">" + newOrf + "\n")
                        faaOut.write(str(faa[gbkey2]) + "\n")
            idxOut.close()
            faaOut.close()

        proteins = fasta(open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i)))
        for key in proteins.keys():
            if key:
                BinDict[cell][key.split(" # ")[0]] = proteins[key]

        if args.gbk:
            coordsByGenome[cell] = renamedCoords
        else:
            # Prodigal keeps start / end / strand in the protein FASTA headers,
            # whether it was just run here or the ORF calls were supplied
            coordsByGenome[cell] = prodigal_coords(
                "%s/ORF_calls/%s-proteins.faa" % (outDirectory, i))
        located = len(coordsByGenome[cell])
        if not located:
            print("  note: no gene coordinates found for %s - strand, start, end and "
                  "gene_length_nt will be NA in the output" % i)

    if not genomes:
        print("No file with the extension '%s' was found in %s" % (args.bin_ext, args.bin_dir))
        raise SystemExit

    # -----------------------------------------------------------------
    # PASS 1: profile every genome with the TXSScan HMM library
    # -----------------------------------------------------------------

    print("")
    print("starting main pipeline...")
    print("Pass 1: profiling secretion system components with the TXSScan HMM library")

    hitsByGenome = {}        # genome -> list of CoreHit (best hit per protein)
    replicon_len = {}        # genome -> {replicon: number of ORFs}
    orfCount = {}            # genome -> number of ORFs

    for genome in genomes:
        faaPath = "%s/ORF_calls/%s-proteins.faa" % (outDirectory, genome)
        hmmDirOut = "%s/HMM_results/%s-HMM" % (outDirectory, genome)
        os.system("mkdir -p %s" % hmmDirOut)

        # ORF -> replicon (contig) and rank of the ORF on that contig
        position_of = {}
        replicon_of = {}
        length_of = {}
        counters = defaultdict(int)
        for orf in BinDict[genome]:
            contig = allButTheLast(orf, "_")
            replicon_of[orf] = contig if contig else orf
            length_of[orf] = len(str(BinDict[genome][orf]).replace("*", ""))
        # ORFs are numbered by Prodigal in genomic order along each contig
        def orf_rank(orf):
            tail = lastItem(orf.split("_"))
            try:
                return int(tail)
            except ValueError:
                return None
        for orf in sorted(BinDict[genome], key=lambda o: (replicon_of[o], orf_rank(o) if
                                                          orf_rank(o) is not None else 0, o)):
            rank = orf_rank(orf)
            if rank is None:
                counters[replicon_of[orf]] += 1
                rank = counters[replicon_of[orf]]
            position_of[orf] = rank
        orfCount[genome] = len(position_of)
        lengths = defaultdict(int)
        for orf, pos in position_of.items():
            if pos > lengths[replicon_of[orf]]:
                lengths[replicon_of[orf]] = pos
        replicon_len[genome] = lengths

        allHits = []

        if args.fast_hmm:
            sys.stdout.write("analyzing " + genome + " (single-pass HMM search)   \r")
            sys.stdout.flush()
            for db_path, kind in hmmDbs:
                threshold = "--cut_ga" if kind == "cut_ga" else "-E %f" % args.e_value_search
                domtbl = "%s/%s.domtblout" % (hmmDirOut, os.path.basename(db_path))
                os.system("hmmsearch --cpu %d %s --domtblout %s -o /dev/null %s %s"
                          % (int(args.t), threshold, domtbl, db_path, faaPath))
                if not os.path.isfile(domtbl):
                    continue
                allHits.extend(parse_domtblout_multi(domtbl, neededGenes,
                                                     args.i_evalue_sel, args.coverage_profile,
                                                     position_of, replicon_of, length_of))
                if not args.keep_hmm:
                    os.system("rm -f %s" % domtbl)
            print("")
            hitsByGenome[genome] = get_best_hits(allHits)
            continue

        count = 0
        total = len(neededGenes)
        for gene_name in sorted(neededGenes):
            count += 1
            perc = (count / float(total)) * 100
            sys.stdout.write("analyzing " + genome + ": %d%%   \r" % perc)
            sys.stdout.flush()

            hmmPath = neededGenes[gene_name]
            profile_len, has_ga = profileMeta[gene_name]
            if args.no_cut_ga or not has_ga:
                threshold = "-E %f" % args.e_value_search
            else:
                threshold = "--cut_ga"
            domtbl = "%s/%s.domtblout" % (hmmDirOut, gene_name)
            os.system("hmmsearch --cpu %d %s --domtblout %s -o /dev/null %s %s"
                      % (int(args.t), threshold, domtbl, hmmPath, faaPath))
            if not os.path.isfile(domtbl):
                continue
            allHits.extend(parse_domtblout(domtbl, gene_name, profile_len,
                                           args.i_evalue_sel, args.coverage_profile,
                                           position_of, replicon_of, length_of))
            if not args.keep_hmm:
                os.system("rm -f %s" % domtbl)
        print("")

        # keep only the best-scoring hit per protein, as MacSyFinder does
        hitsByGenome[genome] = get_best_hits(allHits)

    # a table of all retained hits, for reference
    allHitsOut = open("%s/secretiongenie-all_hits.csv" % outDirectory, "w")
    allHitsOut.write("file,ORF,replicon,position,strand,gene_start,gene_end,gene_length_nt,"
                     "protein_length_aa,gene,i_evalue,bit_score,"
                     "profile_coverage,sequence_coverage,begin,end\n")
    for genome in genomes:
        for hit in sorted(hitsByGenome[genome], key=lambda h: (h.replicon, h.position)):
            allHitsOut.write(",".join([csvSafe(genome), csvSafe(hit.orf), csvSafe(hit.replicon),
                                       str(hit.position)] +
                                      coord_fields(coordsByGenome.get(genome), hit.orf,
                                                   BinDict[genome][hit.orf]) +
                                      [csvSafe(hit.gene_name),
                                       str(hit.i_eval), str(hit.score),
                                       "%.3f" % hit.cov_profile, "%.3f" % hit.cov_seq,
                                       str(hit.begin), str(hit.end)]) + "\n")
    allHitsOut.close()

    # -----------------------------------------------------------------
    # PASS 2: gene neighborhoods / clusters and model quorum
    # -----------------------------------------------------------------

    print("Pass 2: identifying gene neighborhoods (clusters/islands) and validating model quorum")

    systemsByGenome = defaultdict(list)
    rejectedByGenome = defaultdict(list)

    for genome in genomes:
        print("." + " " + genome)
        coreHits = hitsByGenome[genome]
        hitsByReplicon = defaultdict(list)
        for hit in coreHits:
            hitsByReplicon[hit.replicon].append(hit)

        candidates = []

        if args.unordered:
            # MacSyFinder's 'unordered' db-type: no notion of genetic distance, so
            # the clustering step is skipped and a single occurrence of each model
            # is called for the whole dataset, on the gene quorum only
            kept = []
            for model in models:
                geneMap = model.gene_by_name()
                modelHits = [ModelHit(h, geneMap[h.gene_name])
                             for h in coreHits if h.gene_name in geneMap]
                if not modelHits:
                    continue
                pseudo = Cluster(modelHits, model, weights)
                res = match_quorum(model, [pseudo], weights, genome, ignore_forbidden=True)
                if isinstance(res, System):
                    kept.append(res)
                else:
                    rejectedByGenome[genome].append(res)
            kept.sort(key=lambda s: s.model.name)
            systemsByGenome[genome] = kept
            print("   %d system(s) with a genetic potential detected (unordered mode)" % len(kept))
            continue

        for replicon in sorted(hitsByReplicon):
            for model in models:
                geneMap = model.gene_by_name()
                modelHits = [ModelHit(h, geneMap[h.gene_name])
                             for h in hitsByReplicon[replicon] if h.gene_name in geneMap]
                if not modelHits:
                    continue

                rep_len = replicon_len[genome].get(replicon, 0)
                clusters = clusterize_hits(modelHits, model, weights, rep_len,
                                           topology, forced_igms)
                true_clusters, true_loners = get_true_loners(clusters, model, weights)
                combinations = combine_clusters(true_clusters, true_loners,
                                                multi_loci=model.multi_loci,
                                                max_loci=args.max_loci)

                modelSystems = []
                modelRejected = []
                for combination in combinations:
                    if not combination:
                        continue
                    res = match_quorum(model, combination, weights, genome)
                    if isinstance(res, System):
                        modelSystems.append(res)
                    else:
                        modelRejected.append(res)

                # multi_system rescue: hits of multi_system genes belonging to an
                # accepted system may complete an otherwise rejected candidate
                msHits = set()
                for system in modelSystems:
                    msHits |= system.multisystem_hits()
                if msHits and modelRejected:
                    bestByFunction = {}
                    for hit in msHits:
                        func = hit.function()
                        if func not in bestByFunction or hit.score > bestByFunction[func].score:
                            bestByFunction[func] = hit
                    msClusters = []
                    for hit in bestByFunction.values():
                        clone = ModelHit(hit.core, hit.gene)
                        clone.multi_system = True
                        clone.loner = hit.loner
                        msClusters.append(Cluster([clone], model, weights))
                    for combination in combine_multisystems(modelRejected, msClusters):
                        res = match_quorum(model, combination, weights, genome)
                        if isinstance(res, System):
                            modelSystems.append(res)
                        else:
                            modelRejected.append(res)

                candidates.extend(modelSystems)
                rejectedByGenome[genome].extend(modelRejected)

        # keep the best-scoring set of mutually compatible systems, per replicon
        byReplicon = defaultdict(list)
        for system in candidates:
            byReplicon[system.replicon].append(system)
        kept = []
        for replicon in sorted(byReplicon):
            kept.extend(find_best_solution(byReplicon[replicon]))
        kept.sort(key=lambda s: (s.replicon, s.position[0], s.model.name))
        systemsByGenome[genome] = kept
        print("   %d system(s) detected" % len(kept))

    # -----------------------------------------------------------------
    # Summary CSV (gene-level detail, one block per detected system)
    # -----------------------------------------------------------------

    print("..")
    print("...")

    summaryPath = "%s/secretiongenie-summary.csv" % outDirectory
    out = open(summaryPath, "w")
    out.write("file,ORF,gene,function,gene_status,system,system_id,replicon,position,"
              "strand,gene_start,gene_end,gene_length_nt,protein_length_aa,locus,"
              "hit_type,evalue,bit_score,profile_coverage,sequence_coverage,system_score,"
              "system_wholeness,nb_loci,seq\n")
    for genome in genomes:
        for system in systemsByGenome[genome]:
            locusOf = {}
            regular = [c for c in system.clusters if not (c.loner or c.multi_system)]
            regular.sort(key=lambda c: c.hits[0].position)
            for num, cluster in enumerate(regular, start=1):
                for hit in cluster.hits:
                    locusOf[hit] = str(num)
            for cluster in system.clusters:
                if cluster in regular:
                    continue
                for hit in cluster.hits:
                    locusOf[hit] = "-"
            for hit in system.hits:
                seq = BinDict[genome][hit.orf]
                if seq == "EMPTY":
                    seq = ""
                out.write(",".join([csvSafe(genome), csvSafe(hit.orf), csvSafe(hit.gene.name),
                                    csvSafe(hit.function()), csvSafe(hit.status),
                                    csvSafe(system.model.name), csvSafe(system.id),
                                    csvSafe(hit.core.replicon), str(hit.position)] +
                                   coord_fields(coordsByGenome.get(genome), hit.orf, seq) +
                                   ["-" if args.unordered else locusOf.get(hit, "-"),
                                    "unordered" if args.unordered else hit.hit_type(),
                                    str(hit.core.i_eval), str(hit.score),
                                    "%.3f" % hit.core.cov_profile,
                                    "%.3f" % hit.core.cov_seq,
                                    "%.2f" % system.score,
                                    "%.3f" % system.wholeness,
                                    str(system.loci_nb),
                                    csvSafe(str(seq).replace("*", ""))]) + "\n")
            out.write("####################################################\n")
    out.close()

    # rejected candidates, for troubleshooting (as MacSyFinder reports them)
    rejPath = "%s/secretiongenie-rejected_candidates.csv" % outDirectory
    rej = open(rejPath, "w")
    rej.write("file,system,replicon,nb_clusters,genes,positions,reasons\n")
    for genome in genomes:
        seen = set()
        for candidate in rejectedByGenome[genome]:
            hits = candidate.hits
            if not hits:
                continue
            key = (candidate.model.name, tuple(h.position for h in hits))
            if key in seen:
                continue
            seen.add(key)
            rej.write(",".join([csvSafe(genome), csvSafe(candidate.model.name),
                                csvSafe(hits[0].core.replicon), str(len(candidate.clusters)),
                                csvSafe(";".join(h.gene.name for h in hits)),
                                csvSafe(";".join(str(h.position) for h in hits)),
                                csvSafe(" | ".join(candidate.reasons))]) + "\n")
    rej.close()

    if not args.keep_hmm:
        os.system("rm -f %s/ORF_calls/*-prodigal.out" % outDirectory)

    # -----------------------------------------------------------------
    # Heatmap-compatible CSV
    #   Rows: system types (optionally also system|component rows)
    #   Columns: genomes / bins
    # -----------------------------------------------------------------

    print("....")
    print(".....")

    cats = [model.name for model in models]
    if args.genes:
        for model in models:
            for gene in model.genes:
                if gene.status in (MANDATORY, ACCESSORY, NEUTRAL):
                    cats.append("%s|%s" % (model.name, gene.name))

    Dict = defaultdict(lambda: defaultdict(int))
    for genome in genomes:
        for system in systemsByGenome[genome]:
            Dict[genome][system.model.name] += 1
            if args.genes:
                for func in system.functions():
                    Dict[genome]["%s|%s" % (system.model.name, func)] += 1

    outHeat = open("%s/secretiongenie.heatmap.csv" % outDirectory, "w")
    outHeat.write("X,")
    for genome in genomes:
        outHeat.write(csvSafe(genome) + ",")
    outHeat.write("\n")
    for cat in cats:
        outHeat.write(csvSafe(cat) + ",")
        for genome in genomes:
            if args.norm:
                if orfCount.get(genome, 0) > 0:
                    outHeat.write(str((Dict[genome][cat] / float(orfCount[genome])) * float(100)) + ",")
                else:
                    outHeat.write("0,")
            else:
                outHeat.write(str(Dict[genome][cat]) + ",")
        outHeat.write("\n")
    outHeat.close()

    if not args.keep_hmm:
        os.system("rm -rf %s/HMM_results" % outDirectory)

    print("......")
    print(".......")
    print("Finished!")
    print("")
    print("Results are written to %s/secretiongenie-summary.csv and %s/secretiongenie.heatmap.csv"
          % (outDirectory, outDirectory))
    print("(all retained HMM hits: %s/secretiongenie-all_hits.csv; rejected gene clusters: "
          "%s/secretiongenie-rejected_candidates.csv)" % (outDirectory, outDirectory))
    print("Pipeline finished without crashing!!! Thanks for using :)")


if __name__ == '__main__':
    main()
