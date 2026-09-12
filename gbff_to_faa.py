#!/usr/bin/env python3
"""
gbff_to_faa.py - GenBank flat file (.gbff / .gb, optionally gzipped) to a
protein FASTA formatted exactly the way SecretionGenie.py / MagicLamp expect
when ORF calls are supplied instead of being predicted with Prodigal.

Headers are written as   >{replicon}_{n}
where {replicon} is the VERSION accession of the GenBank record (falling back
to the LOCUS name) and {n} is the rank of the CDS along that record, starting
at 1. That is exactly the ORF naming Prodigal produces, and it is what
SecretionGenie uses to derive the replicon and the gene position for the
co-localization ("gene neighborhood") step, so pre-supplying these files lets
you skip Prodigal entirely on RefSeq (GCF_) assemblies that already carry
annotations.

Only stdlib is used (no Biopython), so it runs in any environment.

Usage:
    gbff_to_faa.py input.gb.gz -o out/ORF_calls/ACC.gb-proteins.faa \
                   [--index out/ORF_calls/ACC.gb-proteins.idx] [--min-cds 1]

Exit codes:
    0  success
    3  no CDS translations found in the file (unannotated / RNA-only record)
    4  unreadable / truncated input
"""

import argparse
import gzip
import os
import re
import sys


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


FIRST_INT = re.compile(r"\d+")


def location_start(location):
    """lowest coordinate of a GenBank location string (handles complement/join/<>)"""
    match = FIRST_INT.search(location)
    return int(match.group(0)) if match else 0


def location_end(location):
    ints = FIRST_INT.findall(location)
    return int(ints[-1]) if ints else 0


def parse_records(handle):
    """yield (replicon_name, [cds, ...]) per GenBank record.
       cds = dict(start, end, strand, translation, locus_tag, protein_id, product)"""
    locus = None
    version = None
    cds_list = []
    in_features = False

    current = None            # CDS being accumulated
    qual_name = None          # qualifier currently being continued
    qual_chunks = []

    def flush_qualifier():
        if current is None or qual_name is None:
            return
        value = "".join(qual_chunks).strip().strip('"')
        if qual_name == "translation":
            current["translation"] = re.sub(r"[^A-Za-z*]", "", value)
        elif qual_name in ("locus_tag", "protein_id", "product", "gene"):
            current[qual_name] = value.replace(",", ";")

    def flush_cds():
        flush_qualifier()
        if current is not None and current.get("translation"):
            cds_list.append(current)

    for line in handle:
        line = line.rstrip("\n")

        if line.startswith("LOCUS"):
            locus = (line.split()[1] if len(line.split()) > 1 else None)
            version = None
            in_features = False
            continue

        if line.startswith("VERSION"):
            parts = line.split()
            if len(parts) > 1:
                version = parts[1]
            continue

        if line.startswith("FEATURES"):
            in_features = True
            current = None
            qual_name = None
            qual_chunks = []
            continue

        if line.startswith("ORIGIN") or line.startswith("CONTIG"):
            in_features = False
            flush_cds()
            current = None
            qual_name = None
            qual_chunks = []
            continue

        if line.startswith("//"):
            flush_cds()
            name = version or locus
            if name:
                yield name, cds_list
            locus = version = None
            cds_list = []
            current = None
            qual_name = None
            qual_chunks = []
            in_features = False
            continue

        if not in_features:
            continue

        # feature key lines start at column 6 (5 spaces), qualifiers at column 22
        if len(line) > 5 and line[5] != " " and not line.startswith("      "):
            # a new feature begins: close the previous CDS
            flush_cds()
            current = None
            qual_name = None
            qual_chunks = []
            fields = line.split(None, 1)
            if fields and fields[0] == "CDS":
                location = fields[1].strip() if len(fields) > 1 else ""
                current = {"location": location,
                           "start": location_start(location),
                           "end": location_end(location),
                           "strand": "-" if location.startswith("complement") else "+",
                           "translation": "", "locus_tag": "", "protein_id": "",
                           "product": "", "gene": ""}
            continue

        if current is None:
            continue

        stripped = line.strip()
        if stripped.startswith("/"):
            flush_qualifier()
            if "=" in stripped:
                key, _, value = stripped[1:].partition("=")
                qual_name = key
                qual_chunks = [value]
            else:
                qual_name = stripped[1:]
                qual_chunks = []
        elif qual_name is not None:
            # continuation line of the current qualifier
            qual_chunks.append(stripped if qual_name != "translation" else stripped)
        else:
            # continuation of a multi-line location
            current["location"] += stripped
            current["start"] = location_start(current["location"])
            current["end"] = location_end(current["location"])


def main():
    parser = argparse.ArgumentParser(description="GenBank flat file -> SecretionGenie-ready protein FASTA")
    parser.add_argument("gbff", help="input .gb / .gbff, optionally .gz")
    parser.add_argument("-o", "--out", required=True, help="output protein FASTA")
    parser.add_argument("--index", default=None,
                        help="optional TSV mapping ORF id -> locus_tag/protein_id/product/coords")
    parser.add_argument("--min-cds", type=int, default=1,
                        help="fail (exit 3) if fewer than this many translated CDS are found (default 1)")
    args = parser.parse_args()

    tmp_out = args.out + ".tmp"
    tmp_idx = (args.index + ".tmp") if args.index else None
    total = 0

    try:
        with open_maybe_gzip(args.gbff) as handle, open(tmp_out, "w") as faa:
            idx = open(tmp_idx, "w") if tmp_idx else None
            if idx:
                idx.write("orf_id\treplicon\tposition\tlocus_tag\tprotein_id\tstart\tend\tstrand\tproduct\n")
            for replicon, cds_list in parse_records(handle):
                cds_list.sort(key=lambda c: (c["start"], c["end"]))
                for rank, cds in enumerate(cds_list, start=1):
                    orf_id = "%s_%d" % (replicon, rank)
                    faa.write(">%s\n%s\n" % (orf_id, cds["translation"].replace("*", "")))
                    total += 1
                    if idx:
                        idx.write("\t".join([orf_id, replicon, str(rank),
                                             cds["locus_tag"] or cds["gene"] or "NA",
                                             cds["protein_id"] or "NA",
                                             str(cds["start"]), str(cds["end"]),
                                             cds["strand"],
                                             cds["product"] or "NA"]) + "\n")
            if idx:
                idx.close()
    except (OSError, EOFError, gzip.BadGzipFile) as err:
        sys.stderr.write("unreadable input %s: %s\n" % (args.gbff, err))
        for path in (tmp_out, tmp_idx):
            if path and os.path.exists(path):
                os.remove(path)
        return 4

    if total < args.min_cds:
        sys.stderr.write("no translated CDS found in %s\n" % args.gbff)
        for path in (tmp_out, tmp_idx):
            if path and os.path.exists(path):
                os.remove(path)
        return 3

    os.replace(tmp_out, args.out)
    if tmp_idx:
        os.replace(tmp_idx, args.index)
    sys.stdout.write("%s\t%d\n" % (os.path.basename(args.gbff), total))
    return 0


if __name__ == "__main__":
    sys.exit(main())
