#!/usr/bin/env python3
"""Build candidate VFDB protein-family HMMs and search proteomes (CPU only).

INSTALL
  mamba create -n effectorgenie -c conda-forge -c bioconda python=3.11 mmseqs2 mafft hmmer
  conda activate effectorgenie

BUILD (new output directory required; --prepare-only needs no external tools)
  python effectorgenie_hmms.py build VFDB_setB_pro.fas -o vfdb_hmms -t 24

SEARCH (initial, UNVALIDATED screening settings)
  python effectorgenie_hmms.py search proteins.faa --db vfdb_hmms -o sample_hits -t 24
  python effectorgenie_hmms.py search proteins.faa --db vfdb_hmms -o stricter_hits \
      --seq-bits 40 --dom-bits 40 --coverage 0.8

PER-FAMILY CUTOFFS
  Edit vfdb_hmms/cutoffs.tsv, then add --cutoffs vfdb_hmms/cutoffs.tsv to search.
  Blank cells inherit command-line defaults. All supplied filters must pass.
  seq_bits = full-sequence score; dom_bits = individual domain score;
  coverage = (hmm_to - hmm_from + 1) / profile length, per domain;
  max_i_evalue = domain independent E-value, NOT conditional E-value.
  No score threshold is inferred from training-set self-hits.

DESIGN / LIMITATIONS
  Global MMseqs2 clustering: 30% identity, 80% bidirectional coverage,
  greedy representative clustering. These are heuristic starting settings,
  not proof of common function. VF annotations are preserved, never used as
  mandatory sequence-family boundaries. Exact duplicate sequences collapse
  before clustering; every original record remains in members.tsv.
  MAFFT --auto -> hmmbuild for clusters with >=3 unique sequences.
  Small clusters remain in unmodeled.faa for optional DIAMOND searches.
  All VFDB categories are retained: this is a VF homolog library, NOT an
  effector-specific library. Review annotations/alignments, split ambiguous
  clusters, and curate an effector subset before assigning effector labels.
  Related families may be split; multiple HMMs may match the same protein.
  Domains/rearranged proteins need separate domain-aware curation.
  Near-duplicate sampling can inflate apparent support despite deduplication.
  Building profiles is NOT validating their functional specificity.

CALIBRATION
  For each family, separate reference sequences into sequence-similarity
  clusters BEFORE splitting training and validation sets. Build a development
  HMM on training only. Score independent known positives and hard negatives
  (especially related proteins with different functions). Sweep sequence and
  domain score thresholds, assess sensitivity and false positives, and choose
  a cutoff meeting your specificity goal; assess on a separate test set.
  Score overlap may require better families or domain/context criteria.
  Absence from VFDB does not establish a negative label. Small families may
  lack enough evidence for calibration: retain their provisional status.
  Rebuilding an HMM changes its scores; revalidate after profile changes.
  Higher library size / broader genome screens warrant false-positive checks.

OUTPUTS
  references.tsv: original headers, VF IDs, duplicate representatives
  families.tsv: candidate family size, labels, HMM status
  members.tsv: every original protein -> family
  families/*.faa, alignments/*.afa, profiles/*.hmm, library.hmm
  cutoffs.tsv: blank per-family threshold template (not calibrated)
  run.json: input SHA256, parameters; logs/: exact external commands + logs
  Search: hits.tsv (one row per passing domain), raw.domtbl, search.log,
          settings.json; family annotations carried into hits.tsv.

No automatic resume: use a fresh output directory after failure. Inspect
logs before retrying. No external Python packages required.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
from collections import defaultdict


def fasta(path):
    header, seq = None, []
    # VFDB includes legacy non-UTF8 bytes in some descriptions.
    with open(path, encoding='latin-1') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq)
                header, seq = line[1:], []
            elif header is None:
                raise ValueError('Sequence before first FASTA header')
            else:
                seq.append(''.join(line.split()).upper())
    if header is not None:
        yield header, ''.join(seq)


def write_fasta(path, records):
    with open(path, 'w') as out:
        for name, seq in records:
            out.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + '\n')


def tsv(path, fields, rows):
    with open(path, 'w', newline='') as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def run(cmd, log, output=None):
    cmd = [str(x) for x in cmd]
    with open(log, 'w') as err:
        err.write(json.dumps(cmd) + '\n')
        err.flush()
        if output:
            with open(output, 'w') as out:
                subprocess.run(cmd, stdout=out, stderr=err, check=True)
        else:
            subprocess.run(cmd, stdout=err, stderr=err, check=True)


def require(names):
    missing = [x for x in names if not shutil.which(x)]
    if missing:
        raise ValueError('Missing executables: ' + ', '.join(missing))


def build(a):
    if not a.prepare_only:
        require(['mmseqs', 'mafft', 'hmmbuild'])
    root = Path(a.out).resolve()
    root.mkdir(parents=True, exist_ok=False)
    for sub in ['logs', 'families', 'alignments', 'profiles']:
        (root/sub).mkdir()
    records, sequences, seq_to_id, seen = [], {}, {}, set()
    for header, seq in fasta(a.fasta):
        if not header.strip():
            raise ValueError('Empty FASTA header')
        original = header.split()[0]
        if original in seen:
            raise ValueError(f'Duplicate FASTA ID: {original}')
        seen.add(original)
        seq = seq.rstrip('*')
        if not seq or re.search(r'[^ACDEFGHIKLMNPQRSTVWYBXZJUO]', seq):
            raise ValueError(f'Empty or invalid protein sequence: {original}')
        if seq not in seq_to_id:
            seq_to_id[seq] = f'S{len(sequences)+1:07d}'
            sequences[seq_to_id[seq]] = seq
        records.append(dict(protein_id=original, representative=seq_to_id[seq],
                            vf_ids=';'.join(sorted(set(re.findall(r'\bVF\d+\b', header)))),
                            header=header))
    if not records:
        raise ValueError('No FASTA records')
    tsv(root/'references.tsv', ['protein_id', 'representative', 'vf_ids', 'header'], records)
    write_fasta(root/'unique.faa', sequences.items())
    metadata = dict(vars(a))
    metadata['input_sha256'] = hashlib.sha256(Path(a.fasta).read_bytes()).hexdigest()
    metadata.update(n_records=len(records), n_unique=len(sequences))
    (root/'run.json').write_text(json.dumps(metadata, indent=2))
    print(f'{len(records):,} records; {len(sequences):,} unique sequences', flush=True)
    if a.prepare_only:
        return
    run(['mmseqs', 'easy-cluster', root/'unique.faa', root/'clusters', root/'tmp',
         '--min-seq-id', a.identity, '-c', a.cluster_coverage, '--cov-mode', '0',
         '--cluster-mode', '2', '--threads', a.threads, '-s', '7.5'], root/'logs/clustering.log')
    clusters, assigned = defaultdict(list), set()
    for line in (root/'clusters_cluster.tsv').read_text().splitlines():
        rep, member = line.split('\t')
        if member not in sequences or member in assigned:
            raise ValueError(f'Invalid cluster assignment: {member}')
        assigned.add(member)
        clusters[rep].append(member)
    if assigned != set(sequences):
        raise ValueError('Clustering omitted input sequences')
    refs = defaultdict(list)
    for row in records:
        refs[row['representative']].append(row)
    family_rows, member_rows, small, cutoffs = [], [], [], []
    with open(root/'library.hmm', 'w') as library:
        for i, rep in enumerate(sorted(clusters), 1):
            family = f'EG{i:06d}'
            ids = sorted(clusters[rep])
            originals = [r for sid in ids for r in refs[sid]]
            family_faa = root/'families'/f'{family}.faa'
            write_fasta(family_faa, ((sid, sequences[sid]) for sid in ids))
            vf_ids = sorted({v for r in originals for v in r['vf_ids'].split(';') if v})
            status = 'built' if len(ids) >= a.min_members else 'too_few_unique_sequences'
            family_rows.append(dict(family=family, n_unique=len(ids), n_records=len(originals),
                                    vf_ids=';'.join(vf_ids), status=status))
            member_rows.extend(dict(family=family, **r) for r in originals)
            if status != 'built':
                small.extend((f'{sid} family={family}', sequences[sid]) for sid in ids)
                continue
            aln = root/'alignments'/f'{family}.afa'
            hmm = root/'profiles'/f'{family}.hmm'
            run(['mafft', '--thread', a.threads, '--auto', family_faa],
                root/'logs'/f'{family}.mafft.log', aln)
            run(['hmmbuild', '--amino', '--cpu', a.threads, '-n', family, hmm, aln],
                root/'logs'/f'{family}.hmmbuild.log')
            library.write(hmm.read_text())
            cutoffs.append(dict(family=family, seq_bits='', dom_bits='', coverage='', max_i_evalue=''))
            if i % 50 == 0:
                print(f'Processed {i}/{len(clusters)} clusters', flush=True)
    tsv(root/'families.tsv', ['family', 'n_unique', 'n_records', 'vf_ids', 'status'], family_rows)
    tsv(root/'members.tsv', ['family', 'protein_id', 'representative', 'vf_ids', 'header'], member_rows)
    tsv(root/'cutoffs.tsv', ['family', 'seq_bits', 'dom_bits', 'coverage', 'max_i_evalue'], cutoffs)
    write_fasta(root/'unmodeled.faa', small)
    print(f'Done: {len(clusters)} candidate families, {len(cutoffs)} HMMs. Cutoffs remain unvalidated.')


def domain_rows(path):
    with open(path) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            f = line.split(maxsplit=22)
            if len(f) < 22:
                raise ValueError('Malformed HMMER domtblout row')
            yield dict(protein_id=f[0], family=f[3], seq_bits=float(f[7]),
                       dom_bits=float(f[13]), i_evalue=float(f[12]),
                       coverage=(int(f[16])-int(f[15])+1)/int(f[5]),
                       hmm_from=int(f[15]), hmm_to=int(f[16]),
                       ali_from=int(f[17]), ali_to=int(f[18]))


def passes(row, threshold):
    return (row['seq_bits'] >= threshold['seq_bits'] and
            row['dom_bits'] >= threshold['dom_bits'] and
            row['coverage'] >= threshold['coverage'] and
            row['i_evalue'] <= threshold['max_i_evalue'])


def search(a):
    require(['hmmsearch'])
    db = Path(a.db).resolve()
    if not (db/'library.hmm').is_file() or not (db/'library.hmm').stat().st_size:
        raise ValueError('No HMM library; build profiles first')
    with open(db/'families.tsv') as handle:
        families = {r['family']: r for r in csv.DictReader(handle, delimiter='\t') if r['status']=='built'}
    defaults = dict(seq_bits=a.seq_bits, dom_bits=a.dom_bits, coverage=a.coverage,
                    max_i_evalue=a.max_i_evalue)
    thresholds = {k: defaults.copy() for k in families}
    if a.cutoffs:
        seen = set()
        with open(a.cutoffs) as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                family = row['family']
                if family not in families or family in seen:
                    raise ValueError(f'Unknown or duplicate cutoff family: {family}')
                seen.add(family)
                for key in defaults:
                    value = (row.get(key) or '').strip()
                    if value:
                        thresholds[family][key] = float(value)
    for t in thresholds.values():
        if (not all(math.isfinite(v) for v in t.values()) or t['seq_bits'] < 0 or
                t['dom_bits'] < 0 or not 0 < t['coverage'] <= 1 or t['max_i_evalue'] <= 0):
            raise ValueError('Scores must be >=0; coverage in (0,1]; E-value >0; all finite')
    root = Path(a.out).resolve()
    root.mkdir(parents=True, exist_ok=False)
    settings = dict(vars(a), thresholds=thresholds,
                    library_sha256=hashlib.sha256((db/'library.hmm').read_bytes()).hexdigest())
    (root/'settings.json').write_text(json.dumps(settings, indent=2))
    # Report nonnegative scores first; apply actual thresholds below.
    cmd = ['hmmsearch', '--cpu', a.threads, '--noali', '-T', '0', '--domT', '0',
           '--domtblout', root/'raw.domtbl']
    if a.max:
        cmd.append('--max')  # Slower: bypass HMMER's heuristic filters.
    run(cmd + [db/'library.hmm', Path(a.fasta).resolve()], root/'search.log')
    fields = ['protein_id', 'family', 'vf_ids', 'seq_bits', 'dom_bits', 'i_evalue',
              'coverage', 'hmm_from', 'hmm_to', 'ali_from', 'ali_to']
    count = 0
    with open(root/'hits.tsv', 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for row in domain_rows(root/'raw.domtbl'):
            if passes(row, thresholds[row['family']]):
                row['vf_ids'] = families[row['family']]['vf_ids']
                writer.writerow(row)
                count += 1
    print(f'{count} passing domains written to {root / "hits.tsv"}')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    b = sub.add_parser('build', help='Cluster reference sequences, align, and build HMMs')
    b.add_argument('fasta')
    b.add_argument('-o', '--out', required=True)
    b.add_argument('-t', '--threads', type=int, default=8)
    b.add_argument('--identity', type=float, default=0.30)
    b.add_argument('--cluster-coverage', type=float, default=0.80)
    b.add_argument('--min-members', type=int, default=3)
    b.add_argument('--prepare-only', action='store_true')
    s = sub.add_parser('search', help='Search a proteome with explicit per-family filters')
    s.add_argument('fasta')
    s.add_argument('--db', required=True)
    s.add_argument('-o', '--out', required=True)
    s.add_argument('-t', '--threads', type=int, default=8)
    s.add_argument('--seq-bits', type=float, default=0)
    s.add_argument('--dom-bits', type=float, default=0)
    s.add_argument('--coverage', type=float, default=0.70)
    s.add_argument('--max-i-evalue', type=float, default=1e-3)
    s.add_argument('--cutoffs')
    s.add_argument('--max', action='store_true', help='Disable HMMER heuristic acceleration filters (slow)')
    a = p.parse_args()
    if a.threads < 1:
        p.error('threads must be >=1')
    if a.command == 'build':
        if not (0 < a.identity <= 1 and 0 < a.cluster_coverage <= 1 and a.min_members >= 2):
            p.error('identity and coverage must be in (0,1]; min-members must be >=2')
        build(a)
    else:
        search(a)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        sys.exit(f'ERROR: {exc}')
