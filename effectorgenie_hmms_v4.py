#!/usr/bin/env python3
"""Build candidate VFDB protein-family HMMs and search proteomes (CPU only).

INSTALL
  mamba create -n effectorgenie -c conda-forge -c bioconda python=3.11 mmseqs2 mafft hmmer matplotlib
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
  The new build also writes calibrated/ with provisional embedded GA/TC/NC.
  For an EXISTING build (no clustering/alignment repeated):
    python effectorgenie_hmms.py calibrate --db vfdb_hmms -o vfdb_calibrated -t 24
  Search that library with:
    python effectorgenie_hmms.py search proteins.faa --db vfdb_calibrated \
        -o hits --cut-ga -t 24
  --cut-ga or --cut-nc still applies coverage and independent E-value filters.
  --cut-nc uses embedded NC1 (sequence) and NC2 (domain) scores, inclusively.
  Example: search proteins.faa --db vfdb_calibrated -o hits --cut-nc --evalue 1e-5
  --evalue is an alias for --max-i-evalue. Both mean domain independent E-value.
  A hit must pass BOTH embedded score thresholds AND E-value AND coverage.
  Direct hmmsearch --cut_ga uses embedded scores without those extra filters.

AUTOMATIC CUTOFFS (PROVISIONAL, NOT INDEPENDENT VALIDATION)
  Query each HMM against unique.faa (unaligned raw reference sequences).
  Family members are positives; all other VFDB sequences are BACKGROUND,
  including sequences from small clusters that have no HMM. Background is
  not experimentally established noise. An HMM can recognize real homologs
  in other clusters; this procedure can over-restrict family membership.
  Full sequence score and best domain score are estimated separately.
  NC = maximum observed background score (floored at 0 if none is positive).
  TC = minimum positive score among positives exceeding BOTH background
       maxima; thus the trusted subset excludes overlapping positives.
  GA = midpoint between NC and TC, rounded to 0.01 bits.
  Every embedded pair is: full-sequence score, domain score.
  Profiles without a trusted subset or with missing positive scores are
  withheld in unresolved_profiles/; no invented TC is assigned. Overlap is
  flagged, and positive retention is reported in calibration.tsv.
  Default HMMER acceleration is used. --calibration-max disables heuristic
  filters and can be VERY slow for thousands of profiles; even this is not
  independent validation. Scores are reported down to -1000 bits. No domain
  coverage or E-value filtering is applied during cutoff estimation.
  Original build profiles are preserved. New builds can skip calibration
  with --skip-calibration. Calibration output requires a new directory.
  plots/ contains a two-panel PNG per HMM: full-sequence and best-domain
  histograms, family positives versus background, and GA/TC/NC lines.
  Frequencies use a log(count+1) axis and fixed 5-bit bins. Unreported hits
  are counted in captions, NOT imputed as zero. Plots include unresolved
  models, and use family IDs/VFDB IDs as labels; category here means a
  sequence family, not a broad virulence-function class. --no-plots skips.

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
from collections import defaultdict, Counter

# v3: Excel categorization, Genie-style summaries, optional candidate neighborhoods.
# Install xlrd for .xls and openpyxl for .xlsx: mamba install -c conda-forge xlrd openpyxl
# search ... --vf-categories VFs.xls --sample EXACT_SECRETIONGENIE_FILE_VALUE --ordered-ids
# --ordered-ids ONLY for contig_geneRank IDs (Prodigal or gbff_to_faa output).
# Otherwise provide --gene-map genes.tsv with ORF,replicon,position columns
# (optional start,end,strand) or leave locations unknown. Numeric accession
# suffixes MUST NOT be treated as gene ranks. Use the SAME FAA as SecretionGenie.
# annotate existing_hits.tsv --fasta proteins.faa --db vfdb_calibrated -o annotated \
#   --vf-categories VFs.xls --ordered-ids --sample genome.fna
# Output summary's first 19 columns match SecretionGenie. Secretion-system
# fields are blank. cluster_id identifies proximity neighborhoods, not operons.
# Groups use consecutive hit positions with <= --cluster-gap intervening genes
# (default 5), linear contigs only, ignoring category and strand. Singletons
# remain present. Unplaced proteins receive separate blocks. No topology wrap.
# Summary: one row per protein/family, highest-scoring passing domain shown.
# hits.tsv retains every passing domain. Category heatmap counts distinct ORFs
# per category; a protein in multiple categories contributes once to each.
# Excel annotations are family-level associations, not query-specific proof.
# Missing VFIDs remain unmapped; original FASTA headers are retained separately.

# v4: search assembly.fna --input-type nucleotide --prodigal-mode meta ...
# Requires Prodigal (mamba install -c bioconda prodigal). Accepts uncompressed
# assembled nucleotide FASTA, not raw sequencing reads. Default input remains
# proteins for backward compatibility. Nucleotide mode automatically records
# contig, CDS rank, 1-based inclusive coordinates and strand. Predicted proteins,
# GFF and gene map are retained in ORF_calls/ for reuse with SecretionGenie.
# Use matching Prodigal mode in both tools; reuse the generated FAA to ensure
# exact matching IDs and gene calls. --sample defaults to assembly basename.


def predict_query_proteins(a, root):
    require(['prodigal'])
    if a.gene_map or a.ordered_ids:
        raise ValueError('Nucleotide mode derives locations automatically; omit --gene-map and --ordered-ids')
    original = Path(a.fasta).resolve()
    seen = set()
    for header, seq in fasta(original):
        if not header.strip():
            raise ValueError('Empty nucleotide FASTA header')
        name = header.split()[0]
        if name in seen:
            raise ValueError(f'Duplicate contig ID: {name}')
        seen.add(name)
        if not seq or re.search(r'[^ACGTRYSWKMBDHVN]', seq):
            raise ValueError(f'Invalid nucleotide sequence in {name}; provide assembled DNA FASTA')
    if not seen:
        raise ValueError('No nucleotide FASTA records')
    calls = root/'ORF_calls'
    calls.mkdir()
    proteins = calls/(original.name+'-proteins.faa')
    gff = calls/(original.name+'-genes.gff')
    run(['prodigal','-i',original,'-a',proteins,'-o',gff,'-f','gff',
         '-p',a.prodigal_mode,'-q'], calls/'prodigal.log')
    rows=[]
    for header, seq in fasta(proteins):
        parts=header.split(' # ')
        if len(parts)<4:
            raise ValueError(f'Unexpected Prodigal protein header: {header}')
        name=parts[0].split()[0]
        match=re.fullmatch(r'(.+)_(\d+)',name)
        if not match or match[1] not in seen:
            raise ValueError(f'Cannot map Prodigal ORF to input contig: {name}')
        start,end,strand=int(parts[1]),int(parts[2]),int(parts[3])
        if not 1<=start<=end or strand not in [-1,1]:
            raise ValueError(f'Invalid Prodigal coordinates: {name}')
        rows.append(dict(ORF=name,replicon=match[1],position=int(match[2]),
                         start=start,end=end,strand='+' if strand==1 else '-'))
    if not rows:
        raise ValueError('Prodigal predicted no proteins; inspect ORF_calls/prodigal.log')
    gene_map=calls/'gene_map.tsv'
    tsv(gene_map,['ORF','replicon','position','start','end','strand'],rows)
    result=argparse.Namespace(**vars(a))
    result.input_fasta=str(original)
    result.fasta=str(proteins)
    result.gene_map=str(gene_map)
    result.sample=a.sample or original.name
    print(f'Prodigal predicted {len(rows)} proteins; saved to {proteins}',flush=True)
    return result


def load_vf_categories(path):
    if not path:
        return {}
    suffix = Path(path).suffix.lower()
    tables = []
    if suffix == '.xls':
        try:
            import xlrd
        except ImportError:
            raise ValueError('Reading .xls requires xlrd: mamba install -c conda-forge xlrd')
        book = xlrd.open_workbook(path)
        tables = [[sheet.row_values(i) for i in range(sheet.nrows)] for sheet in book.sheets()]
    elif suffix == '.xlsx':
        try:
            import openpyxl
        except ImportError:
            raise ValueError('Reading .xlsx requires openpyxl')
        book = openpyxl.load_workbook(path, read_only=True, data_only=True)
        tables = [list(sheet.values) for sheet in book.worksheets]
        book.close()
    else:
        with open(path, encoding='utf-8-sig', newline='') as handle:
            tables = [list(csv.reader(handle, delimiter='\t' if suffix=='.tsv' else ','))]
    mapping = defaultdict(list)
    for rows in tables:
        header = None
        for row in rows:
            values = [str(x).strip() if x is not None else '' for x in row]
            if 'VFID' in values and 'VFcategory' in values:
                header = values
                continue
            if header is None:
                continue
            entry = dict(zip(header, values))
            key = entry.get('VFID', '')
            if re.fullmatch(r'VF\d+', key):
                selected = {k:entry.get(k, '') for k in ['VFID','VF_Name','VF_FullName','VFCID','VFcategory','Function','Mechanism']}
                if selected not in mapping[key]:
                    mapping[key].append(selected)
    if not mapping:
        raise ValueError('No VFID/VFcategory table found in category file')
    return dict(mapping)


ANNOTATION_FIELDS = ['VF_Name','VFcategory','Function','Mechanism','mapping_status',
                     'unmapped_vf_ids','vf_annotations_json','reference_annotations']


def family_annotations(db, mapping):
    with open(Path(db)/'families.tsv') as h:
        families = list(csv.DictReader(h, delimiter='\t'))
    references = defaultdict(set)
    members = Path(db)/'members.tsv'
    if members.is_file():
        with open(members) as h:
            for row in csv.DictReader(h, delimiter='\t'):
                # Preserve distinct descriptive headers, minus accession and organism.
                header = row.get('header','')
                desc = header.split(' ',1)[1] if ' ' in header else header
                desc = re.sub(r'\s*\[[^\[\]]*\]\s*$', '', desc)
                if desc:
                    references[row['family']].add(desc)
    result = {}
    for family in families:
        ids = sorted(set(filter(None, family['vf_ids'].split(';'))))
        entries = [entry for vf in ids for entry in mapping.get(vf, [])]
        missing = [vf for vf in ids if vf not in mapping]
        result[family['family']] = {
            **{k:' | '.join(sorted({e[k] for e in entries if e[k]})) for k in ['VF_Name','VFcategory','Function','Mechanism']},
            'mapping_status': ('unmapped' if not entries else 'partial' if missing else 'mapped'),
            'unmapped_vf_ids': ';'.join(missing),
            'vf_annotations_json': json.dumps(entries, ensure_ascii=False),
            'reference_annotations': ' | '.join(sorted(references[family['family']]))}
    return result


def query_metadata(a):
    proteins, loc = {}, {}
    for header, seq in fasta(a.fasta):
        key = header.split()[0]
        if key in proteins:
            raise ValueError(f'Duplicate query protein ID: {key}')
        proteins[key] = seq.rstrip('*')
        loc[key] = dict(replicon='', position='', start='', end='', strand='', location_source='unknown')
        if a.ordered_ids:
            match = re.fullmatch(r'(.+)_(\d+)', key)
            if not match or int(match[2]) < 1:
                raise ValueError(f'--ordered-ids requires contig_positiveRank IDs: {key}')
            loc[key].update(replicon=match[1], position=int(match[2]), location_source='ordered_id')
    if a.gene_map:
        with open(a.gene_map, encoding='utf-8-sig', newline='') as h:
            reader = csv.DictReader(h, delimiter='\t' if Path(a.gene_map).suffix.lower()=='.tsv' else ',')
            if not {'ORF','replicon','position'} <= set(reader.fieldnames or []):
                raise ValueError('Gene map requires ORF, replicon, position columns')
            seen = set()
            for row in reader:
                key = row['ORF']
                if key in seen:
                    raise ValueError(f'Duplicate gene-map ORF: {key}')
                seen.add(key)
                if key not in proteins:
                    continue
                pos = int(row['position'])
                if pos < 1 or not row['replicon']:
                    raise ValueError(f'Invalid location: {key}')
                loc[key] = dict(replicon=row['replicon'], position=pos,
                                start=row.get('start',''), end=row.get('end',''), strand=row.get('strand',''),
                                location_source='gene_map')
    placed = [(x['replicon'],x['position']) for x in loc.values() if x['position']!='']
    if len(set(placed)) != len(placed):
        raise ValueError('Multiple ORFs have the same replicon/gene-rank location')
    return proteins, loc


def write_genie_outputs(a, root, annotations, proteins, locations):
    sample = a.sample or Path(a.fasta).name
    best = {}
    with open(root/'hits.tsv') as h:
        for row in csv.DictReader(h, delimiter='\t'):
            key = row['protein_id'], row['family']
            if key[0] not in proteins:
                raise ValueError(f'Hit protein absent from query FASTA: {key[0]}')
            rank = (float(row['dom_bits']), -float(row['i_evalue']), float(row['coverage']))
            if key not in best or rank > best[key][0]:
                best[key] = rank,row
    ids = sorted({key[0] for key in best}, key=lambda p:(locations[p]['replicon'],
                 locations[p]['position'] if locations[p]['position']!='' else math.inf,p))
    clusters, last = {}, None
    block = 0
    for key in ids:
        loc = locations[key]
        if (last is None or loc['position']=='' or last['position']=='' or
            loc['replicon']!=last['replicon'] or loc['position']-last['position']-1 > a.cluster_gap):
            block += 1
        clusters[key] = f'{sample}:EGC{block:06d}'
        last = loc
    fields = ['file','ORF','gene','function','gene_status','system','system_id','replicon','position','locus',
              'hit_type','evalue','bit_score','profile_coverage','sequence_coverage','system_score',
              'system_wholeness','nb_loci','seq','cluster_id','cluster_kind','start','end','strand','location_source',
              'vf_ids','seq_bits','hmm_from','hmm_to','ali_from','ali_to'] + ANNOTATION_FIELDS
    rows, counts = [], defaultdict(set)
    by_protein = defaultdict(list)
    for (protein,family), (_,hit) in best.items():
        by_protein[protein].append((family,hit))
    for protein in ids:
        loc = locations[protein]
        for _, hit in sorted(by_protein[protein]):
            annotation = annotations[hit['family']]
            row = dict.fromkeys(fields,'')
            row.update(file=sample, ORF=protein, gene=hit['family'], function=annotation['VF_Name'] or hit['family'],
                       gene_status='candidate',
                       locus=clusters[protein], hit_type='family_homolog', evalue=hit['i_evalue'],
                       bit_score=hit['dom_bits'], profile_coverage=hit['coverage'],
                       sequence_coverage=(int(hit['ali_to'])-int(hit['ali_from'])+1)/len(proteins[protein]),
                       seq=proteins[protein], cluster_id=clusters[protein],
                       cluster_kind='proximity_neighborhood' if loc['position']!='' else 'unplaced',
                       **loc, **annotation)
            for name in ['vf_ids','seq_bits','hmm_from','hmm_to','ali_from','ali_to']:
                row[name] = hit[name]
            rows.append(row)
            categories = {e['VFcategory'] for e in json.loads(annotation['vf_annotations_json']) if e['VFcategory']}
            if annotation['unmapped_vf_ids'] or not categories:
                categories.add('Unmapped')
            for category in categories:
                counts[category].add(protein)
    # Valid rectangular CSV for joins, plus SecretionGenie's visual separators.
    for filename,separators in [('effectorgenie-all_hits.csv',False),('effectorgenie-summary.csv',True)]:
        with open(root/filename,'w',newline='') as out:
            writer=csv.DictWriter(out,fieldnames=fields); writer.writeheader()
            previous=None
            for row in rows:
                if separators and previous is not None and previous!=row['cluster_id']:
                    out.write('####################################################\n')
                writer.writerow(row); previous=row['cluster_id']
            if separators and rows:
                out.write('####################################################\n')
    all_categories = {'Unmapped'} | {e['VFcategory'] for ann in annotations.values()
                       for e in json.loads(ann['vf_annotations_json']) if e['VFcategory']}
    with open(root/'effectorgenie.heatmap.csv','w',newline='') as out:
        writer=csv.writer(out); writer.writerow(['X',sample])
        for category in sorted(all_categories):
            writer.writerow([category,len(counts[category])])
    print(f'Genie summary: {len(ids)} proteins, {block} candidate neighborhoods/unplaced blocks; file={sample}')


def annotate_existing(a):
    mapping = load_vf_categories(a.vf_categories)
    annotations = family_annotations(a.db,mapping)
    proteins, locations = query_metadata(a)
    root = Path(a.out).resolve(); root.mkdir(parents=True,exist_ok=False)
    with open(a.hits) as inp, open(root/'hits.tsv','w',newline='') as out:
        reader=csv.DictReader(inp,delimiter='\t')
        required={'protein_id','family','vf_ids','seq_bits','dom_bits','i_evalue','coverage','hmm_from','hmm_to','ali_from','ali_to'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError('Input must be EffectorGenie hits.tsv')
        fields=[f for f in reader.fieldnames if f not in ANNOTATION_FIELDS]+ANNOTATION_FIELDS
        writer=csv.DictWriter(out,fieldnames=fields,delimiter='\t');writer.writeheader()
        for row in reader:
            if row['family'] not in annotations:
                raise ValueError(f'Unknown family: {row["family"]}')
            writer.writerow(dict(row,**annotations[row['family']]))
    (root/'annotation_settings.json').write_text(json.dumps(vars(a),indent=2))
    write_genie_outputs(a,root,annotations,proteins,locations)


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
        require(['mmseqs', 'mafft', 'hmmbuild'] + ([] if a.skip_calibration else ['hmmsearch']))
        if not a.skip_calibration and not a.no_plots:
            try:
                import matplotlib
            except ImportError:
                raise ValueError('Install matplotlib before building, or use --no-plots')
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
    if not a.skip_calibration and cutoffs:
        calibrate(argparse.Namespace(db=str(root), out=str(root/'calibrated'),
                                    threads=a.threads, calibration_max=a.calibration_max,
                                    no_plots=a.no_plots))


def estimate_thresholds(positive, background, expected):
    """positive: ID -> (sequence, best domain); background: maxima or None."""
    nc = tuple(max(0.0, x) if x is not None else 0.0 for x in background)
    base = dict(n_positive_expected=expected, n_positive_scored=len(positive),
                nc_seq=nc[0], nc_dom=nc[1], tc_seq='', tc_dom='', ga_seq='', ga_dom='',
                n_positive_retained=0, n_positive_trusted=0,
                background_seq_observed=background[0] is not None,
                background_dom_observed=background[1] is not None)
    if len(positive) != expected:
        return dict(base, status='unresolved_missing_positive_scores'), None
    # HMMER describes TC as a true-positive score above the known noise.
    trusted = [s for s in positive.values() if s[0] > nc[0] and s[1] > nc[1]]
    if not trusted:
        return dict(base, status='unresolved_no_positive_above_background'), None
    tc = (min(s[0] for s in trusted), min(s[1] for s in trusted))
    ga = tuple(round((n+t)/2, 2) for n, t in zip(nc, tc))
    if any(not n < g <= t for n, g, t in zip(nc, ga, tc)):
        return dict(base, status='unresolved_rounding_gap'), None
    retained = sum(s[0] >= ga[0] and s[1] >= ga[1] for s in positive.values())
    base.update(tc_seq=tc[0], tc_dom=tc[1], ga_seq=ga[0], ga_dom=ga[1],
                n_positive_retained=retained, n_positive_trusted=len(trusted))
    status = 'separated' if len(trusted)==expected else 'overlap_restrictive'
    return dict(base, status=status), {'GA': ga, 'TC': tc, 'NC': nc}


def embedded_profile(text, tags):
    """Insert HMMER ASCII optional header fields before the HMM matrix."""
    lines = text.splitlines(keepends=True)
    result, inserted = [], False
    for line in lines:
        if not inserted and re.match(r'^(GA|TC|NC)\s', line):
            continue
        if line.startswith('HMM ') and not inserted:
            result.extend(f'{tag:<6}{tags[tag][0]:.2f} {tags[tag][1]:.2f};\n'
                          for tag in ['GA', 'TC', 'NC'])
            inserted = True
        result.append(line)
    if not inserted or not text.rstrip().endswith('//'):
        raise ValueError('Invalid or incomplete HMM profile')
    return ''.join(result)


def calibrate(a):
    require(['hmmsearch'])
    if not a.no_plots:
        try:
            import matplotlib
        except ImportError:
            raise ValueError('Install matplotlib for histograms, or use --no-plots')
    db, root = Path(a.db).resolve(), Path(a.out).resolve()
    for name in ['unique.faa', 'members.tsv', 'families.tsv', 'library.hmm']:
        if not (db/name).is_file():
            raise ValueError(f'Missing calibration input: {db/name}')
    with open(db/'families.tsv') as handle:
        families = [r for r in csv.DictReader(handle, delimiter='\t') if r['status']=='built']
    if not families:
        raise ValueError('No built profiles')
    membership = defaultdict(set)
    with open(db/'members.tsv') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            membership[row['family']].add(row['representative'])
    for f in families:
        if len(membership[f['family']]) != int(f['n_unique']):
            raise ValueError(f'Inconsistent membership: {f["family"]}')
    root.mkdir(parents=True, exist_ok=False)
    (root/'profiles').mkdir()
    (root/'unresolved_profiles').mkdir()
    config = dict(vars(a), source_library_sha256=hashlib.sha256((db/'library.hmm').read_bytes()).hexdigest(),
                  source_sequences_sha256=hashlib.sha256((db/'unique.faa').read_bytes()).hexdigest(),
                  cutoff_method='NC_background_max_TC_trusted_subset_GA_midpoint', provisional=True)
    (root/'calibration_settings.json').write_text(json.dumps(config, indent=2))
    print('Scoring HMM library against raw VFDB references; this can take time.', flush=True)
    cmd = ['hmmsearch', '--cpu', a.threads, '--noali', '-o', '/dev/null',
           '-T', '-1000', '--domT', '-1000', '--tblout', root/'reference.tbl',
           '--domtblout', root/'reference.domtbl']
    if a.calibration_max:
        cmd.append('--max')
    run(cmd + [db/'library.hmm', db/'unique.faa'], root/'calibration.log')
    positive_seq, positive_dom = defaultdict(dict), defaultdict(dict)
    negative_seq, negative_dom = {}, {}
    hist_seq, hist_dom = defaultdict(Counter), defaultdict(Counter)
    total_references = sum(1 for _ in fasta(db/'unique.faa'))
    # Sequence scores are read independently: a background sequence need not
    # have a reported domain to contribute to the sequence noise maximum.
    with open(root/'reference.tbl') as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            f = line.split()
            sid, family, score = f[0], f[2], float(f[5])
            if sid in membership[family]:
                positive_seq[family][sid] = score
            else:
                negative_seq[family] = max(negative_seq.get(family, -math.inf), score)
                hist_seq[family][math.floor(score/5)] += 1
    # Reduce domain rows to one best domain per protein, for correct histogram
    # counts. HMMER emits each model-target's domain records contiguously.
    previous, best = None, -math.inf
    def record_background_domain(key, value):
        if key is not None:
            family, sid = key
            if sid not in membership[family]:
                hist_dom[family][math.floor(value/5)] += 1
    for row in domain_rows(root/'reference.domtbl'):
        sid, family, score = row['protein_id'], row['family'], row['dom_bits']
        key = (family, sid)
        if key != previous:
            record_background_domain(previous, best)
            previous, best = key, score
        else:
            best = max(best, score)
        if sid in membership[family]:
            positive_dom[family][sid] = max(positive_dom[family].get(sid, -math.inf), score)
        else:
            negative_dom[family] = max(negative_dom.get(family, -math.inf), score)
    record_background_domain(previous, best)
    reports, accepted, cutoffs = [], [], []
    with open(root/'library.hmm', 'w') as library:
        for family_number, f in enumerate(families, 1):
            family = f['family']
            positives = {sid: (score, positive_dom[family][sid])
                         for sid, score in positive_seq[family].items() if sid in positive_dom[family]}
            report, tags = estimate_thresholds(positives,
                (negative_seq.get(family), negative_dom.get(family)), len(membership[family]))
            reports.append(dict(family=family, **report))
            if family_number % 50 == 0:
                print(f'Calibrating/plotting {family_number}/{len(families)} profiles', flush=True)
            if not a.no_plots:
                plot_scores(root, f, report, tags,
                            positive_seq[family], positive_dom[family],
                            hist_seq[family], hist_dom[family], total_references,
                            a.calibration_max)
            source = db/'profiles'/f'{family}.hmm'
            if tags is None:
                shutil.copy2(source, root/'unresolved_profiles'/source.name)
                continue
            text = embedded_profile(source.read_text(), tags)
            (root/'profiles'/source.name).write_text(text)
            library.write(text)
            accepted.append(f)
            cutoffs.append(dict(family=family, seq_bits=tags['GA'][0], dom_bits=tags['GA'][1],
                                coverage='', max_i_evalue=''))
    tsv(root/'families.tsv', list(families[0]), accepted)
    tsv(root/'calibration.tsv', list(reports[0]), reports)
    tsv(root/'cutoffs.tsv', ['family', 'seq_bits', 'dom_bits', 'coverage', 'max_i_evalue'], cutoffs)
    # Keep annotation provenance available to downstream joins.
    shutil.copy2(db/'members.tsv', root/'members.tsv')
    print(f'Embedded cutoffs: {len(accepted)} profiles; unresolved: {len(families)-len(accepted)}. '
          f'See {root / "calibration.tsv"}. These are provisional training/reference cutoffs.')


def plot_scores(root, family, report, tags, pos_seq, pos_dom, neg_seq, neg_dom, total, exhaustive):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    folder = root/'plots'
    folder.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    expected = int(family['n_unique'])
    for index, (ax, pos, neg, label) in enumerate(zip(axes, [pos_seq, pos_dom],
            [neg_seq, neg_dom], ['Full-sequence bit score', 'Best-domain bit score'])):
        positives = Counter(math.floor(v/5) for v in pos.values())
        bins = set(positives) | set(neg)
        if bins:
            low, high = min(bins), max(bins)
            edges = [5*x for x in range(low, high+2)]
            ax.stairs([neg.get(x, 0) for x in range(low, high+1)], edges,
                      color='#d97732', label='Other VFDB proteins', linewidth=1.4)
            ax.stairs([positives.get(x, 0) for x in range(low, high+1)], edges,
                      color='#176a9b', label='Family members', linewidth=1.7)
        if tags:
            for tag, color, style in [('GA', '#49318c', '-'), ('TC', '#227342', '--'), ('NC', '#a52c35', ':')]:
                ax.axvline(tags[tag][index], color=color, linestyle=style,
                           label=f'{tag} {tags[tag][index]:.2f}', linewidth=1.2)
        else:
            ax.axvline(report['nc_seq' if index==0 else 'nc_dom'], color='#a52c35',
                       linestyle=':', label='Background NC estimate')
        ax.set_ylim(bottom=0)
        ax.set_yscale('function', functions=(np.log1p, np.expm1))
        ax.set_ylabel('Protein count (log1p scale)')
        ax.set_xlabel(label)
        ax.set_title(f'Reported: {len(pos)}/{expected} family; {sum(neg.values())}/{total-expected} other', fontsize=10)
        ax.grid(axis='y', alpha=.2)
        ax.legend(fontsize=7)
    fig.suptitle(f'{family["family"]} | {family["vf_ids"][:110]}\n{report["status"]} — provisional reference-derived thresholds', fontsize=11)
    mode = 'Heuristic filters disabled' if exhaustive else 'Heuristic filters enabled'
    fig.text(.5, .02, f'{mode}; 5-bit bins. Unreported proteins are not plotted. Self-hits are not independent validation.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .05, 1, .90))
    fig.savefig(folder/f'{family["family"]}.png', dpi=130)
    plt.close(fig)


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
    if a.input_type=='nucleotide':
        require(['prodigal'])
        if a.gene_map or a.ordered_ids:
            raise ValueError('Nucleotide mode derives locations automatically; omit --gene-map and --ordered-ids')
    mapping = load_vf_categories(a.vf_categories)
    annotations = family_annotations(a.db, mapping)
    db = Path(a.db).resolve()
    if not (db/'library.hmm').is_file() or not (db/'library.hmm').stat().st_size:
        raise ValueError('No HMM library; build profiles first')
    with open(db/'families.tsv') as handle:
        families = {r['family']: r for r in csv.DictReader(handle, delimiter='\t') if r['status']=='built'}
    defaults = dict(seq_bits=a.seq_bits, dom_bits=a.dom_bits, coverage=a.coverage,
                    max_i_evalue=a.max_i_evalue)
    thresholds = {k: defaults.copy() for k in families}
    embedded_tag = 'GA' if a.cut_ga else 'NC' if a.cut_nc else None
    if embedded_tag:
        if a.cutoffs or a.seq_bits != 0 or a.dom_bits != 0:
            raise ValueError('Embedded cutoffs cannot be combined with score overrides or --cutoffs')
        for family in families:
            text = (db/'profiles'/f'{family}.hmm').read_text()
            match = re.search(r'^'+embedded_tag+r'\s+([\d.eE+-]+)\s+([\d.eE+-]+);', text, re.M)
            if not match:
                raise ValueError(f'Missing {embedded_tag}: {family}')
            thresholds[family].update(seq_bits=float(match[1]), dom_bits=float(match[2]))
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
    if a.input_type=='nucleotide':
        a=predict_query_proteins(a,root)
    proteins, locations = query_metadata(a)
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
              'coverage', 'hmm_from', 'hmm_to', 'ali_from', 'ali_to'] + ANNOTATION_FIELDS
    count = 0
    with open(root/'hits.tsv', 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for row in domain_rows(root/'raw.domtbl'):
            if passes(row, thresholds[row['family']]):
                row['vf_ids'] = families[row['family']]['vf_ids']
                row.update(annotations[row['family']])
                writer.writerow(row)
                count += 1
    print(f'{count} passing domains written to {root / "hits.tsv"}')
    write_genie_outputs(a, root, annotations, proteins, locations)


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
    b.add_argument('--skip-calibration', action='store_true')
    b.add_argument('--calibration-max', action='store_true', help='Disable heuristic scoring filters (very slow)')
    b.add_argument('--no-plots', action='store_true')
    c = sub.add_parser('calibrate', help='Embed provisional cutoffs in copies of existing profiles; make histograms')
    c.add_argument('--db', required=True)
    c.add_argument('-o', '--out', required=True)
    c.add_argument('-t', '--threads', type=int, default=8)
    c.add_argument('--calibration-max', action='store_true')
    c.add_argument('--no-plots', action='store_true')
    s = sub.add_parser('search', help='Search a proteome with explicit per-family filters')
    s.add_argument('fasta', help='Protein FASTA, or assembled DNA FASTA with --input-type nucleotide')
    s.add_argument('--input-type',choices=['proteins','nucleotide'],default='proteins')
    s.add_argument('--prodigal-mode',choices=['meta','single'],default='meta',
                   help='Gene prediction mode for nucleotide input (default meta; single for suitable isolate genomes)')
    s.add_argument('--db', required=True)
    s.add_argument('-o', '--out', required=True)
    s.add_argument('-t', '--threads', type=int, default=8)
    s.add_argument('--seq-bits', type=float, default=0)
    s.add_argument('--dom-bits', type=float, default=0)
    s.add_argument('--coverage', type=float, default=0.70)
    s.add_argument('--max-i-evalue', '--evalue', dest='max_i_evalue', type=float, default=1e-3,
                   help='Maximum domain independent E-value; additional to score/coverage checks (default 0.001)')
    s.add_argument('--cutoffs')
    embedded = s.add_mutually_exclusive_group()
    embedded.add_argument('--cut-ga', action='store_true', help='Use embedded GA scores plus E-value/coverage filters')
    embedded.add_argument('--cut-nc', action='store_true', help='Use embedded NC scores plus E-value/coverage filters')
    s.add_argument('--max', action='store_true', help='Disable HMMER heuristic acceleration filters (slow)')
    n = sub.add_parser('annotate', help='Categorize existing hits and create Genie summaries without HMMER')
    n.add_argument('hits')
    n.add_argument('--fasta', required=True)
    n.add_argument('--db', required=True)
    n.add_argument('-o','--out',required=True)
    for parser in [s,n]:
        parser.add_argument('--vf-categories', help='VFDB category table: .xls, .xlsx, .csv or .tsv')
        parser.add_argument('--sample', help='Exact file value used in SecretionGenie summary; default FASTA basename')
        positions=parser.add_mutually_exclusive_group()
        positions.add_argument('--ordered-ids', action='store_true', help='Explicitly treat contig_geneRank IDs as locations')
        positions.add_argument('--gene-map', help='TSV/CSV with ORF,replicon,position and optional start,end,strand')
        parser.add_argument('--cluster-gap',type=int,default=5,help='Maximum intervening genes between successive VF hits (default 5)')
    a = p.parse_args()
    if getattr(a,'cluster_gap',0)<0:
        p.error('cluster-gap must be >=0')
    if getattr(a,'threads',1) < 1:
        p.error('threads must be >=1')
    if a.command == 'build':
        if not (0 < a.identity <= 1 and 0 < a.cluster_coverage <= 1 and a.min_members >= 2):
            p.error('identity and coverage must be in (0,1]; min-members must be >=2')
        build(a)
    elif a.command == 'calibrate':
        calibrate(a)
    elif a.command == 'annotate':
        annotate_existing(a)
    else:
        search(a)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        sys.exit(f'ERROR: {exc}')
