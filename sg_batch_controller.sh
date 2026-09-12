#!/usr/bin/env bash
#
# sg_batch_controller.sh - batch controller for screening NCBI RefSeq (GCF_)
# assemblies with SecretionGenie.py, 1000 genomes at a time.
#
# For every batch it will:
#   1. download GenBank flat files (.gb.gz) with bit-dl-ncbi-assemblies  [bit2 env]
#   2. convert each .gb.gz to a SecretionGenie-ready protein FASTA and delete
#      the .gb.gz immediately (so peak disk stays low, and Prodigal is skipped)
#   3. run SecretionGenie in N parallel shards, one thread each   [magiclamp env]
#   4. merge the shard results into the master result files
#   5. delete the batch working directory and move on
#
# It is resumable: finished batches are marked with a .done file and skipped,
# so you can Ctrl-C or lose the machine and just re-run the same command.
#
# Usage:
#   ./sg_batch_controller.sh prep      # build the GCF accession list
#   ./sg_batch_controller.sh run       # process all batches (resumable)
#   ./sg_batch_controller.sh status    # progress + throughput + ETA
#   ./sg_batch_controller.sh pivot     # long counts -> wide heatmap CSV
#   ./sg_batch_controller.sh estimate  # time/disk projection for your machine
#
# Everything below can also be overridden from the environment, e.g.
#   BATCH_SIZE=500 JOBS=24 ./sg_batch_controller.sh run
#
set -euo pipefail

# =====================================================================
# CONFIGURATION
# =====================================================================

# --- inputs ---
ASM_INFO="${ASM_INFO:-$HOME/databases/ncbi_assembly_info.tsv}"
SG_SCRIPT="${SG_SCRIPT:-$HOME/bin/SecretionGenie.py}"
GBFF2FAA="${GBFF2FAA:-$HOME/bin/gbff_to_faa.py}"
TXSS_DIR="${TXSS_DIR:-$HOME/databases/TXSScan}"     # must hold definitions/ and profiles/

# --- where everything goes ---
PROJECT="${PROJECT:-$HOME/databases/secretion_screen}"
WORK="$PROJECT/work"                                 # transient per-batch scratch
RESULTS="$PROJECT/results"                           # kept forever
LOGS="$PROJECT/logs"
HMM_DB="$PROJECT/hmm_db"                             # shared concatenated HMM library

# --- scale knobs ---
BATCH_SIZE="${BATCH_SIZE:-1000}"                     # genomes per batch
JOBS="${JOBS:-24}"                                   # parallel SecretionGenie shards
THREADS_PER_JOB="${THREADS_PER_JOB:-1}"              # hmmsearch --cpu per shard
DL_JOBS="${DL_JOBS:-20}"                             # bit concurrent downloads (NCBI caps at 20)
CONVERT_JOBS="${CONVERT_JOBS:-$JOBS}"
PREFETCH="${PREFETCH:-1}"                            # 1 = download batch N+1 while batch N runs

# --- SecretionGenie options ---
TOPOLOGY="${TOPOLOGY:-linear}"                       # 'circular' only if every replicon is closed
MODELS="${MODELS:-all}"
FAST_HMM="${FAST_HMM:-1}"                            # 1 = one HMM library per genome (identical results)
EXTRA_SG_ARGS="${EXTRA_SG_ARGS:-}"                   # e.g. "-i_evalue_sel 0.001 --genes"
KEEP_ORF_INDEX="${KEEP_ORF_INDEX:-1}"                # 1 = keep ORF -> locus_tag/protein_id maps
                                                     # (~100 MB per 1000 genomes; set 0 to skip)
KEEP_SEQS="${KEEP_SEQS:-1}"                          # 1 = keep the protein sequence column in the
                                                     # master summaries (set 0 to shrink them ~5x)

# --- accession filtering (applied by 'prep') ---
KEEP_GROUPS="${KEEP_GROUPS:-bacteria,archaea}"       # 'all' to disable the taxonomic filter
LATEST_ONLY="${LATEST_ONLY:-1}"                      # keep version_status == latest
ANNOTATED_ONLY="${ANNOTATED_ONLY:-1}"                # require protein_coding_gene_count > 0

# --- conda environments ---
BIT_ENV="${BIT_ENV:-bit2}"
SG_ENV="${SG_ENV:-magiclamp}"
SKIP_CONDA="${SKIP_CONDA:-0}"                        # 1 = run commands directly (envs already active)

# =====================================================================
# derived paths / helpers
# =====================================================================

ACC_LIST="$PROJECT/gcf_accessions.txt"
BATCH_DIR="$PROJECT/batches"
MASTER_COUNTS="$RESULTS/system_counts.tsv.gz"        # accession <tab> system <tab> n_systems
MASTER_STATS="$RESULTS/genome_stats.tsv.gz"          # accession <tab> n_proteins <tab> n_systems
SUMMARY_DIR="$RESULTS/summaries"                     # per-batch gene-level summaries (gz)
INDEX_DIR="$RESULTS/orf_index"                       # per-batch ORF -> locus_tag maps (gz)
PROBLEMS="$RESULTS/problem_accessions.tsv"           # accession <tab> reason

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOGS/controller.log"; }
die() { log "ERROR: $*"; exit 1; }

in_bit() {
    if [[ "$SKIP_CONDA" == "1" ]]; then "$@"; else
        conda run --no-capture-output -n "$BIT_ENV" "$@"; fi
}
in_sg() {
    if [[ "$SKIP_CONDA" == "1" ]]; then "$@"; else
        conda run --no-capture-output -n "$SG_ENV" "$@"; fi
}

secs_to_hms() { printf '%dh%02dm%02ds' $(( $1/3600 )) $(( ($1%3600)/60 )) $(( $1%60 )); }

preflight() {
    local problems=0

    if [[ ! -f "$SG_SCRIPT" ]]; then
        log "ERROR: SecretionGenie not found at: $SG_SCRIPT"
        log "       set SG_SCRIPT=/full/path/to/SecretionGenie.py (or edit it at the top of this script)"
        problems=1
    elif [[ "$FAST_HMM" == "1" ]] && ! grep -q "fast_hmm" "$SG_SCRIPT"; then
        log "ERROR: $SG_SCRIPT is an older SecretionGenie without the --fast_hmm / -hmm_db_dir"
        log "       / --build_hmm_db options, so the shared HMM library cannot be used."
        log "       Argparse would abbreviate -hmm_db_dir to -h and die with"
        log "       'argument -h/--help: ignored explicit argument mm_db_dir'."
        log "       Fix it either way:"
        log "         - install the updated SecretionGenie.py at $SG_SCRIPT   (recommended), or"
        log "         - re-run with FAST_HMM=0, which uses one hmmsearch per profile"
        log "           (identical results, about 12% slower)"
        problems=1
    fi

    if [[ ! -f "$GBFF2FAA" ]]; then
        log "ERROR: gbff_to_faa.py not found at: $GBFF2FAA  (set GBFF2FAA=/full/path)"
        problems=1
    fi

    if [[ ! -d "$TXSS_DIR/definitions" || ! -d "$TXSS_DIR/profiles" ]]; then
        log "ERROR: $TXSS_DIR must contain both definitions/ and profiles/"
        log "       git clone https://github.com/macsy-models/TXSScan.git \"$TXSS_DIR\""
        problems=1
    fi

    [[ "$problems" -eq 0 ]] || exit 1
    log "preflight ok: SecretionGenie=$SG_SCRIPT  fast_hmm=$FAST_HMM  TXSScan=$TXSS_DIR"
}

# =====================================================================
# prep : build the list of RefSeq (GCF_) accessions to screen
# =====================================================================

cmd_prep() {
    mkdir -p "$PROJECT" "$LOGS" "$RESULTS" "$SUMMARY_DIR" "$INDEX_DIR" "$BATCH_DIR"
    [[ -s "$ASM_INFO" ]] || die "assembly info file not found: $ASM_INFO"
    preflight

    log "building GCF accession list from $ASM_INFO"
    log "  filters: groups=$KEEP_GROUPS latest_only=$LATEST_ONLY annotated_only=$ANNOTATED_ONLY"

    # Column map of the NCBI assembly summary table:
    #   1 assembly_accession   11 version_status   18 gbrs_paired_asm
    #  25 group                36 protein_coding_gene_count
    # The GenBank-flavoured table carries GCA_ in column 1 and the paired RefSeq
    # GCF_ accession in column 18; the RefSeq-flavoured table already has GCF_ in
    # column 1. Both are handled here.
    awk -F'\t' -v groups="$KEEP_GROUPS" -v latest="$LATEST_ONLY" -v annot="$ANNOTATED_ONLY" '
        BEGIN {
            OFS="\t"
            if (groups != "all") {
                n = split(groups, g, ",")
                for (i = 1; i <= n; i++) { gsub(/^[ \t]+|[ \t]+$/, "", g[i]); keep[g[i]] = 1 }
                use_groups = 1
            }
        }
        /^#/ { next }
        {
            acc = ""
            if ($1  ~ /^GCF_/) acc = $1
            else if ($18 ~ /^GCF_/) acc = $18
            if (acc == "") next
            if (latest == 1 && $11 != "latest") next
            if (use_groups && !($25 in keep)) next
            if (annot == 1 && $36 != "" && $36 !~ /^[0-9]+$/) next
            if (annot == 1 && $36 ~ /^[0-9]+$/ && $36 + 0 == 0) next
            print acc
        }
    ' "$ASM_INFO" | LC_ALL=C sort -u > "$ACC_LIST.tmp"

    mv "$ACC_LIST.tmp" "$ACC_LIST"
    local total; total=$(wc -l < "$ACC_LIST")
    log "wrote $total GCF accessions to $ACC_LIST"

    # split into fixed-size batches, zero-padded so they sort naturally
    rm -f "$BATCH_DIR"/batch_*.acc
    split -l "$BATCH_SIZE" -d -a 5 --additional-suffix=.acc "$ACC_LIST" "$BATCH_DIR/batch_"
    local nbatches; nbatches=$(ls "$BATCH_DIR"/batch_*.acc | wc -l)
    log "split into $nbatches batches of up to $BATCH_SIZE genomes"

    if [[ "$FAST_HMM" == "1" ]]; then
        log "building the shared concatenated HMM library once (used by all shards)"
        mkdir -p "$HMM_DB"
        local probe="$WORK/hmmdb_probe"
        mkdir -p "$probe/bins"
        in_sg python3 "$SG_SCRIPT" -bin_dir "$probe/bins" -bin_ext gb -out "$probe/out" \
            -txss_dir "$TXSS_DIR" -models "$MODELS" -hmm_db_dir "$HMM_DB" --build_hmm_db \
            >> "$LOGS/hmm_db_build.log" 2>&1 || die "could not build the HMM library (see $LOGS/hmm_db_build.log)"
        rm -rf "$probe"
        log "HMM library ready in $HMM_DB"
    else
        log "FAST_HMM=0: skipping the shared HMM library, each shard will search one profile at a time"
    fi
    log "now run:  $0 run"
}

# =====================================================================
# per-batch steps
# =====================================================================

download_batch() {                      # $1 = batch acc file, $2 = destination dir
    local acc_file="$1" dest="$2"
    mkdir -p "$dest"
    log "  downloading $(wc -l < "$acc_file") GenBank files -> $dest"
    in_bit bit-dl-ncbi-assemblies -w "$acc_file" -f genbank -j "$DL_JOBS" -o "$dest" \
        >> "$LOGS/download.log" 2>&1 || log "  WARNING: downloader returned non-zero (continuing with what arrived)"
    if [[ -s "$dest/ncbi-accessions-not-found.txt" ]]; then
        awk '{print $1"\tnot_found_at_ncbi"}' "$dest/ncbi-accessions-not-found.txt" >> "$PROBLEMS"
    fi
    touch "$dest/.downloaded"
}

process_batch() {                       # $1 = batch id, $2 = acc file, $3 = download dir
    local bid="$1" acc_file="$2" dl_dir="$3"
    local bwork="$WORK/$bid"
    rm -rf "$bwork"; mkdir -p "$bwork"

    # ---------- 2. GenBank -> protein FASTA, then drop the GenBank file -------
    local ncpu_conv="$CONVERT_JOBS"
    log "  converting GenBank -> proteins with $ncpu_conv workers (Prodigal skipped)"
    mkdir -p "$bwork/faa" "$bwork/idx"
    : > "$bwork/converted.tsv"
    find "$dl_dir" -maxdepth 1 -name '*.gb.gz' -print0 |
        xargs -0 -P "$ncpu_conv" -I{} bash -c '
            gb="$1"; faa_dir="$2"; idx_dir="$3"; conv="$4"; probl="$5"; script="$6"
            acc=$(basename "$gb" .gb.gz)
            if python3 "$script" "$gb" -o "$faa_dir/$acc.gb-proteins.faa" \
                                       --index "$idx_dir/$acc.idx" >> "$conv" 2>/dev/null; then
                rm -f "$gb"
            else
                printf "%s\tno_translated_CDS_or_unreadable\n" "$acc" >> "$probl"
                rm -f "$gb"
            fi
        ' _ {} "$bwork/faa" "$bwork/idx" "$bwork/converted.tsv" "$PROBLEMS" "$GBFF2FAA"

    local n_faa; n_faa=$(find "$bwork/faa" -name '*-proteins.faa' | wc -l)
    log "  $n_faa proteomes ready"
    [[ "$n_faa" -gt 0 ]] || { log "  nothing to process in $bid"; rm -rf "$bwork" "$dl_dir"; return 0; }

    # ---------- 3. shard and run SecretionGenie in parallel ------------------
    log "  running SecretionGenie: $JOBS shards x $THREADS_PER_JOB thread(s)"
    local i=0 shard
    while IFS= read -r faa; do
        shard=$(( i % JOBS )); i=$(( i + 1 ))
        mkdir -p "$bwork/shard$shard/bins" "$bwork/shard$shard/out/ORF_calls"
        local base; base=$(basename "$faa" -proteins.faa)      # e.g. GCF_000005845.2.gb
        mv "$faa" "$bwork/shard$shard/out/ORF_calls/$base-proteins.faa"
        : > "$bwork/shard$shard/bins/$base"                     # placeholder: ORF calls already exist
    done < <(find "$bwork/faa" -name '*-proteins.faa' | LC_ALL=C sort)

    local fast_flag="" ; [[ "$FAST_HMM" == "1" ]] && fast_flag="--fast_hmm"
    seq 0 $(( JOBS - 1 )) | xargs -P "$JOBS" -I{} bash -c '
        s="$1"; bwork="$2"; sg="$3"; txss="$4"; models="$5"; topo="$6"; thr="$7"
        fast="$8"; hmmdb="$9"; extra="${10}"; logs="${11}"; skip="${12}"; env="${13}"
        [[ -d "$bwork/shard$s/bins" ]] || exit 0
        ls -A "$bwork/shard$s/bins" | grep -q . || exit 0
        cmd=(python3 "$sg" -bin_dir "$bwork/shard$s/bins" -bin_ext gb
             -out "$bwork/shard$s/out" -txss_dir "$txss" -models "$models"
             -topology "$topo" -t "$thr")
        # -hmm_db_dir only exists in the SecretionGenie versions that support
        # --fast_hmm, so it is passed only when fast mode is actually requested
        [[ -n "$fast" ]] && cmd+=("$fast" -hmm_db_dir "$hmmdb")
        [[ -n "$extra" ]] && read -r -a extra_arr <<< "$extra" && cmd+=("${extra_arr[@]}")
        if [[ "$skip" == "1" ]]; then "${cmd[@]}"; else conda run --no-capture-output -n "$env" "${cmd[@]}"; fi \
            > "$logs/shard$s.log" 2>&1 || echo "shard $s failed, see $logs/shard$s.log" >&2
    ' _ {} "$bwork" "$SG_SCRIPT" "$TXSS_DIR" "$MODELS" "$TOPOLOGY" "$THREADS_PER_JOB" \
      "$fast_flag" "$HMM_DB" "$EXTRA_SG_ARGS" "$LOGS" "$SKIP_CONDA" "$SG_ENV"

    # ---------- 4. merge shard results into the master files -----------------
    log "  merging results"
    local merged="$bwork/merged-summary.csv"
    local header="file,ORF,gene,function,gene_status,system,system_id,replicon,position,locus,hit_type,evalue,bit_score,profile_coverage,sequence_coverage,system_score,system_wholeness,nb_loci,seq"
    printf '%s\n' "$header" > "$merged"
    local sfile
    for sfile in "$bwork"/shard*/out/secretiongenie-summary.csv; do
        [[ -f "$sfile" ]] || continue
        # drop the per-shard header, strip the '.gb' placeholder suffix from the
        # genome column so the master tables carry clean GCF_ accessions
        awk -F, 'NR>1 { if ($0 ~ /^#/) { print; next } sub(/\.gb$/, "", $1); sub(/\.gb_/, "_", $7); print }' OFS=, "$sfile" >> "$merged"
    done
    if [[ "$KEEP_SEQS" == "1" ]]; then
        gzip -c "$merged" > "$SUMMARY_DIR/$bid.secretiongenie-summary.csv.gz"
    else
        awk -F, 'BEGIN{OFS=","} { if ($0 ~ /^#/) { print; next } NF=NF-1; print }' "$merged" |
            gzip -c > "$SUMMARY_DIR/$bid.secretiongenie-summary.csv.gz"
    fi

    # long-format counts: accession <tab> system <tab> number of systems
    awk -F, 'NR>1 && $0 !~ /^#/ && $1 != "" { key = $1 FS $6 FS $7; if (!(key in seen)) { seen[key]=1; n[$1 FS $6]++ } }
             END { for (k in n) { split(k, a, FS); print a[1] "\t" a[2] "\t" n[k] } }' \
        "$merged" | LC_ALL=C sort > "$bwork/counts.tsv"
    gzip -c "$bwork/counts.tsv" >> "$MASTER_COUNTS"

    # per-genome stats, including the genomes where nothing was found (0 systems)
    : > "$bwork/stats.tsv"
    while IFS= read -r base; do
        local acc="${base%.gb}"
        local faa="$bwork"/shard*/out/ORF_calls/"$base"-proteins.faa
        local nprot=0
        for f in $faa; do [[ -f "$f" ]] && nprot=$(grep -c '^>' "$f" || true); done
        local nsys; nsys=$(awk -F'\t' -v a="$acc" '$1==a {s+=$3} END {print s+0}' "$bwork/counts.tsv")
        printf '%s\t%s\t%s\n' "$acc" "$nprot" "$nsys" >> "$bwork/stats.tsv"
    done < <(find "$bwork"/shard*/bins -maxdepth 1 -type f -printf '%f\n' 2>/dev/null | LC_ALL=C sort)
    gzip -c "$bwork/stats.tsv" >> "$MASTER_STATS"

    # keep the ORF -> locus_tag / protein_id map so hits can be traced back
    if [[ "$KEEP_ORF_INDEX" == "1" ]] && compgen -G "$bwork/idx/*.idx" > /dev/null; then
        tar -czf "$INDEX_DIR/$bid.orf_index.tar.gz" -C "$bwork/idx" . 2>/dev/null || true
    fi

    # ---------- 5. clean up --------------------------------------------------
    local n_sys; n_sys=$(awk -F'\t' '{s+=$3} END {print s+0}' "$bwork/counts.tsv")
    log "  batch $bid: $n_faa genomes, $n_sys systems detected"
    rm -rf "$bwork" "$dl_dir"
    touch "$BATCH_DIR/$bid.done"
}

# =====================================================================
# run : the main loop
# =====================================================================

cmd_run() {
    [[ -s "$ACC_LIST" ]] || die "no accession list yet - run '$0 prep' first"
    preflight
    mkdir -p "$WORK" "$RESULTS" "$LOGS" "$SUMMARY_DIR" "$INDEX_DIR"
    touch "$PROBLEMS"

    local batches=() f
    while IFS= read -r f; do batches+=("$f"); done < <(ls "$BATCH_DIR"/batch_*.acc | LC_ALL=C sort)
    local total=${#batches[@]}
    log "=== starting run: $total batches, BATCH_SIZE=$BATCH_SIZE JOBS=$JOBS PREFETCH=$PREFETCH ==="

    local run_start=$SECONDS done_now=0
    local idx bid acc_file dl_dir next_bid next_acc next_dl prefetch_pid=""

    for (( idx = 0; idx < total; idx++ )); do
        acc_file="${batches[$idx]}"
        bid=$(basename "$acc_file" .acc)
        dl_dir="$WORK/dl_$bid"

        if [[ -f "$BATCH_DIR/$bid.done" ]]; then
            continue
        fi

        local t0=$SECONDS
        log "--- $bid ($((idx+1))/$total) ---"

        # wait for a prefetched download, or download now
        if [[ -n "$prefetch_pid" ]]; then
            wait "$prefetch_pid" || log "  WARNING: prefetch download exited non-zero"
            prefetch_pid=""
        fi
        [[ -f "$dl_dir/.downloaded" ]] || download_batch "$acc_file" "$dl_dir"

        # kick off the next batch's download in the background
        if [[ "$PREFETCH" == "1" && $(( idx + 1 )) -lt $total ]]; then
            next_acc="${batches[$((idx+1))]}"
            next_bid=$(basename "$next_acc" .acc)
            next_dl="$WORK/dl_$next_bid"
            if [[ ! -f "$BATCH_DIR/$next_bid.done" && ! -f "$next_dl/.downloaded" ]]; then
                download_batch "$next_acc" "$next_dl" &
                prefetch_pid=$!
            fi
        fi

        process_batch "$bid" "$acc_file" "$dl_dir"

        done_now=$(( done_now + 1 ))
        local dt=$(( SECONDS - t0 ))
        local elapsed=$(( SECONDS - run_start ))
        local remaining=$(( total - idx - 1 ))
        local eta=$(( remaining * elapsed / done_now ))
        log "  batch time $(secs_to_hms $dt) | elapsed $(secs_to_hms $elapsed) | $remaining batches left | ETA $(secs_to_hms $eta)"
    done

    [[ -n "$prefetch_pid" ]] && wait "$prefetch_pid" || true
    log "=== all batches finished in $(secs_to_hms $(( SECONDS - run_start ))) ==="
    log "run '$0 pivot' to build the wide heatmap CSV"
}

# =====================================================================
# status / pivot / clean / estimate
# =====================================================================

cmd_status() {
    local total done_n
    total=$(ls "$BATCH_DIR"/batch_*.acc 2>/dev/null | wc -l)
    done_n=$(ls "$BATCH_DIR"/*.done 2>/dev/null | wc -l)
    echo "batches:        $done_n / $total done"
    if [[ -f "$MASTER_STATS" ]]; then
        local ngen nsys
        ngen=$(zcat "$MASTER_STATS" | wc -l)
        nsys=$(zcat "$MASTER_STATS" | awk -F'\t' '{s+=$3} END {print s+0}')
        echo "genomes done:   $ngen"
        echo "systems found:  $nsys"
        echo "systems/genome: $(zcat "$MASTER_STATS" | awk -F'\t' '{s+=$3; n++} END {if(n) printf "%.2f\n", s/n}')"
    fi
    [[ -f "$PROBLEMS" ]] && echo "problem accs:   $(wc -l < "$PROBLEMS")"
    echo "disk in results: $(du -sh "$RESULTS" 2>/dev/null | cut -f1)"
    if [[ -f "$MASTER_COUNTS" ]]; then
        echo
        echo "top systems so far:"
        zcat "$MASTER_COUNTS" | awk -F'\t' '{n[$2]+=$3} END {for (s in n) printf "  %-14s %d\n", s, n[s]}' | sort -k2,2nr
    fi
}

cmd_pivot() {
    [[ -f "$MASTER_COUNTS" ]] || die "no counts yet: $MASTER_COUNTS"
    local out="$RESULTS/secretiongenie.heatmap.csv"
    log "pivoting long counts -> $out"
    zcat "$MASTER_COUNTS" | awk -F'\t' '
        { n[$1 SUBSEP $2] = $3; genomes[$1] = 1; systems[$2] = 1 }
        END {
            ns = 0; for (s in systems) sys_list[++ns] = s
            # deterministic system order
            for (i = 1; i < ns; i++) for (j = i+1; j <= ns; j++)
                if (sys_list[j] < sys_list[i]) { t = sys_list[i]; sys_list[i] = sys_list[j]; sys_list[j] = t }
            printf "X"
            ng = 0; for (g in genomes) g_list[++ng] = g
            for (i = 1; i < ng; i++) for (j = i+1; j <= ng; j++)
                if (g_list[j] < g_list[i]) { t = g_list[i]; g_list[i] = g_list[j]; g_list[j] = t }
            for (j = 1; j <= ng; j++) printf ",%s", g_list[j]
            printf "\n"
            for (i = 1; i <= ns; i++) {
                printf "%s", sys_list[i]
                for (j = 1; j <= ng; j++) {
                    k = g_list[j] SUBSEP sys_list[i]
                    printf ",%d", (k in n) ? n[k] : 0
                }
                printf "\n"
            }
        }' > "$out"
    log "wrote $out ($(wc -l < "$out") rows)"
    log "note: with >~20k genomes prefer the long table $MASTER_COUNTS for plotting in R/pandas"
}

cmd_clean() {
    log "removing transient work directories under $WORK (results are untouched)"
    rm -rf "$WORK"/dl_* "$WORK"/batch_* 2>/dev/null || true
}

cmd_estimate() {
    local per_genome_cpu="${PER_GENOME_CPU:-8}"    # CPU-seconds per genome, 280 profiles, ~4300 proteins
    local n
    n=$( [[ -s "$ACC_LIST" ]] && wc -l < "$ACC_LIST" || echo "${N_GENOMES:-1000}" )
    local compute=$(( n * per_genome_cpu / JOBS ))
    local dl_mb=$(( n * 4 ))
    echo "genomes to screen:        $n"
    echo "assumed CPU-s per genome: $per_genome_cpu (measured on a 2.9 GHz core, 280 profiles)"
    echo "parallel shards:          $JOBS"
    echo "compute time:             ~$(secs_to_hms $compute)"
    echo "download volume:          ~$(( dl_mb / 1024 )) GB of .gb.gz"
    echo "peak disk per batch:      ~$(( BATCH_SIZE * 4 / 1024 + BATCH_SIZE * 2 / 1024 )) GB (GenBank + proteomes)"
    echo
    echo "per batch of $BATCH_SIZE: compute ~$(secs_to_hms $(( BATCH_SIZE * per_genome_cpu / JOBS )))"
    echo "                         download ~$(secs_to_hms $(( BATCH_SIZE * 4 / 15 )))  (assuming ~15 MB/s aggregate)"
    echo "with PREFETCH=1 the two overlap, so expect the slower of the two per batch."
}

# =====================================================================
# entry point
# =====================================================================

mkdir -p "$LOGS"
case "${1:-}" in
    prep)     cmd_prep ;;
    run)      cmd_run ;;
    status)   cmd_status ;;
    pivot)    cmd_pivot ;;
    clean)    cmd_clean ;;
    estimate) cmd_estimate ;;
    *)
        sed -n '2,40p' "$0"
        echo
        echo "commands: prep | run | status | pivot | clean | estimate"
        exit 1 ;;
esac
