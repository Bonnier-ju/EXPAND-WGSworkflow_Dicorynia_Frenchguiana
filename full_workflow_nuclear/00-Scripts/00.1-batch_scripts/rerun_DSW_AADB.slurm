#!/bin/bash
#SBATCH -p unlimitq
#SBATCH --time=100:00:00
#SBATCH -J rerun-DSW_AADB
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=48G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# =============================================================================
# Single-sample reprocessing: DSW_AADB
# Steps: fastp (02.1) -> BWA-MEM2 + samtools (03.1) -> HaplotypeCaller (05.1)
# Raw FASTQs must exist before running this script.
# =============================================================================

SAMPLE="DSW_AADB"

# --- Paths -------------------------------------------------------------------
BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
REF="/home/jbonnier/work/chapitre_5/raw_data/reference_genome/nuclear/Dgu_HS1_HYBRID_SCAFFOLD.fa"

RAW_R1="/home/jbonnier/work/chapitre_5/raw_data/ind_dicorynia/EXPAND/DSW_AADBOSDA_7_1_22VLJNLT4.UDI382_clean.fastq.gz"
RAW_R2="/home/jbonnier/work/chapitre_5/raw_data/ind_dicorynia/EXPAND/DSW_AADBOSDA_7_2_22VLJNLT4.UDI382_clean.fastq.gz"

TRIM_DIR="${BASE_DIR}/02-filtering_indexing/workflow_full/cleaned_reads"
TRIM_R1="${TRIM_DIR}/${SAMPLE}_R1.trimmed.fastq.gz"
TRIM_R2="${TRIM_DIR}/${SAMPLE}_R2.trimmed.fastq.gz"
FASTP_LOG="${TRIM_DIR}/logs/${SAMPLE}.fastp.log"

ALIGN_DIR="${BASE_DIR}/03-reading_sorting/workflow_full/alignments"
BAM="${ALIGN_DIR}/${SAMPLE}.sorted.bam"
BAI="${BAM}.bai"
ALIGN_LOG="${ALIGN_DIR}/logs/${SAMPLE}.alignment.log"

GVCF_DIR="${BASE_DIR}/05-variant_calling/workflow_full/haplotypecaller_gvcf"
GVCF="${GVCF_DIR}/${SAMPLE}.g.vcf.gz"
GVCF_LOG="${GVCF_DIR}/logs/${SAMPLE}.haplotypecaller.log"

TMP_DIR="/tmp/${USER}_${SAMPLE}_${SLURM_JOB_ID}"
mkdir -p "$TMP_DIR" "${TRIM_DIR}/logs" "${ALIGN_DIR}/logs" "${GVCF_DIR}/logs"

# --- Thread allocation -------------------------------------------------------
TOTAL="${SLURM_CPUS_PER_TASK:-8}"
SORT_T=$(( TOTAL / 4 ))   # 16/4=4 — more threads left for BWA
[ "$SORT_T" -lt 1 ] && SORT_T=1
VIEW_T=1
BWA_T=$(( TOTAL - SORT_T - VIEW_T ))
[ "$BWA_T" -lt 1 ] && BWA_T=1
PAIRHMM_T=$(( TOTAL - 1 ))
[ "$PAIRHMM_T" -lt 1 ] && PAIRHMM_T=1

echo "======================================================"
echo "START $(date -Iseconds) sample=${SAMPLE}"
echo "  SLURM job ID : ${SLURM_JOB_ID}"
echo "  CPUs         : ${SLURM_CPUS_PER_TASK}"
echo "  Threads BWA  : ${BWA_T} | sort: ${SORT_T} | view: ${VIEW_T}"
echo "  Threads GATK : ${PAIRHMM_T} (PairHMM)"
echo "======================================================"

# =============================================================================
# STEP 1: fastp — quality trimming
# =============================================================================
echo ""
echo "[STEP 1/3] fastp — quality trimming"
echo "  START $(date -Iseconds)"
echo "  Input R1 : $RAW_R1"
echo "  Input R2 : $RAW_R2"
echo "  Output R1: $TRIM_R1"
echo "  Output R2: $TRIM_R2"

if [ -s "$TRIM_R1" ] && [ -s "$TRIM_R2" ]; then
    echo "  SKIP: trimmed FASTQs already exist ($(du -sh "$TRIM_R1" | cut -f1) / $(du -sh "$TRIM_R2" | cut -f1))"
else
    if [ ! -s "$RAW_R1" ] || [ ! -s "$RAW_R2" ]; then
        echo "ERROR: raw FASTQ not found — download incomplete?"
        echo "  R1: $RAW_R1"
        echo "  R2: $RAW_R2"
        exit 1
    fi
    echo "  [OK] raw FASTQs found ($(du -sh "$RAW_R1" | cut -f1) / $(du -sh "$RAW_R2" | cut -f1))"
    module purge
    module load bioinfo/fastp/0.23.4
    echo "  [OK] fastp module loaded"
    fastp \
        --in1  "$RAW_R1" \
        --in2  "$RAW_R2" \
        --out1 "$TRIM_R1" \
        --out2 "$TRIM_R2" \
        --detect_adapter_for_pe \
        --thread "${SLURM_CPUS_PER_TASK:-4}" \
        --trim_poly_g \
        --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --length_required 80 \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 30 \
        --dup_calc_accuracy 2 \
        --html "${TRIM_DIR}/logs/${SAMPLE}.fastp.html" \
        --json "${TRIM_DIR}/logs/${SAMPLE}.fastp.json" \
        2> "$FASTP_LOG"
    echo "  DONE fastp $(date -Iseconds)"
    echo "  Trimmed R1: $(du -sh "$TRIM_R1" | cut -f1)"
    echo "  Trimmed R2: $(du -sh "$TRIM_R2" | cut -f1)"
    echo "  Log: $FASTP_LOG"
fi

# =============================================================================
# STEP 2: BWA-MEM2 + samtools — mapping, sorting, indexing
# =============================================================================
echo ""
echo "[STEP 2/3] BWA-MEM2 + samtools — alignment, sorting, indexing"
echo "  START $(date -Iseconds)"
echo "  Reference : $REF"
echo "  Output BAM: $BAM"
echo "  Threads   : bwa=${BWA_T} view=${VIEW_T} sort=${SORT_T}"

if [ -s "$BAM" ] && [ -s "$BAI" ]; then
    echo "  SKIP: BAM already exists ($(du -sh "$BAM" | cut -f1))"
else
    module purge
    module load bioinfo/bwa-mem2/2.2.1
    module load bioinfo/samtools/1.21
    echo "  [OK] bwa-mem2 + samtools modules loaded"
    RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:EXPAND\tPL:ILLUMINA"
    echo "  Read group: $RG"
    echo "START $(date -Iseconds) bwa=${BWA_T} view=${VIEW_T} sort=${SORT_T}" > "$ALIGN_LOG"
    echo "  Aligning reads..."
    bwa-mem2 mem -M -R "$RG" -t "$BWA_T" "$REF" "$TRIM_R1" "$TRIM_R2" | \
        samtools view -b -@ "$VIEW_T" - | \
        samtools sort -@ "$SORT_T" -m 4G -T "${TMP_DIR}/sort" -o "$BAM" - 2>> "$ALIGN_LOG"
    echo "  Indexing BAM..."
    samtools index "$BAM" 2>> "$ALIGN_LOG"
    echo "DONE $(date -Iseconds) bam=${BAM}" >> "$ALIGN_LOG"
    echo "  DONE BWA-MEM2 $(date -Iseconds)"
    echo "  BAM size  : $(du -sh "$BAM" | cut -f1)"
    echo "  BAI exists: $([ -s "$BAI" ] && echo yes || echo NO)"
    # Quick mapping stats
    module load bioinfo/samtools/1.21
    echo "  Mapping stats:"
    samtools flagstat "$BAM" | grep -E "mapped|total" | sed 's/^/    /'
fi

# =============================================================================
# STEP 3: GATK HaplotypeCaller — gVCF
# =============================================================================
echo ""
echo "[STEP 3/3] GATK HaplotypeCaller — gVCF mode"
echo "  START $(date -Iseconds)"
echo "  BAM input : $BAM"
echo "  GVCF out  : $GVCF"
echo "  Java mem  : 20G | PairHMM threads: ${PAIRHMM_T}"

module purge
module load devel/python/Python-3.11.1
module load devel/java/17.0.6
module load bioinfo/GATK/4.2.6.1
module load bioinfo/samtools/1.21
echo "  [OK] Python + GATK + samtools modules loaded"

# Ensure reference indexes exist
if [ ! -s "${REF}.fai" ]; then
    echo "  INFO: creating reference .fai index..."
    samtools faidx "$REF"
fi
REF_DICT="${REF%.fa}.dict"
if [ ! -s "$REF_DICT" ]; then
    echo "  INFO: creating reference sequence dictionary..."
    gatk CreateSequenceDictionary -R "$REF" -O "$REF_DICT"
fi
echo "  [OK] reference indexes verified"

echo "  Running HaplotypeCaller (this step takes several hours)..."
echo "START $(date -Iseconds) sample=${SAMPLE}" > "$GVCF_LOG"

gatk --java-options "-Xms40G -Xmx40G -XX:ParallelGCThreads=4" HaplotypeCaller \
    -R "$REF" \
    -I "$BAM" \
    -O "$GVCF" \
    -ERC GVCF \
    --native-pair-hmm-threads "$PAIRHMM_T" \
    --create-output-variant-index true \
    --tmp-dir "$TMP_DIR" >> "$GVCF_LOG" 2>&1

echo "DONE $(date -Iseconds) gvcf=${GVCF}" >> "$GVCF_LOG"
echo "  DONE HaplotypeCaller $(date -Iseconds)"
echo "  gVCF size : $(du -sh "$GVCF" | cut -f1)"
echo "  Index     : $([ -s "${GVCF}.tbi" ] && echo yes || echo NO)"

# --- Cleanup -----------------------------------------------------------------
echo ""
echo "Cleaning up tmp directory: $TMP_DIR"
rm -rf "$TMP_DIR"

echo ""
echo "======================================================"
echo "DONE $(date -Iseconds) all steps completed for ${SAMPLE}"
echo "  fastp log  : $FASTP_LOG"
echo "  align log  : $ALIGN_LOG"
echo "  GATK log   : $GVCF_LOG"
echo "  gVCF       : $GVCF"
echo "======================================================"
