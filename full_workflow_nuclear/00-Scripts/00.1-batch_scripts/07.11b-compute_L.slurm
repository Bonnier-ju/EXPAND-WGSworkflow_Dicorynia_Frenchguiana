#!/bin/bash
# 07.11b-compute_L.slurm
#
# Estimate L (callable non-TE nucleotide sites) from 10 representative GVCFs
# for use in Stairway Plot 2 blueprint files (07.11-stairway_plot.slurm).
#
# L = number of sites genotyped with GQ >= 20 outside TE-masked regions.
# Ten individuals spread across the 4 K=4 ADMIXTURE clusters (3+3+2+2),
# chosen by evenly-spaced sampling of sorted sample IDs within each cluster.
#
# Pipeline per GVCF:
#   bcftools view -R non_te.bed -i 'GQ>=20'    → restrict + quality filter
#   awk                                          → convert to 0-based BED
#   bedtools intersect -a stdin -b non_te.bed   → clip reference blocks at
#                                                  TE boundaries
#   awk sum                                      → callable sites for this sample
#
# Output:
#   07.11-stairway_plot/07.11b-compute_L/L_per_sample.tsv
#   07.11-stairway_plot/07.11b-compute_L/07.11b-compute_L.log
#   STDOUT: recommended L value to replace L_CALLABLE in 07.11-stairway_plot.slurm

#SBATCH -p workq
#SBATCH --time=08:00:00
#SBATCH -J 07.11b-compute_L
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment
# -------------------------------------------------------------------
module purge
module load bioinfo/Bcftools/1.9
module load bioinfo/bedtools/2.31.1

# -------------------------------------------------------------------
# STEP 1: Parameters and paths
# -------------------------------------------------------------------
K_GLOBAL=4
GQ_MIN=20

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
GVCF_DIR="${BASE}/05-variant_calling/workflow_full/haplotypecaller_gvcf"
CLUSTER_ASSIGNMENTS="${BASE}/metadata/cluster_assignments.tsv"
TE_GFF="/home/jbonnier/work/chapitre_5/raw_data/reference_genome/nuclear/Dgu_HS1_HYBRID_SCAFFOLD.fa.out.gff.gz"
GENOME_FAI="/home/jbonnier/work/chapitre_5/raw_data/reference_genome/nuclear/Dgu_HS1_HYBRID_SCAFFOLD.fa.fai"

OUT_DIR="${BASE}/07-population_structure/workflow_full/07.11-stairway_plot/07.11b-compute_L"
LOG="${OUT_DIR}/07.11b-compute_L.log"

mkdir -p "${OUT_DIR}"
echo "START $(date -Iseconds) — 07.11b compute_L" | tee "${LOG}"

# -------------------------------------------------------------------
# STEP 2: Build non-TE BED from RepeatMasker GFF + genome FAI
# RepeatMasker GFF: 1-based inclusive coords (cols 4 and 5).
# BED is 0-based half-open → start = min(col4,col5) - 1, end = max(col4,col5).
# -------------------------------------------------------------------
TE_BED="${OUT_DIR}/te_merged.bed"
NON_TE_BED="${OUT_DIR}/non_te.bed"
GENOME_SIZES="${OUT_DIR}/genome.sizes"

echo "STEP 2: Building non-TE BED ..." | tee -a "${LOG}"

# Genome sizes from FAI (chr \t length), sorted lexicographically
awk '{print $1"\t"$2}' "${GENOME_FAI}" | sort -k1,1 > "${GENOME_SIZES}"

TOTAL_GENOME=$(awk '{sum+=$2} END{print sum}' "${GENOME_SIZES}")
echo "  Genome total: ${TOTAL_GENOME} bp" | tee -a "${LOG}"

# TE regions: sort, merge overlapping intervals
zcat "${TE_GFF}" | awk '!/^#/ {
    s = ($4 < $5) ? $4 - 1 : $5 - 1
    e = ($4 < $5) ? $5     : $4
    print $1"\t"s"\t"e
}' | sort -k1,1 -k2,2n | bedtools merge -i - > "${TE_BED}"

N_TE=$(awk '{sum+=$3-$2} END{print sum+0}' "${TE_BED}")
echo "  TE-masked bases: ${N_TE} bp" | tee -a "${LOG}"

# Non-TE complement (requires TE_BED and GENOME_SIZES sorted consistently)
# Sort TE_BED by chr then start to match GENOME_SIZES lexicographic order
sort -k1,1 -k2,2n "${TE_BED}" | bedtools complement -i - -g "${GENOME_SIZES}" > "${NON_TE_BED}"

N_NON_TE=$(awk '{sum+=$3-$2} END{print sum+0}' "${NON_TE_BED}")
echo "  Non-TE callable genome: ${N_NON_TE} bp" | tee -a "${LOG}"

# -------------------------------------------------------------------
# STEP 3: Select 10 representative GVCFs (3+3+2+2 across K=4 clusters)
# For each cluster, samples are sorted alphabetically and evenly spaced.
# -------------------------------------------------------------------
SELECTED_SAMPLES="${OUT_DIR}/selected_samples.tsv"
echo "STEP 3: Selecting 10 representative GVCFs ..." | tee -a "${LOG}"

K_COL="K${K_GLOBAL}"
> "${SELECTED_SAMPLES}"

for POP in Pop_1 Pop_2 Pop_3 Pop_4; do
    case "${POP}" in
        Pop_1|Pop_2) N_PICK=3 ;;
        *)            N_PICK=2 ;;
    esac

    # Extract sorted sample IDs for this cluster, then pick N_PICK evenly spaced
    awk -v kcol="${K_COL}" -v pop="${POP}" \
        'NR==1 { for(i=1;i<=NF;i++) if($i==kcol) col=i; next }
         $col==pop { print $1 }' \
        "${CLUSTER_ASSIGNMENTS}" | sort | \
    awk -v n="${N_PICK}" -v pop="${POP}" '
    { samples[++cnt]=$0 }
    END {
        if (cnt <= n) {
            for (i=1; i<=cnt; i++) print pop"\t"samples[i]
        } else {
            for (i=0; i<n; i++) {
                idx = (n == 1) ? 1 : int(i * (cnt-1) / (n-1) + 0.5) + 1
                if (idx > cnt) idx = cnt
                print pop"\t"samples[idx]
            }
        }
    }' >> "${SELECTED_SAMPLES}"

    N=$(awk -v p="${POP}" '$1==p' "${SELECTED_SAMPLES}" | wc -l)
    echo "  ${POP}: ${N} samples selected" | tee -a "${LOG}"
done

echo "  --- Selected samples ---" | tee -a "${LOG}"
cat "${SELECTED_SAMPLES}" | tee -a "${LOG}"
echo "  ------" | tee -a "${LOG}"

# -------------------------------------------------------------------
# STEP 4: Count callable non-TE sites per GVCF
# bcftools -R restricts to records overlapping non-TE regions.
# Reference blocks (INFO/END) may span TE/non-TE boundaries:
# awk expands them to 0-based BED intervals; bedtools intersect clips
# them at non-TE boundaries; final awk sums the clipped lengths.
# -------------------------------------------------------------------
L_VALUES="${OUT_DIR}/L_per_sample.tsv"
echo -e "population\tsample\tL_callable" > "${L_VALUES}"
echo "STEP 4: Counting callable non-TE sites per GVCF ..." | tee -a "${LOG}"

while IFS=$'\t' read -r POP SAMPLE; do
    GVCF="${GVCF_DIR}/${SAMPLE}.g.vcf.gz"
    if [ ! -s "${GVCF}" ]; then
        echo "  WARN: GVCF not found, skipping: ${GVCF}" | tee -a "${LOG}"
        continue
    fi

    echo -n "  ${POP}  ${SAMPLE}: " | tee -a "${LOG}"

    L=$(bcftools view \
            -R "${NON_TE_BED}" \
            -i "GQ>=${GQ_MIN}" \
            "${GVCF}" | \
        awk '!/^#/ {
            pos = $2 - 1
            if (match($8, /END=([0-9]+)/, m)) { end = m[1]+0 }
            else { end = pos + 1 }
            print $1"\t"pos"\t"end
        }' | \
        bedtools intersect -a stdin -b "${NON_TE_BED}" | \
        awk '{sum += $3 - $2} END {print sum+0}')

    echo "${L} bp" | tee -a "${LOG}"
    echo -e "${POP}\t${SAMPLE}\t${L}" >> "${L_VALUES}"
done < "${SELECTED_SAMPLES}"

# -------------------------------------------------------------------
# STEP 5: Compute mean and median L across all processed samples
# -------------------------------------------------------------------
echo "STEP 5: Summary ..." | tee -a "${LOG}"

awk -F'\t' '
NR==1 { next }
$3+0 > 0 {
    vals[++n] = $3+0
    sum += $3+0
}
END {
    if (n == 0) { print "ERROR: no valid L values computed"; exit 1 }
    mean = sum / n

    # Sort for median (bubble sort, n <= 10)
    for (i=1; i<=n; i++)
        for (j=i+1; j<=n; j++)
            if (vals[j] < vals[i]) { t=vals[i]; vals[i]=vals[j]; vals[j]=t }

    if (n % 2 == 0) median = (vals[n/2] + vals[n/2+1]) / 2
    else             median = vals[int(n/2)+1]

    printf "\n========= L estimation from %d GVCFs (GQ >= 20, non-TE) =========\n", n
    printf "  Mean   L = %d bp\n", int(mean)
    printf "  Median L = %d bp\n", int(median)
    printf "  Min    L = %d bp\n", vals[1]
    printf "  Max    L = %d bp\n", vals[n]
    printf "\n  --> To update 07.11-stairway_plot.slurm, set:\n"
    printf "      L_CALLABLE=%d\n", int(median)
    printf "=================================================================\n\n"
}' "${L_VALUES}" | tee -a "${LOG}"

echo "DONE $(date -Iseconds) — outputs in ${OUT_DIR}" | tee -a "${LOG}"
