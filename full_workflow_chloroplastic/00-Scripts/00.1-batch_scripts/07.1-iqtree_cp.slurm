#!/bin/bash
#SBATCH -p workq
#SBATCH --time=06:00:00
#SBATCH -J cp-07.1-iqtree
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load compilers/gcc/12.2.0
module load bioinfo/IQ-TREE/3.0.1
module load statistics/R/4.3.0

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# To run on the test dataset, change WORKFLOW:
#   WORKFLOW="workflow_test_9ind"
# -------------------------------------------------------------------
WORKFLOW="workflow_full"

BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic"
INPUT_ALN="${BASE_DIR}/07-phylogeny_phylogeography_structure/${WORKFLOW}/alignments_phylo/alignment_trimmed.fasta"
RESULTS_DIR="${BASE_DIR}/07-phylogeny_phylogeography_structure/${WORKFLOW}/07.1-iqtree_phylogeny"
LOG_FILE="${RESULTS_DIR}/iqtree_cp.log"
PREFIX="${RESULTS_DIR}/cp_phylo"
METADATA="${BASE_DIR}/metadata/sample_sheet_chloroplastic_complete_filtered.csv"
COLORS="${BASE_DIR}/metadata/sites_couleurs.csv"
R_SCRIPT="${BASE_DIR}/00-Scripts/00.2-R_scripts/07.1-iqtree_viz_cp.R"

# Outgroup: 4 Cameroun/Benin individuals (outgroup declared in sample sheet)
# Used to root the tree — geographically and genetically distant from French Guiana
OUTGROUP="DBR_AGN,DBR_AGP,DBR_AGQ,DBR_AGR"

THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p "$RESULTS_DIR"

if [ ! -s "$INPUT_ALN" ]; then
    echo "ERROR: input alignment not found or empty: $INPUT_ALN" >&2; exit 1
fi
if [ -s "${PREFIX}.treefile" ]; then
    echo "SKIP: output already exists: ${PREFIX}.treefile"; exit 0
fi

echo "START $(date -Iseconds) iqtree cp" > "$LOG_FILE"
echo "INFO: alignment=${INPUT_ALN}" >> "$LOG_FILE"
echo "INFO: outgroup=${OUTGROUP}" >> "$LOG_FILE"
echo "INFO: threads=${THREADS}" >> "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: Model selection + ML tree + bootstrap
#   -s   : input alignment (FASTA)
#   -m MFP : ModelFinder Plus — selects best substitution model
#            (BIC criterion) and runs tree inference with that model
#   -B 1000 : ultrafast bootstrap (UFBoot2) — 1000 replicates
#             values >= 95 considered well-supported
#   -alrt 1000 : SH-aLRT test — 1000 replicates
#                complementary support metric (>= 80% considered good)
#   -o   : outgroup sample IDs (comma-separated, no spaces)
#   -T   : number of threads
#   --prefix : output file prefix
#   --redo   : overwrite previous run if partial files exist
# -------------------------------------------------------------------
echo "INFO: running IQ-TREE3 ModelFinder + ML + UFBoot" >> "$LOG_FILE"

iqtree3 \
    -s "$INPUT_ALN" \
    -m MFP \
    -B 1000 \
    -alrt 1000 \
    -o "$OUTGROUP" \
    -T "$THREADS" \
    --prefix "$PREFIX" \
    --redo \
    2>> "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Validate outputs
# -------------------------------------------------------------------
if [ -s "${PREFIX}.treefile" ]; then
    best_model=$(grep "^Best-fit model" "${PREFIX}.iqtree" 2>/dev/null | awk '{print $NF}' || echo "unknown")
    log_l=$(grep "^Log-likelihood of the tree" "${PREFIX}.iqtree" 2>/dev/null | awk '{print $NF}' || echo "unknown")
    n_leaves=$(grep -c "^>" "$INPUT_ALN" || echo "unknown")
    echo "DONE $(date -Iseconds) model=${best_model} logL=${log_l} leaves=${n_leaves} output=${PREFIX}.treefile" | tee -a "$LOG_FILE"
else
    echo "ERROR: IQ-TREE output treefile not found" | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# STEP 4: Tree visualization with R (ggtree)
# -------------------------------------------------------------------
echo "INFO: running R visualization script" >> "$LOG_FILE"

Rscript --vanilla "$R_SCRIPT" \
    "${PREFIX}.treefile" \
    "$METADATA" \
    "$COLORS" \
    "$RESULTS_DIR" \
    >> "$LOG_FILE" 2>&1

if ls "${RESULTS_DIR}/cp_phylo_all_rect.png" &>/dev/null; then
    echo "DONE $(date -Iseconds) tree plots saved to ${RESULTS_DIR}" | tee -a "$LOG_FILE"
else
    echo "ERROR: R visualization failed — check log" | tee -a "$LOG_FILE"; exit 1
fi


# -------------------------------------------------------------------
# Description
# Maximum likelihood phylogenetic inference of the chloroplast alignment.
# Step 1: ModelFinder (MFP) selects the best-fit substitution model
#         using BIC — avoids under/over-parameterized models.
# Step 2: ML tree inference with the selected model.
# Step 3: UFBoot2 (1000 rep.) + SH-aLRT (1000 rep.) for node support.
#         UFBoot >= 95 and SH-aLRT >= 80 indicate well-supported nodes.
# Step 4: Tree rooted with 4 Cameroun/Benin outgroup individuals.
# All outputs in 07.1-iqtree_phylogeny/:
#   cp_phylo.treefile              : ML tree with bootstrap support (Newick)
#   cp_phylo.iqtree                : full run report (model, logL, parameters)
#   cp_phylo.log                   : IQ-TREE run log
#   cp_phylo.mldist                : pairwise ML distance matrix
#   cp_phylo_all_rect.png          : all samples, rectangular
#   cp_phylo_all_rect_sqrt.png     : all samples, rectangular, sqrt branch lengths
#   cp_phylo_all_circ.png          : all samples, circular
#   cp_phylo_all_circ_sqrt.png     : all samples, circular, sqrt branch lengths
#   cp_phylo_fg_rect.png           : French Guiana only, rectangular
#   cp_phylo_sites_rect.png        : site-level, rectangular (true branch lengths)
#   cp_phylo_sites_rect_clado.png  : site-level, rectangular, cladogram (topology only)
#   cp_phylo_sites_circ.png        : site-level, circular (true branch lengths)
#   cp_phylo_sites_circ_clado.png  : site-level, circular, cladogram (topology only)
#   cp_phylo_node_support.tsv      : node support table (alrt, ufboot)
# Haploid chloroplast: no ploidy correction needed
# -------------------------------------------------------------------
