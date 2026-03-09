#!/bin/bash
#SBATCH -p workq
#SBATCH --time=01:00:00
#SBATCH -J cp-06.6-alignment_prep
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/trimAl/1.4.1

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# Pour tester sur le jeu de données test, changer WORKFLOW :
#   WORKFLOW="workflow_test_9ind"
# -------------------------------------------------------------------
WORKFLOW="workflow_full"

BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic"
INPUT_FASTA="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/consensus_sequences/all_samples_consensus.fasta"
OUTPUT_DIR="${BASE_DIR}/07-phylogeny_phylogeography_structure/${WORKFLOW}/alignments_phylo"
LOG_FILE="${OUTPUT_DIR}/alignment_prep_cp.log"

OUTPUT_FASTA="${OUTPUT_DIR}/alignment_trimmed.fasta"
STATS_FILE="${OUTPUT_DIR}/trimal_stats.txt"

# Seuil trimAl : conserver les colonnes présentes dans >= 80% des séquences
# (élimine les colonnes avec >20% de N ou de gaps)
GT_THRESHOLD=0.8

mkdir -p "$OUTPUT_DIR"

if [ ! -s "$INPUT_FASTA" ]; then
    echo "ERROR: input FASTA not found or empty: $INPUT_FASTA" >&2; exit 1
fi
if [ -s "$OUTPUT_FASTA" ]; then
    echo "SKIP: output already exists: $OUTPUT_FASTA"; exit 0
fi

echo "START $(date -Iseconds) alignment_prep cp" > "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: Stats de l'alignement d'entrée
# -------------------------------------------------------------------
echo "INFO: input alignment stats" >> "$LOG_FILE"
n_seq_in=$(grep -c "^>" "$INPUT_FASTA")
seq_len_in=$(awk '/^>/{next} {len+=length($0)} END{print len/NR}' "$INPUT_FASTA" | awk '{printf "%d\n", $1}')
echo "INFO: sequences=${n_seq_in} mean_length_before=${seq_len_in}" >> "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Nettoyage de l'alignement avec trimAl
#   -gt 0.8 : conserver colonnes présentes dans >= 80% des séquences
#             (supprime colonnes avec >20% de gaps/N)
#   -fasta   : format de sortie FASTA
#   -colnumbering : afficher les positions conservées
# Note : les séquences consensus sur le même génome de référence sont
#        déjà alignées en coordonnées → trimAl retire uniquement les
#        colonnes trop incomplètes (beaucoup de N)
# -------------------------------------------------------------------
echo "INFO: running trimAl -gt ${GT_THRESHOLD}" >> "$LOG_FILE"

trimal \
    -in  "$INPUT_FASTA" \
    -out "$OUTPUT_FASTA" \
    -gt  "$GT_THRESHOLD" \
    -fasta 2>> "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 4: Stats de l'alignement de sortie
# -------------------------------------------------------------------
if [ -s "$OUTPUT_FASTA" ]; then
    n_seq_out=$(grep -c "^>" "$OUTPUT_FASTA")
    seq_len_out=$(awk '/^>/{next} {len+=length($0)} END{print len/NR}' "$OUTPUT_FASTA" | awk '{printf "%d\n", $1}')

    # Stats détaillées avec trimAl
    trimal \
        -in  "$OUTPUT_FASTA" \
        -sident 2>/dev/null > "$STATS_FILE" || true

    echo "DONE $(date -Iseconds) sequences=${n_seq_out} positions_before=${seq_len_in} positions_after=${seq_len_out} output=${OUTPUT_FASTA}" | tee -a "$LOG_FILE"

    if [ "$n_seq_out" -ne "$n_seq_in" ]; then
        echo "WARNING: sequence count changed (in=${n_seq_in} out=${n_seq_out}) — check for all-gap sequences" | tee -a "$LOG_FILE"
    fi
else
    echo "ERROR: trimAl output is empty" | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# Description
# Step 1: Vérification du multi-FASTA d'entrée (06.5 consensus_sequences)
# Step 2: trimAl -gt 0.8 → supprime colonnes avec >20% de gaps/N
#         Les séquences sont déjà pseudo-alignées sur le même génome
#         de référence → pas besoin de MAFFT
# Step 3: Le fichier de sortie est prêt pour IQ-TREE2 (07.1)
# Sorties :
#   alignment_trimmed.fasta : alignement nettoyé pour phylogénie
#   trimal_stats.txt        : identité/conservation des colonnes
# Outgroup : individus Cameroun_Benin déjà présents dans l'alignement,
#            seront déclarés outgroup dans IQ-TREE (07.1)
# -------------------------------------------------------------------
