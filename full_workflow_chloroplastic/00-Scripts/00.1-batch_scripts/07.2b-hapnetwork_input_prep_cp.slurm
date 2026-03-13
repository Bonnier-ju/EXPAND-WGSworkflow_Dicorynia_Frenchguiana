#!/bin/bash
#SBATCH -p workq
#SBATCH --time=00:15:00
#SBATCH -J cp-07.2b-hapnet-input
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# =============================================================================
# STEP 07.2b — Prepare input files for HapNetworkView
# Produces:
#   - alignment_trimmed.phy     : PHYLIP sequential format (all 164 sequences)
#   - alignment_trimmed_fg.phy  : PHYLIP sequential format (159 FG samples only,
#                                 no Cameroun_Benin, no Dicorynia paranensis)
#   - group_file_all.txt        : sample_id -> site (all samples, tab-separated)
#   - group_file_fg.txt         : sample_id -> site (French Guiana only)
# =============================================================================

# --- Paths -------------------------------------------------------------------
WORKFLOW="workflow_full"
BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic"

ALN="${BASE_DIR}/07-phylogeny_phylogeography_structure/${WORKFLOW}/alignments_phylo/alignment_trimmed.fasta"
META="${BASE_DIR}/metadata/sample_sheet_chloroplastic_complete_filtered.csv"
OUTPUT_DIR="${BASE_DIR}/07-phylogeny_phylogeography_structure/${WORKFLOW}/07.2-haplotype_networks"

mkdir -p "${OUTPUT_DIR}"

LOG="${OUTPUT_DIR}/hapnetwork_input_prep_cp.log"
echo "START $(date -Iseconds) HapNetworkView input preparation" | tee "$LOG"

# --- Step 1: FASTA -> PHYLIP sequential format --------------------------------
echo "INFO: converting FASTA to PHYLIP..." | tee -a "$LOG"

PHY="${OUTPUT_DIR}/alignment_trimmed.phy"

python3 - <<'PYEOF'
import sys

fasta_file = "/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/07-phylogeny_phylogeography_structure/workflow_full/alignments_phylo/alignment_trimmed.fasta"
phy_file   = "/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/07-phylogeny_phylogeography_structure/workflow_full/07.2-haplotype_networks/alignment_trimmed.phy"

# Parse FASTA
seqs = {}
order = []
name = None
with open(fasta_file) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            order.append(name)
            seqs[name] = ""
        else:
            seqs[name] += line

n_seq = len(order)
seq_len = len(seqs[order[0]])

# Validate all sequences have same length
for s in order:
    if len(seqs[s]) != seq_len:
        sys.exit(f"ERROR: sequence {s} has length {len(seqs[s])}, expected {seq_len}")

# Write PHYLIP sequential
with open(phy_file, "w") as out:
    out.write(f" {n_seq} {seq_len}\n")
    for s in order:
        # Name truncated/padded to 10 characters
        name_field = s[:10].ljust(10)
        out.write(f"{name_field} {seqs[s]}\n")

print(f"INFO: PHYLIP written: {n_seq} sequences x {seq_len} bp -> {phy_file}")
PYEOF

echo "INFO: PHYLIP file created: ${PHY}" | tee -a "$LOG"

# --- Step 2: Group file — all samples ----------------------------------------
echo "INFO: generating group file (all samples)..." | tee -a "$LOG"

GROUP_ALL="${OUTPUT_DIR}/group_file_all.txt"

# Header-free tab-separated: sample_id <TAB> site
# Treemutation -> Angela, Herbier -> Dicorynia paranensis
tail -n +2 "${META}" | awk -F',' '
{
    sample = $1
    project = $3
    site    = $4
    if (project == "Treemutation") site = "Angela"
    if (site == "Herbier") site = "Dicorynia_paranensis"
    print sample "\t" site
}' > "${GROUP_ALL}"

n_all=$(wc -l < "${GROUP_ALL}")
echo "INFO: group_file_all.txt: ${n_all} samples" | tee -a "$LOG"

# --- Step 3: Group file — French Guiana only ---------------------------------
echo "INFO: generating group file (French Guiana only)..." | tee -a "$LOG"

GROUP_FG="${OUTPUT_DIR}/group_file_fg.txt"

# Exclude Cameroun_Benin and Dicorynia paranensis
tail -n +2 "${META}" | awk -F',' '
{
    sample = $1
    project = $3
    site    = $4
    if (project == "Treemutation") site = "Angela"
    if (site == "Herbier") site = "Dicorynia_paranensis"
    if (site != "Cameroun_Benin" && site != "Dicorynia_paranensis")
        print sample "\t" site
}' > "${GROUP_FG}"

n_fg=$(wc -l < "${GROUP_FG}")
echo "INFO: group_file_fg.txt: ${n_fg} samples" | tee -a "$LOG"

# --- Step 4: PHYLIP — French Guiana only (no Cameroun_Benin, no D.paranensis) -
echo "INFO: generating PHYLIP (French Guiana only)..." | tee -a "$LOG"

PHY_FG="${OUTPUT_DIR}/alignment_trimmed_fg.phy"

# Build list of FG sample IDs from group_file_fg.txt
FG_SAMPLES=$(cut -f1 "${GROUP_FG}")

python3 - <<PYEOF
import sys

fasta_file = "/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/07-phylogeny_phylogeography_structure/workflow_full/alignments_phylo/alignment_trimmed.fasta"
phy_file   = "${PHY_FG}"
fg_samples = set("""${FG_SAMPLES}""".split())

# Parse FASTA — keep only FG samples
seqs  = {}
order = []
name  = None
with open(fasta_file) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
        else:
            if name in fg_samples:
                seqs.setdefault(name, "")
                seqs[name] += line
                if name not in order:
                    order.append(name)

n_seq   = len(order)
seq_len = len(seqs[order[0]])

with open(phy_file, "w") as out:
    out.write(f" {n_seq} {seq_len}\n")
    for s in order:
        out.write(f"{s[:10].ljust(10)} {seqs[s]}\n")

print(f"INFO: PHYLIP FG written: {n_seq} sequences x {seq_len} bp -> {phy_file}")
PYEOF

echo "INFO: PHYLIP FG file created: ${PHY_FG}" | tee -a "$LOG"

# --- Check ------------------------------------------------------------------
echo "INFO: output files:" | tee -a "$LOG"
ls -lh "${OUTPUT_DIR}/alignment_trimmed.phy" \
        "${OUTPUT_DIR}/alignment_trimmed_fg.phy" \
        "${GROUP_ALL}" \
        "${GROUP_FG}" | tee -a "$LOG"

echo "DONE $(date -Iseconds) output=${OUTPUT_DIR}" | tee -a "$LOG"
