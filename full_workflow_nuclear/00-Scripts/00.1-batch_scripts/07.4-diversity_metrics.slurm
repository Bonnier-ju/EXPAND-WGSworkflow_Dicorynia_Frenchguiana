#!/bin/bash
#SBATCH -p workq
#SBATCH --time=24:00:00
#SBATCH -J 07.4-diversity
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=20G
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module purge
module load bioinfo/Bcftools/1.9
module load bioinfo/VCFtools/0.1.16
module load bioinfo/PLINK/1.90b7

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
VCF_66="${BASE}/06-variant_filtering/workflow_test/06.6-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"
VCF_65="${BASE}/06-variant_filtering/workflow_test/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"

IN_DIR="${BASE}/07-population_structure/workflow_test/07.3-population_inputs"
OUT_DIR="${BASE}/07-population_structure/workflow_test/07.4-diversity_metrics"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR"

LOG_FILE="${LOG_DIR}/07.4-diversity_metrics.log"
SUMMARY_FILE="${OUT_DIR}/07.4-population_diversity_summary.tsv"
FST_SUMMARY="${OUT_DIR}/07.4-pairwise_fst_summary.tsv"

MAP_FILE="${IN_DIR}/sample_population_map.tsv"
KEEP_FILE="${IN_DIR}/keep.samples.txt"
CLUSTER_FILE="${IN_DIR}/plink.cluster.txt"
POPS_LIST="${IN_DIR}/populations.list"
SAMPLES_DIR="${IN_DIR}/population_samples"

if [ -s "$VCF_66" ]; then
    INPUT_VCF="$VCF_66"
    INPUT_STEP="06.6"
else
    INPUT_VCF="$VCF_65"
    INPUT_STEP="06.5"
fi

echo "START $(date -Iseconds) 07.4 diversity metrics" > "$LOG_FILE"
echo "Input VCF: $INPUT_VCF (from $INPUT_STEP)" | tee -a "$LOG_FILE"
echo "Input population map: $MAP_FILE" | tee -a "$LOG_FILE"

for f in "$INPUT_VCF" "$MAP_FILE" "$KEEP_FILE" "$CLUSTER_FILE" "$POPS_LIST"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: required file missing: $f" | tee -a "$LOG_FILE"
        exit 1
    fi
done

# vcftools 0.1.16 cannot read VCFv4.3. If needed, downgrade header to VCFv4.2.
VCFTOOLS_VCF="$INPUT_VCF"
VCF_FMT=$(bcftools view -h "$INPUT_VCF" | awk -F= '/^##fileformat=/{print $2; exit}')
if [ "$VCF_FMT" = "VCFv4.3" ]; then
    echo "INFO: Input is VCFv4.3, creating VCFv4.2-compatible copy for vcftools." | tee -a "$LOG_FILE"
    VCFTOOLS_VCF="${TMP_DIR}/input_for_vcftools.vcf.gz"
    VCFTOOLS_HDR="${TMP_DIR}/input_for_vcftools.header.txt"
    bcftools view -h "$INPUT_VCF" | sed 's/^##fileformat=VCFv4\.3$/##fileformat=VCFv4.2/' > "$VCFTOOLS_HDR"
    bcftools reheader -h "$VCFTOOLS_HDR" -o "$VCFTOOLS_VCF" "$INPUT_VCF"
    bcftools index -f -t "$VCFTOOLS_VCF"
    echo "INFO: VCFv4.2-compatible file ready: $VCFTOOLS_VCF" | tee -a "$LOG_FILE"
fi

# 1) Per-population metrics (Ho, He, FIS, selfing proxy, pi, TajimaD)
echo -e "population\tn_individuals\tmean_Ho\tmean_He\tmean_Fis\tselfing_proxy\tmean_pi\tmean_tajimaD" > "$SUMMARY_FILE"

while read -r pop; do
    [ -z "$pop" ] && continue
    safe=$(echo "$pop" | sed 's/[^A-Za-z0-9._-]/_/g')
    keep_pop="${SAMPLES_DIR}/${safe}.keep.txt"
    if [ ! -s "$keep_pop" ]; then
        echo "WARN: skip pop ${pop}, missing keep file" | tee -a "$LOG_FILE"
        continue
    fi

    n_ind=$(wc -l < "$keep_pop")
    prefix="${OUT_DIR}/pop_${safe}"

    vcftools --gzvcf "$VCFTOOLS_VCF" --keep "$keep_pop" --het --out "$prefix" >> "$LOG_FILE" 2>&1 || true
    vcftools --gzvcf "$VCFTOOLS_VCF" --keep "$keep_pop" --site-pi --out "$prefix" >> "$LOG_FILE" 2>&1 || true
    vcftools --gzvcf "$VCFTOOLS_VCF" --keep "$keep_pop" --TajimaD 10000 --out "$prefix" >> "$LOG_FILE" 2>&1 || true

    het_file="${prefix}.het"
    pi_file="${prefix}.sites.pi"
    taj_file="${prefix}.Tajima.D"

    mean_Ho="NA"
    mean_He="NA"
    mean_Fis="NA"
    selfing="NA"
    mean_pi="NA"
    mean_taj="NA"

    if [ -s "$het_file" ]; then
        mean_Ho=$(awk 'NR>1 && $3>0{ho=($3-$1)/$3; s+=ho; n++} END{if(n>0) printf "%.6f", s/n; else print "NA"}' "$het_file")
        mean_He=$(awk 'NR>1 && $3>0{he=($3-$2)/$3; s+=he; n++} END{if(n>0) printf "%.6f", s/n; else print "NA"}' "$het_file")
        mean_Fis=$(awk 'NR>1{s+=$4; n++} END{if(n>0) printf "%.6f", s/n; else print "NA"}' "$het_file")
        selfing=$(awk 'NR>1{s+=$4; n++} END{if(n>0){f=s/n; den=1+f; if(den!=0) printf "%.6f", (2*f)/den; else print "NA"} else print "NA"}' "$het_file")
    fi

    if [ -s "$pi_file" ]; then
        mean_pi=$(awk 'NR>1 && $3!="nan"{s+=$3; n++} END{if(n>0) printf "%.6f", s/n; else print "NA"}' "$pi_file")
    fi

    if [ -s "$taj_file" ]; then
        mean_taj=$(awk 'NR>1 && $4!="nan"{s+=$4; n++} END{if(n>0) printf "%.6f", s/n; else print "NA"}' "$taj_file")
    fi

    echo -e "${pop}\t${n_ind}\t${mean_Ho}\t${mean_He}\t${mean_Fis}\t${selfing}\t${mean_pi}\t${mean_taj}" >> "$SUMMARY_FILE"
done < "$POPS_LIST"

# 2) Allelic richness proxy + private polymorphism proxy from plink strat freq
plink_prefix="${OUT_DIR}/for_freq"
plink --vcf "$INPUT_VCF" \
    --allow-extra-chr \
    --double-id \
    --keep "$KEEP_FILE" \
    --make-bed \
    --out "$plink_prefix" >> "$LOG_FILE" 2>&1

plink --bfile "$plink_prefix" \
    --allow-extra-chr \
    --freq \
    --within "$CLUSTER_FILE" \
    --out "${OUT_DIR}/pop_freq" >> "$LOG_FILE" 2>&1

FRQ_STRAT="${OUT_DIR}/pop_freq.frq.strat"
AR_OUT="${OUT_DIR}/07.4-allelic_richness_proxy.tsv"
PRIV_OUT="${OUT_DIR}/07.4-private_polymorphism_proxy.tsv"

if [ -s "$FRQ_STRAT" ]; then
    awk '
    NR==1{
      for(i=1;i<=NF;i++){
        if($i=="SNP") s=i;
        if($i=="CLST") c=i;
        if($i=="MAF") m=i;
      }
      next
    }
    {
      snp=$s; pop=$c; maf=$m+0;
      total[pop]++;
      if(maf>0){poly[pop]++; alleles[pop]+=2; snp_poly[snp]++; last_pop[snp]=pop}
      else {alleles[pop]+=1}
    }
    END{
      print "population\tloci_total\tpolymorphic_loci\tallelic_richness_proxy_mean_alleles_per_locus";
      for(p in total){
        ar=(total[p]>0)?alleles[p]/total[p]:"NA";
        printf "%s\t%d\t%d\t%.6f\n", p, total[p], poly[p]+0, ar;
      }
    }' "$FRQ_STRAT" > "$AR_OUT"

    awk '
    NR==1{
      for(i=1;i<=NF;i++){
        if($i=="SNP") s=i;
        if($i=="CLST") c=i;
        if($i=="MAF") m=i;
      }
      next
    }
    {
      snp=$s; pop=$c; maf=$m+0;
      if(maf>0){snp_poly[snp]++; last_pop[snp]=pop}
    }
    END{
      for(snp in snp_poly){
        if(snp_poly[snp]==1){priv[last_pop[snp]]++}
      }
      print "population\tprivate_polymorphic_loci_proxy";
      for(p in priv) printf "%s\t%d\n", p, priv[p];
    }' "$FRQ_STRAT" > "$PRIV_OUT"
fi

# 3) Pairwise FST (vcftools)
echo -e "pop1\tpop2\tmean_fst" > "$FST_SUMMARY"
mapfile -t pops < "$POPS_LIST"
for ((i=0; i<${#pops[@]}; i++)); do
    for ((j=i+1; j<${#pops[@]}; j++)); do
        p1="${pops[$i]}"
        p2="${pops[$j]}"
        s1=$(echo "$p1" | sed 's/[^A-Za-z0-9._-]/_/g')
        s2=$(echo "$p2" | sed 's/[^A-Za-z0-9._-]/_/g')
        k1="${SAMPLES_DIR}/${s1}.samples.txt"
        k2="${SAMPLES_DIR}/${s2}.samples.txt"
        outp="${OUT_DIR}/fst_${s1}_vs_${s2}"
        vcftools --gzvcf "$VCFTOOLS_VCF" \
            --weir-fst-pop "$k1" \
            --weir-fst-pop "$k2" \
            --out "$outp" >> "$LOG_FILE" 2>&1 || true
        fst_file="${outp}.weir.fst"
        mean_fst="NA"
        if [ -s "$fst_file" ]; then
            mean_fst=$(awk '$3!="nan"{s+=$3; n++} END{if(n>0) printf "%.6f", s/n; else print "NA"}' "$fst_file")
        fi
        echo -e "${p1}\t${p2}\t${mean_fst}" >> "$FST_SUMMARY"
    done
done

echo "DONE $(date -Iseconds) 07.4 diversity metrics completed" | tee -a "$LOG_FILE"
echo "Main outputs:"
echo "  - $SUMMARY_FILE"
echo "  - $AR_OUT"
echo "  - $PRIV_OUT"
echo "  - $FST_SUMMARY"
