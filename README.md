# EXPAND-WGSworkflow_Dicorynia_Frenchguiana
This repository contains reproducible pipelines to analyze 10X whole-genome sequencing data from 180 individuals sampled across ~20 sites in French Guiana. Using a reference genome and a GATK-based workflow on HPC cluster, the project investigates genetic diversity, phylogeography, and genotype–environment associations.
Nuclear WGS workflow for *Dicorynia* on HPC (Slurm), from raw FASTQ QC to variant discovery.  
Pipeline includes FastQC/MultiQC, fastp cleaning, bwa-mem2 alignment, samtools/qualimap coverage QC, and downstream variant calling/filtering with GATK.

