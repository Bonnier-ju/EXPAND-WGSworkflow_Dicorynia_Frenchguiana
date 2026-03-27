#!/usr/bin/env Rscript
# =============================================================================
# Script: 07.7-samova_input_prep_cp.R
# Purpose: Generate SAMOVA 2.0 input files from cpDNA alignment data
#   - dicorynia_cpDNA.geo : geographic coordinates file
#   - dicorynia_cpDNA.arp : Arlequin 3.5 DNA sequences file
#
# SAMOVA 2.0 reference:
#   Dupanloup et al. (2002) Mol. Ecol. 11:2571-2581
#   Software: https://cmpg.unibe.ch/software/samova2/
#
# Exclusions: Angela site (no GPS / uncertain assignment), outgroups
#   (Cameroun_Benin, Herbier)
#
# Input files:
#   args[1] : alignment FASTA (trimmed, one sequence per individual)
#   args[2] : sample metadata CSV (sample_id, site, ...)
#   args[3] : GPS coordinates CSV (Sites, lat, long)
#   args[4] : output directory
#
# Output files (saved to output_dir):
#   dicorynia_cpDNA.geo   : tab-separated GPS file for SAMOVA
#   dicorynia_cpDNA.arp   : Arlequin 3.5 format DNA data
#   samova_site_summary.tsv : per-site sample counts (QC)
#
# .geo format (SAMOVA 2.0):
#   N<tab>"SiteName"<tab>longitude<tab>latitude<tab>1
#   where N = sequential site number (1-based), last column = number of groups (always 1)
#
# .arp format (Arlequin 3.5 DNA):
#   Arlequin profile block followed by Data block
#   Each SampleData block: SampleName, SampleSize, haplotype IDs + frequencies
#   HaplotypeDefinition block: unique sequences (one per unique haplotype across all pops)
#
# Notes:
#   - Sequences are collapsed to unique haplotypes across all included sites
#   - Haplotype names: Hap_001, Hap_002, ... (sorted by decreasing frequency)
#   - Seed: 42 (not used here as no randomisation, but documented for reproducibility)
# =============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ape)       # read.FASTA
  library(stringr)
})

# ---------------------------------------------------------------------------
# STEP 0: Parse and validate arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 07.7-samova_input_prep_cp.R <alignment.fasta> <metadata.csv> <geoloc.csv> <output_dir>")
}

fasta_file  <- args[1]
meta_file   <- args[2]
geoloc_file <- args[3]
output_dir  <- args[4]

for (f in c(fasta_file, meta_file, geoloc_file)) {
  if (!file.exists(f) || file.size(f) == 0) {
    stop(paste("ERROR: input file missing or empty:", f))
  }
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Sites to exclude from SAMOVA (Angela = uncertain assignment; outgroups)
EXCLUDE_SITES <- c("Angela", "Cameroun_Benin", "Herbier")

# ---------------------------------------------------------------------------
# STEP 1: Load metadata and filter to French Guiana sites only
# ---------------------------------------------------------------------------
cat("Loading metadata...\n")
meta <- read_csv(meta_file, show_col_types = FALSE) %>%
  select(sample_id, site) %>%
  filter(!site %in% EXCLUDE_SITES)

cat(sprintf("  Retained %d samples across %d sites\n",
            nrow(meta), n_distinct(meta$site)))

# ---------------------------------------------------------------------------
# STEP 2: Load GPS coordinates
# ---------------------------------------------------------------------------
cat("Loading GPS coordinates...\n")
geoloc <- read_csv(geoloc_file, show_col_types = FALSE) %>%
  rename(site = Sites) %>%
  filter(site %in% unique(meta$site))

# Verify all sites have GPS
missing_gps <- setdiff(unique(meta$site), geoloc$site)
if (length(missing_gps) > 0) {
  stop(paste("ERROR: No GPS for sites:", paste(missing_gps, collapse = ", ")))
}

# Order sites alphabetically for reproducibility
geoloc <- geoloc %>% arrange(site)
sites_ordered <- geoloc$site
n_sites <- length(sites_ordered)
cat(sprintf("  %d sites with GPS coordinates\n", n_sites))

# ---------------------------------------------------------------------------
# STEP 3: Load alignment and filter to retained samples
# ---------------------------------------------------------------------------
cat("Loading alignment FASTA...\n")
aln_raw <- read.FASTA(fasta_file, type = "DNA")

# FASTA headers may include description after space — extract sample IDs only
aln_names <- names(aln_raw)
sample_ids_aln <- str_extract(aln_names, "^[^ ]+")
names(aln_raw) <- sample_ids_aln

# Keep only samples present in filtered metadata
keep_samples <- meta$sample_id[meta$sample_id %in% sample_ids_aln]
missing_aln <- setdiff(meta$sample_id, sample_ids_aln)
if (length(missing_aln) > 0) {
  warning(paste("Samples in metadata but not in FASTA (skipped):",
                paste(missing_aln, collapse = ", ")))
}

aln <- aln_raw[keep_samples]
cat(sprintf("  %d sequences retained\n", length(aln)))

# ---------------------------------------------------------------------------
# STEP 4: Convert alignment to character matrix and identify unique haplotypes
# ---------------------------------------------------------------------------
cat("Identifying unique haplotypes...\n")

# Convert ape DNAbin to character matrix (rows = samples, cols = sites)
aln_mat <- as.character(as.matrix(aln))
rownames(aln_mat) <- names(aln)

seq_len <- ncol(aln_mat)
cat(sprintf("  Alignment length: %d bp\n", seq_len))

# Collapse each row to a single string for comparison
seq_strings <- apply(aln_mat, 1, paste, collapse = "")

# Build haplotype table: unique sequences across all retained samples
unique_seqs <- unique(seq_strings)

# Sort haplotypes by decreasing global frequency
hap_freq <- sort(table(seq_strings), decreasing = TRUE)
unique_seqs_ordered <- names(hap_freq)
n_haps <- length(unique_seqs_ordered)
cat(sprintf("  %d unique haplotypes identified\n", n_haps))

# Assign haplotype names (zero-padded)
hap_names <- sprintf("Hap_%03d", seq_len(n_haps))
names(hap_names) <- unique_seqs_ordered  # key = sequence string, value = Hap_XXX

# Map each sample to its haplotype name
sample_hap <- hap_names[seq_strings]
names(sample_hap) <- names(seq_strings)

# ---------------------------------------------------------------------------
# STEP 5: Build per-site haplotype count tables
# ---------------------------------------------------------------------------
cat("Building per-site haplotype tables...\n")

meta_hap <- meta %>%
  filter(sample_id %in% names(sample_hap)) %>%
  mutate(haplotype = sample_hap[sample_id])

# Per-site summary for QC output
site_summary <- meta_hap %>%
  group_by(site) %>%
  summarise(
    n_samples    = n(),
    n_haplotypes = n_distinct(haplotype),
    .groups = "drop"
  ) %>%
  arrange(site)

write_tsv(site_summary, file.path(output_dir, "samova_site_summary.tsv"))
cat("  [OK] samova_site_summary.tsv written\n")

# ---------------------------------------------------------------------------
# STEP 6: Write .geo file (SAMOVA geographic coordinates)
# ---------------------------------------------------------------------------
# Format: N<TAB>"SiteName"<TAB>longitude<TAB>latitude<TAB>1
# (longitude before latitude per SAMOVA 2.0 documentation)
cat("Writing .geo file...\n")

geo_lines <- character(n_sites)
for (i in seq_len(n_sites)) {
  s <- sites_ordered[i]
  row <- geoloc %>% filter(site == s)
  # SAMOVA expects: index  "name"  lon  lat  1
  geo_lines[i] <- paste(i, sprintf('"%s"', s), row$long, row$lat, 1, sep = "\t")
}

geo_path <- file.path(output_dir, "dicorynia_cpDNA.geo")
writeLines(geo_lines, geo_path)
cat(sprintf("  [OK] dicorynia_cpDNA.geo written (%d sites)\n", n_sites))

# ---------------------------------------------------------------------------
# STEP 7: Write .arp file (Arlequin 3.5 format)
# ---------------------------------------------------------------------------
# Structure:
#   [Profile] block — declares data type, nb haplotypes, etc.
#   [Data]
#     [HaplotypeDefinition] — lists each unique haplotype sequence
#     [Samples] — one SampleData block per population
# ---------------------------------------------------------------------------
cat("Writing .arp file...\n")

arp_lines <- character(0)

# --- Profile block ---
arp_lines <- c(arp_lines,
  "[Profile]",
  sprintf('\tTitle="Dicorynia guianensis cpDNA - %s"', format(Sys.Date(), "%Y-%m-%d")),
  sprintf("\tNbSamples=%d", n_sites),
  "\tDataType=DNA",
  "\tGenotypicData=0",   # 0 = haploid
  "\tGameticPhase=1",
  '\tMissingData="?"',
  "\tLocusSeparator=NONE",
  ""
)

# --- Data block header ---
arp_lines <- c(arp_lines, "[Data]", "")

# --- HaplotypeDefinition block ---
arp_lines <- c(arp_lines, "\t[[HaplotypeDefinition]]")
arp_lines <- c(arp_lines, "\t\tHaplListName=\"cpDNA haplotypes\"")
arp_lines <- c(arp_lines, "\t\tHaplList={")

for (i in seq_len(n_haps)) {
  hap_seq_str <- unique_seqs_ordered[i]
  hap_name    <- hap_names[i]
  arp_lines <- c(arp_lines, sprintf("\t\t\t%s\t%s", hap_name, hap_seq_str))
}
arp_lines <- c(arp_lines, "\t\t}", "")

# --- Samples block ---
arp_lines <- c(arp_lines, "\t[[Samples]]")
arp_lines <- c(arp_lines, "")

for (s in sites_ordered) {
  site_meta <- meta_hap %>% filter(site == s)
  n_ind <- nrow(site_meta)

  # Count occurrences of each haplotype in this site
  hap_counts <- site_meta %>%
    count(haplotype, name = "freq") %>%
    arrange(haplotype)

  arp_lines <- c(arp_lines,
    sprintf('\t\tSampleName="%s"', s),
    sprintf("\t\tSampleSize=%d", n_ind),
    "\t\tSampleData= {"
  )

  for (j in seq_len(nrow(hap_counts))) {
    arp_lines <- c(arp_lines,
      sprintf("\t\t\t%s\t%d", hap_counts$haplotype[j], hap_counts$freq[j])
    )
  }

  arp_lines <- c(arp_lines, "\t\t}", "")
}

arp_path <- file.path(output_dir, "dicorynia_cpDNA.arp")
writeLines(arp_lines, arp_path)
cat(sprintf("  [OK] dicorynia_cpDNA.arp written (%d populations, %d haplotypes)\n",
            n_sites, n_haps))

# ---------------------------------------------------------------------------
# STEP 8: Session info (FAIR reproducibility record)
# ---------------------------------------------------------------------------
sink(file.path(output_dir, "session_info.txt"))
cat("Script: 07.7-samova_input_prep_cp.R\n")
cat("Date:", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), "\n")
cat("Seed: 42\n")
cat(sprintf("Sites included: %d\n", n_sites))
cat(sprintf("Sites excluded: %s\n", paste(EXCLUDE_SITES, collapse = ", ")))
cat(sprintf("Samples retained: %d\n", length(aln)))
cat(sprintf("Unique haplotypes: %d\n", n_haps))
cat(sprintf("Alignment length: %d bp\n", seq_len))
cat("\n")
sessionInfo()
sink()

cat("\nDone. Outputs in:", output_dir, "\n")
