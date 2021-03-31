# Script to run get_gene_id and get_mlst_id and return matrix with full and partial mlst matches

source("lib/functions_lib.R")
source("lib/save_CD630_mlst_gene_sequences.R")
source('~/thesis/cdifficile_mixed_infection_recurrence/lib/create_genomic_pos_df.R')

gene_key <- read.csv("data/gene_key.csv")
genomic_pos_df <- read.csv("data/genomic_pos.csv")

# reference genome
cd_630 <- read.fasta("data/cdiff_630.fasta", as.string = TRUE, seqonly = TRUE)
cd_630 <- unlist(cd_630) %>% toString()

# mlst profiles that list the sequence of genes ID's the define an MLST 
mlst_profiles <- read_tsv("data/mlst_profiles.txt")

PSM001_mlst <- vcf_to_mlst("data/PSM001__aln_mpileup_raw.vcf", gene_key, genomic_pos_df, mlst_profiles)

PSM170_mlst <- vcf_to_mlst("data/PSM170__aln_mpileup_raw.vcf", gene_key, genomic_pos_df, mlst_profiles)
PSM007_mlst <- vcf_to_mlst("data/PSM007__aln_mpileup_raw.vcf", gene_key, genomic_pos_df, mlst_profiles)

PSM008_mlst <- vcf_to_mlst("data/PSM008__aln_mpileup_raw.vcf", gene_key, genomic_pos_df, mlst_profiles)

PSM016_mlst <- vcf_to_mlst("data/PSM016__aln_mpileup_raw.vcf", gene_key, genomic_pos_df, mlst_profiles)

# NOTES: 
#psm170 and show that they are different.
# write a read_me to explain how this works at the moment and outline how to use this for the bigger picture 