# Script to run get_gene_id and get_mlst_id and return matrix with full and partial mlst matches
source("lib/functions_lib.R")
source("lib/save_CD630_mlst_gene_sequences.R")
source('lib/create_genomic_pos_df.R')

gene_key <- read.csv("data/gene_key.csv")
genomic_pos_df <- read.csv("data/genomic_pos.csv")

# mlst profiles that list the sequence of genes ID's the define an MLST 
mlst_profiles <- read_tsv("data/mlst_profiles.txt")

#PSM001_mlst <- vcf_to_mlst("data/PSM001__aln_mpileup_raw.vcf", gene_key, genomic_pos_df, mlst_profiles)
PSM002_mlst <- vcf_to_mlst("data/PSM002__aln_mpileup_raw.vcf.gz", gene_key, genomic_pos_df, mlst_profiles)

files <- list.files(path = "/nfs/esnitkin/Project_Cdiff/Sequence_data/Project_propensity_score_match/consensus/2019_07_11_09_41_20_core_results/data_matrix/snpEff_results/", pattern = '__aln_mpileup_raw.vcf.gz$', full.names = TRUE)

predicted_MLST <- sapply(files, function(f){
  vcf_to_mlst(f, gene_key, genomic_pos_df, mlst_profiles)
})
sample_names <- gsub("/nfs/esnitkin/Project_Cdiff/Sequence_data/Project_propensity_score_match/consensus/2019_07_11_09_41_20_core_results/data_matrix/snpEff_results//|__aln_mpileup_raw.vcf.gz" , "", files)
write_csv(data.frame(sample_names, predicted_MLST), "data/predicted_MLST.csv")
