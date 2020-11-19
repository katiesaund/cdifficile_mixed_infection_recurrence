# Define starting positions for MLST genes 

# TODO add comment to describe script goals/inputs/outputs

library(tidyverse)

adk_pos <- 113228
atpA_pos <- 3432464
dxr_pos <- 2466836
glyA_pos <- 3162737
recA_pos <- 1539910
sodA_pos <- 1889811
tpi_pos <- 3706953

#Make positions df with make of gene corresponding to the starting position of the gene in cd630 genome.
#gene_id is the gene name and gene_pos is the starting position. 
genomic_pos_df <- as.data.frame(matrix(0, nrow = 7, ncol = 2))
genomic_pos_df[, 1] <- c("adk", "atpA", "dxr", "glyA", "recA", "sodA", "tpi")
genomic_pos_df[, 2] <- c(adk_pos, 
                         atpA_pos, 
                         dxr_pos, 
                         glyA_pos, 
                         recA_pos, 
                         sodA_pos, 
                         tpi_pos) 
colnames(genomic_pos_df) <- c("gene_id", "gene_pos")
write_csv(genomic_pos_df, "data/genomic_pos.csv")
