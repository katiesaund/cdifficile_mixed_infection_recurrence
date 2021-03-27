# Define starting positions for MLST genes 

#This script takes in the positions within the Cdiff genome where every MLST gene starts and outputs a csv files 
#that has gene_id and gene_position information for each gene. 

library(tidyverse)
# Cd 630: https://www.ncbi.nlm.nih.gov/nuccore/AM180355.1
source('~/thesis/cdifficile_mixed_infection_recurrence/lib/gene_positions.R')


#Make positions df with make of gene corresponding to the starting position of the gene in cd630 genome.
#gene_id is the gene name and gene_pos is the starting position. 
genomic_pos_df <- as.data.frame(matrix(0, nrow = 7, ncol = 2))
genomic_pos_df[, 1] <- c("adk", "atpA", "dxr", "glyA", "recA", "sodA", "tpi")
genomic_pos_df[, 2] <- c(adk_start, 
                         atpA_start, 
                         dxr_start, 
                         glyA_start, 
                         recA_start, 
                         sodA_start, 
                         tpi_start) 
colnames(genomic_pos_df) <- c("gene_id", "gene_pos")
write_csv(genomic_pos_df, "data/genomic_pos.csv")
