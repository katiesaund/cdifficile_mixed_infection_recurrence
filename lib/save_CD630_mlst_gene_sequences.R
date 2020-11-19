library(tidyverse)
mlst_profiles <- read_tsv("data/mlst_profiles.txt")

# matrix with the gene, gene id, and sequence for each of the MLST genes 
gene_key <- read.csv("data/gene_key.csv")

# establish the cd630 mlst genes knowing that cd630 is mlst 54 
mlst_profiles_cd630 <- mlst_profiles %>% filter(ST == 54)
cd_630_adk <- gene_key %>% filter(`ID` == mlst_profiles_cd630$adk, gene == "adk") 
cd_630_atpA <- gene_key %>% filter(`ID` == mlst_profiles_cd630$atpA, gene == "atpA")
cd_630_dxr <- gene_key %>% filter(`ID` == mlst_profiles_cd630$dxr, gene == "dxr")
cd_630_glyA <- gene_key %>% filter(`ID` == mlst_profiles_cd630$glyA, gene == "glyA")
cd_630_recA <- gene_key %>% filter(`ID` == mlst_profiles_cd630$recA, gene == "recA")
cd_630_sodA <- gene_key %>% filter(`ID` == mlst_profiles_cd630$sodA, gene == "sodA")
cd_630_tpi <- gene_key %>% filter(`ID` == mlst_profiles_cd630$tpi, gene == "tpi")
