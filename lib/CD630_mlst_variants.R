# Bailey Garb 
# 2020-11-17

### This is a script to identify variants for each gene key vs. CD630. Store data in matrix with 
### same format as new_genes_matrices. This will be used to compare new variant fulfillment object
### vs. new_genes_matrices. 


### Set up  ### 
library(vcfR)
library(tidyverse)
library(seqinr)

# Load in cd_630 fasta file. 
cd_630 <- read.fasta("data/cdiff_630.fasta", as.string = TRUE, seqonly = TRUE)

# mlst profiles that list the sequence of genes ID's the define an MLST 
# mlst_profiles <- read_tsv("data/mlst_profiles.txt")
# 
# # matrix with the gene, gene id, and sequence for each of the MLST genes 
# gene_key <- read.csv("data/gene_key.csv")
# 
# # establish the cd630 mlst genes knowing that cd630 is mlst 54 
# mlst_profiles_cd630 <- mlst_profiles %>% filter(ST == 54)
# cd_630_adk <- gene_key %>% filter(`ID` == mlst_profiles_cd630$adk, gene == "adk") 
# cd_630_atpA <- gene_key %>% filter(`ID` == mlst_profiles_cd630$atpA, gene == "atpA")
# cd_630_dxr <- gene_key %>% filter(`ID` == mlst_profiles_cd630$dxr, gene == "dxr")
# cd_630_glyA <- gene_key %>% filter(`ID` == mlst_profiles_cd630$glyA, gene == "glyA")
# cd_630_recA <- gene_key %>% filter(`ID` == mlst_profiles_cd630$recA, gene == "recA")
# cd_630_sodA <- gene_key %>% filter(`ID` == mlst_profiles_cd630$sodA, gene == "sodA")
# cd_630_tpi <- gene_key %>% filter(`ID` == mlst_profiles_cd630$tpi, gene == "tpi")

# Define starting positions for MLST genes 
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

#combine to make a list of all genes of cd630 that make up it's mlst 
mlst_cd630 <- rbind(cd_630_adk, cd_630_atpA, cd_630_dxr, cd_630_glyA, cd_630_recA, cd_630_sodA, cd_630_tpi) 
mlst_cd630 <- mlst_cd630[, 2:4]

# function taken from https://www.r-bloggers.com/extract-different-characters-between-two-strings-of-equal-length/ 
#that takes two strings of equal length and returns the nucleotide and position where they differ 
list.string.diff <- function(a, b, exclude = c("-", "?"), ignore.case = TRUE, show.excluded = FALSE)
{
  if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
  if(ignore.case)
  {
    a <- toupper(a)
    b <- toupper(b)
  }
  split_seqs <- strsplit(c(a, b), split = "")
  only.diff <- (split_seqs[[1]] != split_seqs[[2]])
  only.diff[
    (split_seqs[[1]] %in% exclude) |
      (split_seqs[[2]] %in% exclude)
  ] <- NA
  diff.info<-data.frame(which(is.na(only.diff)|only.diff),
                        split_seqs[[1]][only.diff],split_seqs[[2]][only.diff])
  names(diff.info)<-c("position","cd630","mlst")
  if(!show.excluded) diff.info<-na.omit(diff.info)
  diff.info
}

### make variant matrices ###
# Variant matrices are tables that contain the positions where cd630 varies from the gene key for adk. 
#the result is a long matrix of a lot of positions and gene_ids that define the possible variations from cd630 that would identify a gene_ID

# cutoff is a matrix of the positions in the gene key where the genes switch from one to the next
# Not used in functionized version but still need for sodA
cutoff <- gene_key[match(unique(gene_key$gene), gene_key$gene),1] %>% t() 
cutoff[length(cutoff)+1] <- nrow(gene_key)
cutoff <- cutoff %>% as.matrix()
rownames(cutoff) <- c("adk", 'atpA', "dxr", 'glyA', 'recA','sodA','tpi', "end")
cutoff <- cutoff %>% as.data.frame() %>% t()
## the numbers added after the variant_matrix is made is used to give the position in the whole genome and now just in reference to the gene. 
genes <- c("adk", 'atpA', "dxr", 'glyA', 'recA','sodA','tpi')
gene_pos <-list(tibble(1:57), tibble(58:118), tibble(119:179), tibble(180:269), tibble(270:319), tibble(320:390), tibble(391:473))
positions_df <- tibble(genes, gene_pos)

make_variant_matrix <- function(gene_name, g_key, gene_pos_df, genomic_location)
{
  variant_matrix <- NULL 
  # for i in the length of each mlst gene
  for(i in gene_pos_df$gene_pos[gene_pos_df$genes==gene_name] %>% unlist() %>% unname()){
    #get sequence for cd630 gene 
    temp <- eval(as.name(paste("cd_630_",gene_name, sep = "")))$sequence
    # get string difference between cd630 gene and test sequence 
    current <- list.string.diff(toString(temp), toString(g_key[i, 4])) %>% as.matrix()
    # add gene key information to the difference in current
    current <- cbind(current, as.matrix(rep(g_key$ID[i], nrow(current))))
    # add current to variant matrix
    variant_matrix <- rbind(variant_matrix, as.matrix(current)) %>% as.data.frame()
  }
  #print(as.vector(as.numeric(as.character(variant_matrix$position))) + genomic_location %>% filter(genes == gene_name) %>% pull(gene_pos))
  variant_matrix$position <- as.vector(as.numeric(as.character(variant_matrix$position))) + genomic_location %>% filter(genes == gene_name) %>% pull(gene_pos)
  colnames(variant_matrix)[4] <- "gene_ID"
  return(variant_matrix)
}

cd630_adk_matrix <- make_variant_matrix("adk", gene_key, positions_df, genomic_pos_df)
cd630_atpA_matrix <- make_variant_matrix("atpA", gene_key, positions_df, genomic_pos_df)
cd630_dxr_matrix <- make_variant_matrix("dxr", gene_key, positions_df, genomic_pos_df)
cd630_glyA_matrix <- make_variant_matrix("glyA", gene_key, positions_df, genomic_pos_df)
cd630_recA_matrix <- make_variant_matrix("recA", gene_key, positions_df, genomic_pos_df)
cd630_tpi_matrix <- make_variant_matrix("tpi", gene_key, positions_df, genomic_pos_df)

#sodA
## issue, number 340 has a different length 
cd630_sodA_matrix <- NULL
for(i in cutoff[[6]]:339){
  current <- list.string.diff(toString(cd_630_sodA$sequence), toString(gene_key[i,4])) %>% as.matrix()
  gene_id <- as.matrix(rep(gene_key$ID[i], nrow(current)))
  current <- cbind(current, gene_id)
  cd630_sodA_matrix <- rbind(cd630_sodA_matrix, as.matrix(current)) %>% as.data.frame()
}

##340 issue 
# the function used to identify the differences won't work unless the sequences are the same length. 
# for sequence 340 of sodA the 299th position is missing compared to cd630. 

#a
cd_630_short_340_a <- "gatgcacttgaaccttatatagataaagaaacaatgaaactgcatcatgataagcattatcaagcttatgttgataaattaaatgctgctcttgaaaaatatcctgagctttat
aattattctttatgtgaattattgcaaaatttagattctttacctaaagatattgctacaactgtaagaaataatgcaggtggagctt
ataatcataaattcttttttgatataatgacgccagaaaaaaccataccttctgaatctttaaaagaagctattgatagagactttggttcttttg"
gene_key_short_a <- "gatgcacttgaaccttatatagataaagaaacaatgaaactgcatcatgataagcattatcaagcttatgttgataaat
taaatgctgctcttgaaaaatatcctgagctttataattattctttatgtgaattattgcaaaatttagattctttacctaaagatattgctacaactgtaagaaataatgcaggtggagcttataatcataaattct
tttttgatataatgacgccagaaaaaaccataccttctgaatctttaaaagaagctattgatagagactttggttcttttg"
current <- list.string.diff(toString(cd_630_short_340_a), toString(gene_key_short_a)) %>% as.matrix()
gene_id <- as.matrix(rep(gene_key$ID[340], nrow(current)))
current <- cbind(current, gene_id)
cd630_sodA_matrix <- rbind(cd630_sodA_matrix, as.matrix(current)) %>% as.data.frame()

#gap
issue_340 <- data.frame(position = (299), cd630 = "A", mlst ="-", V4 = "21")
cd630_sodA_matrix <- rbind(cd630_sodA_matrix, as.matrix(issue_340)) %>% as.data.frame()

#b
cd_630_short_340_b <- "aaaatttaagcaagagttccaaaaatctgctttagatgtctttggttctggttgggcttggcttgtagctactaaagatgggaaattatctattatgactactccaaatcaggatagccctgtaagtaaaaacctaactcctataatagga"
gene_key_short_b <- "aaaatttaagcaagagttccaaaaatctgctttagatgtctttggttctggttgggcttggcttgtagctactaaagatggcaaattatctattatgactactccaaatcaggatagccctgtaagtaaaaacctaactcctataatagga"
current <- list.string.diff(toString(cd_630_short_340_b), toString(gene_key_short_b)) %>% as.matrix()
gene_id <- as.matrix(rep(gene_key$ID[340], nrow(current)))
current <- cbind(current, gene_id)
cd630_sodA_matrix <- rbind(cd630_sodA_matrix, as.matrix(current)) %>% as.data.frame()


for(i in 341:384){
  current <- list.string.diff(toString(cd_630_sodA$sequence), toString(gene_key[i,4])) %>% as.matrix()
  gene_id <- as.matrix(rep(gene_key$ID[i], nrow(current)))
  current <- cbind(current, gene_id)
  cd630_sodA_matrix <- rbind(cd630_sodA_matrix, as.matrix(current)) %>% as.data.frame()
}

### 385 issue 
#sequence 385 is shorter than cd630 
cd_630_shorted_385 <- cd_630_sodA$sequence %>% substring(1, str_length(cd_630_sodA$sequence)-1)
current <- list.string.diff(toString(cd_630_shorted_385), toString(gene_key[385,4])) %>% as.matrix()
gene_id <- as.matrix(rep(gene_key$ID[385], nrow(current)))
current <- cbind(current, gene_id)
cd630_sodA_matrix <- rbind(cd630_sodA_matrix, as.matrix(current)) %>% as.data.frame()

### gap issue 
issue_385 <- data.frame(position = (450), cd630 = "A", mlst ="-", V4 = "66")

cd630_sodA_matrix <- rbind(cd630_sodA_matrix, as.matrix(issue_385)) %>% as.data.frame()


for(i in 386:(cutoff[[7]]-1)){
  current <- list.string.diff(toString(cd_630_sodA$sequence), toString(gene_key[i,4])) %>% as.matrix()
  current <- cbind(current, as.matrix(rep(gene_key$ID[i], nrow(current))))
  cd630_sodA_matrix <- rbind(cd630_sodA_matrix, as.matrix(current)) %>% as.data.frame()
}
cd630_sodA_matrix$position <- as.vector(as.numeric(as.character(cd630_sodA_matrix$position))) + sodA_pos
colnames(cd630_sodA_matrix)[4] <- "gene_ID"

#fill in gaps 
sodA_gene_key <- gene_key %>% filter(gene == "sodA") 
sodA_gene_key$length = sodA_gene_key$sequence %>% str_length()
sodA_gene_key$length = as.numeric(sodA_gene_key$length)
sodA_gaps <- filter(sodA_gene_key, length == 449)