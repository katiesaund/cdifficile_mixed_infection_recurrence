# Bailey Garb 
# 19-11-15 
### Reminders ###
# PSM
# The mlst definition of the gene is not the whole genome
# MLST is defines by a unique set of 7 gene ID 
# A gene ID is a specific unique squence for one of the seven MLST genes 
#
#
# Current issue. Assumed that all variantes mapped to a MLST but that doesn't seem to true. It is possible that there is a non-informant variant that doesn't define the mslt. Example: adk for PSM170 is defined as the same as 
# CD630 but there is a variant that we found suggesting that there is a non MLST gene ID defining variant in the genome. 
# Next steps: write a function that gives  gene_ID numbers that are possible. The same gene_ID as CD630 if variant isn't found or there are no variants. It should work for all the genes and not 7 versions of the function. 

### Set up  ### 
source("lib/functions_lib.R")
source("lib/save_CD630_mlst_gene_sequences.R")

# load test sample which should be a population sample or a mixed sample. The sample must be in vcf formatting. Vcf is a file with snp information in comparision to the referene genome CD630. 

# PSM170 is now called sample_vcf
sample_vcf <- read.vcfR("data/PSM170__aln_mpileup_raw.vcf")
# PSM001 <- read.vcfR("data/PSM001__aln_mpileup_raw.vcf")
# fix is just the part of the VCF file that we need which is the information about the variants and not the information about how it was sequenced. 
# Fix is now called sample_variant_df
sample_variant_df <- as.data.frame(getFIX(sample_vcf))
sample_variant_df$POS <- as.numeric(as.character(sample_variant_df$POS))

# establish the mlst variance in the sample compared to cd630 
# these have a list of all the variation between the gene for the test sample and cd630 and the position in the genome where that happens. 
# the numbers denote the positions in the genome where the NCBI defines the genes in Cdif. Example adk: https://www.ncbi.nlm.nih.gov/nuccore/NC_009089.1?report=fasta&from=113228&to=113878
sample_adk <- sample_variant_df %>% filter(POS >= 113228 & POS <= 113878)
sample_atpA <- sample_variant_df %>% filter(POS >= 3432464 & POS <= 3434242)
sample_dxr <- sample_variant_df %>% filter(POS >=  2466836 & POS <= 2467990)
sample_glyA <- sample_variant_df %>% filter(POS >= 3162737 & POS <= 3163981)
sample_recA <- sample_variant_df %>% filter(POS >= 1539910 & POS <= 1540956)
sample_sodA <- sample_variant_df %>% filter(POS >= 1889811 & POS <= 1890515)
sample_tpi <- sample_variant_df %>% filter(POS >= 3706953 & POS <= 3707696)


# matrix with the gene, gene id, and sequence for each of the MLST genes 
gene_key <- read.csv("data/gene_key.csv")
genomic_pos_df <- read.csv("data/genomic_pos.csv")

# reference genome
cd_630 <- read.fasta("data/cdiff_630.fasta", as.string = TRUE, seqonly = TRUE)
cd_630 <- unlist(cd_630) %>% toString()

# mlst profiles that list the sequence of genes ID's the define an MLST 
mlst_profiles <- read_tsv("data/mlst_profiles.txt")

### make variant matrices ###
# Variant matrices are tables that contain the positions where cd630 varies from the gene key for adk. 
#the result is a long matrix of a lot of positions and gene_ids that define the possible variations from cd630 that would identify a gene_ID

# cutoff is a matrix of the positions in the gene key where the genes switch from one to the next
cutoff <- gene_key[match(unique(gene_key$gene), gene_key$gene), 1] %>% t() 
cutoff[length(cutoff) + 1] <- nrow(gene_key)
cutoff <- cutoff %>% as.matrix()
rownames(cutoff) <- c("adk", 'atpA', "dxr", 'glyA', 'recA','sodA','tpi', "end")
cutoff <- cutoff %>% as.data.frame() %>% t()
## the numbers added after the variant_matrix is made is used to give the position in the whole genome and now just in reference to the gene. 
genes <- c("adk", 'atpA', "dxr", 'glyA', 'recA','sodA','tpi')
gene_pos <- list(tibble(1:57), tibble(58:118), tibble(119:179), tibble(180:269), tibble(270:319), tibble(320:390), tibble(391:473))
positions_df <- tibble(genes, gene_pos)

# Generate Variant Matrices for each gene
new_atpA <- get_gene_id("atpA", gene_key, positions_df, genomic_pos_df, sample_adk, cd_630_atpA)
new_adk <- get_gene_id("adk", gene_key, positions_df, genomic_pos_df, sample_adk, cd_630_adk)
# don't forget to add sample and cd630 version of each gene as input variables below: 
new_dxr <- get_gene_id("dxr", gene_key, positions_df, genomic_pos_df)
new_glyA <- get_gene_id("glyA", gene_key, positions_df, genomic_pos_df)
new_recA <- get_gene_id("recA", gene_key, positions_df, genomic_pos_df)
new_tpi <- get_gene_id("tpi", gene_key, positions_df, genomic_pos_df)

#sodA
## issue, number 340 has a different length 
sodA_variant_matrix <- NULL
for(i in cutoff[[6]]:339){
  current <- list.string.diff(toString(cd_630_sodA$sequence), toString(gene_key[i,4])) %>% as.matrix()
  gene_id <- as.matrix(rep(gene_key$ID[i], nrow(current)))
  current <- cbind(current, gene_id)
  sodA_variant_matrix <- rbind(sodA_variant_matrix, as.matrix(current)) %>% as.data.frame()
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
sodA_variant_matrix <- rbind(sodA_variant_matrix, as.matrix(current)) %>% as.data.frame()

#gap
issue_340 <- data.frame(position = (299), cd630 = "A", mlst ="-", V4 = "21")
sodA_variant_matrix <- rbind(sodA_variant_matrix, as.matrix(issue_340)) %>% as.data.frame()

#b
cd_630_short_340_b <- "aaaatttaagcaagagttccaaaaatctgctttagatgtctttggttctggttgggcttggcttgtagctactaaagatgggaaattatctattatgactactccaaatcaggatagccctgtaagtaaaaacctaactcctataatagga"
gene_key_short_b <- "aaaatttaagcaagagttccaaaaatctgctttagatgtctttggttctggttgggcttggcttgtagctactaaagatggcaaattatctattatgactactccaaatcaggatagccctgtaagtaaaaacctaactcctataatagga"
current <- list.string.diff(toString(cd_630_short_340_b), toString(gene_key_short_b)) %>% as.matrix()
gene_id <- as.matrix(rep(gene_key$ID[340], nrow(current)))
current <- cbind(current, gene_id)
sodA_variant_matrix <- rbind(sodA_variant_matrix, as.matrix(current)) %>% as.data.frame()


for(i in 341:384){
  current <- list.string.diff(toString(cd_630_sodA$sequence), toString(gene_key[i,4])) %>% as.matrix()
  gene_id <- as.matrix(rep(gene_key$ID[i], nrow(current)))
  current <- cbind(current, gene_id)
  sodA_variant_matrix <- rbind(sodA_variant_matrix, as.matrix(current)) %>% as.data.frame()
}

### 385 issue 
#sequence 385 is shorter than cd630 
cd_630_shorted_385 <- cd_630_sodA$sequence %>% substring(1, str_length(cd_630_sodA$sequence)-1)
current <- list.string.diff(toString(cd_630_shorted_385), toString(gene_key[385,4])) %>% as.matrix()
gene_id <- as.matrix(rep(gene_key$ID[385], nrow(current)))
current <- cbind(current, gene_id)
sodA_variant_matrix <- rbind(sodA_variant_matrix, as.matrix(current)) %>% as.data.frame()

### gap issue 
issue_385 <- data.frame(position = (450), cd630 = "A", mlst ="-", V4 = "66")

sodA_variant_matrix <- rbind(sodA_variant_matrix, as.matrix(issue_385)) %>% as.data.frame()


for(i in 386:(cutoff[[7]]-1)){
  current <- list.string.diff(toString(cd_630_sodA$sequence), toString(gene_key[i,4])) %>% as.matrix()
  current <- cbind(current, as.matrix(rep(gene_key$ID[i], nrow(current))))
  sodA_variant_matrix <- rbind(sodA_variant_matrix, as.matrix(current)) %>% as.data.frame()
}
sodA_variant_matrix$position <- as.vector(as.numeric(as.character(sodA_variant_matrix$position))) + sodA_pos
colnames(sodA_variant_matrix)[4] <- "gene_ID"

#fill in gaps 
sodA_gene_key <- gene_key %>% filter(gene == "sodA") 
sodA_gene_key$length = sodA_gene_key$sequence %>% str_length()
sodA_gene_key$length = as.numeric(sodA_gene_key$length)
sodA_gaps <- filter(sodA_gene_key, length == 449)



