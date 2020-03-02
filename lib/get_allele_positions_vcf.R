# Bailey Garb 
# 19-11-15 

### Reminders ###
# The mlst definition of the gene is not the whole genome
# MLST is defines by a unique set of 7 gene ID 
# A gene ID is a specific unique squence for one of the seven MLST genes 
#
#
# Current issue. Assumed that all variantes mapped to a MLST but that doesn't seem to true. It is possible that there is a non-informant variant that doesn't define the mslt. Example: adk for PSM170 is defined as the same as CD630 but there is a variant that we found suggesting that there is a non MLST gene ID defining variant in the genome. 
# Next steps: write a function that gives  gene_ID numbers that are possible. The same gene_ID as CD630 if variant isn't found or there are no variants. It should work for all the genes and not 7 versions of the function. 

### Set up  ### 
library(vcfR)
library(tidyverse)
library(seqinr)

# load test sample which should be a population sample or a mized sample. The sample must be in vcf formatting. Vcf is a file with snp information in comparision to the referene genome CD630. 
PSM170 <- read.vcfR("../data/PSM170__aln_mpileup_raw.vcf")

# fix is just the part of the VCF file that we need which is the information about the variants and not the information about how it was sequenced. 
fix <- as.data.frame(getFIX(PSM170))
fix$POS <- as.numeric(as.character(fix$POS))

# establish the mlst variance in the sample compared to cd630 
# these have a list of all the variation between the gene for the test sample and cd630 and the position in the genome where that happens. 
# the numbers denote the positions in the genome where the NCBI defines the genes in Cdif. Example adk: https://www.ncbi.nlm.nih.gov/nuccore/NC_009089.1?report=fasta&from=113228&to=113878
adk <- fix %>% filter(POS >= 113228 & POS <= 113878)
atpA <- fix %>% filter(POS >= 3432464 & POS <= 3434242)
dxr <- fix %>% filter(POS >=  2466836 & POS <= 2467990)
glyA <- fix %>% filter(POS >= 3162737 & POS <= 3163981)
recA <- fix %>% filter(POS >= 1539910 & POS <= 1540956)
sodA <- fix %>% filter(POS >= 1889811 & POS <= 1890515)
tpi <- fix %>% filter(POS >= 3706953 & POS <= 3707696)

# matrix with the gene, gene id, and sequence for each of the MLST genes 
gene_key <- read.csv("/Users/baileygarb/Desktop/cdifficile_mixed_infection_recurrence/data/gene_key.csv")

# reference genome
cd_630 <- read.fasta("/Users/baileygarb/Desktop/cdifficile_mixed_infection_recurrence/data/cdiff_630.fasta", as.string = TRUE, seqonly = TRUE)
cd_630 <-  unlist(cd_630) %>% toString()

# mlst profiles that list the sequence of genes ID's the define an MLST 
mlst_profiles <- read_tsv("/Users/baileygarb/Desktop/greatlakes_mount/Project_Cdiff/Analysis/mixed_infection_recurrence/2019-11-08_make_gene_keys/data/mlst_profiles.txt")

# establish the cd630 mlst genes knowing that cd630 is mlst 54 
mlst_profiles_cd630 <- mlst_profiles %>% filter(ST == 54)
cd_630_adk <- gene_key %>% filter(`ID` == mlst_profiles_cd630$adk, gene == "adk") 
cd_630_atpA <- gene_key %>% filter(`ID` == mlst_profiles_cd630$atpA, gene == "atpA")
cd_630_dxr <- gene_key %>% filter(`ID` == mlst_profiles_cd630$dxr, gene == "dxr")
cd_630_glyA <- gene_key %>% filter(`ID` == mlst_profiles_cd630$glyA, gene == "glyA")
cd_630_recA <- gene_key %>% filter(`ID` == mlst_profiles_cd630$recA, gene == "recA")
cd_630_sodA <- gene_key %>% filter(`ID` == mlst_profiles_cd630$sodA, gene == "sodA")
cd_630_tpi <- gene_key %>% filter(`ID` == mlst_profiles_cd630$tpi, gene == "tpi")

#combine to make a list of all genes of cd630 that make up it's mlst 
mlst_cd630 <- rbind(cd_630_adk, cd_630_atpA, cd_630_dxr, cd_630_glyA, cd_630_recA, cd_630_sodA, cd_630_tpi) 
mlst_cd630 <- mlst_cd630[, 2:4]

# function taken from https://www.r-bloggers.com/extract-different-characters-between-two-strings-of-equal-length/ that takes two strings of equal length and returns the nucleotide and postion where they differ 
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
# Variant matrices are tables that contains the positions where cd630 varies from the gene key for adk. 
#the result is a long matrix of a lot of positions and gene_ids that define the possible variations from cd630 that would identify a gene_ID

# cutoff is a matrix of the positions in the gene key where the genes switch from one to the next
cutoff <- gene_key[match(unique(gene_key$gene), gene_key$gene),1]

## the numbers added after the variant_matrix is made is used to give the position in the whole genome and now just in reference to the gene. 

#adk
adk_variant_matrix <- NULL
for(i in cutoff[1]:(cutoff[2]-1)){
  current <- list.string.diff(toString(cd_630_adk$sequence), toString(gene_key[i,4])) %>% as.matrix()
  current <- cbind(current, as.matrix(rep(gene_key$ID[i], nrow(current))))
  adk_variant_matrix <- rbind(adk_variant_matrix, as.matrix(current)) %>% as.data.frame()
}

adk_variant_matrix$position <- as.vector(as.numeric(as.character(adk_variant_matrix$position))) + 113228
colnames(adk_variant_matrix)[4] <- "gene_ID"

filter(adk_variant_matrix, position == as.numeric(as.character(adk$POS)))
filter(adk_variant_matrix, cd630 == as.character(adk$REF), mlst == as.character(adk$ALT))

#atpA
atpA_variant_matrix <- NULL
for(i in cutoff[2]:(cutoff[3]-1)){
  current <- list.string.diff(toString(cd_630_atpA$sequence), toString(gene_key[i,4])) %>% as.matrix()
  current <- cbind(current, as.matrix(rep(gene_key$ID[i], nrow(current))))
  atpA_variant_matrix <- rbind(atpA_variant_matrix, as.matrix(current)) %>% as.data.frame()
}

atpA_variant_matrix$position <- as.vector(as.numeric(as.character(atpA_variant_matrix$position))) + 3432464
colnames(atpA_variant_matrix)[4] <- "gene_ID"

#dxr
dxr_variant_matrix <- NULL
for(i in cutoff[3]:(cutoff[4]-1)){
  current <- list.string.diff(toString(cd_630_dxr$sequence), toString(gene_key[i,4])) %>% as.matrix()
  current <- cbind(current, as.matrix(rep(gene_key$ID[i], nrow(current))))
  dxr_variant_matrix <- rbind(dxr_variant_matrix, as.matrix(current)) %>% as.data.frame()
}

dxr_variant_matrix$position <- as.vector(as.numeric(as.character(dxr_variant_matrix$position))) + 2466836
colnames(dxr_variant_matrix)[4] <- "gene_ID"

#glyA

glyA_variant_matrix <- NULL
for(i in cutoff[4]:(cutoff[5]-1)){
  current <- list.string.diff(toString(cd_630_glyA$sequence), toString(gene_key[i,4])) %>% as.matrix()
  current <- cbind(current, as.matrix(rep(gene_key$ID[i], nrow(current))))
  glyA_variant_matrix <- rbind(glyA_variant_matrix, as.matrix(current)) %>% as.data.frame()
}

glyA_variant_matrix$position <- as.vector(as.numeric(as.character(glyA_variant_matrix$position))) + 3162737
colnames(glyA_variant_matrix)[4] <- "gene_ID"

#recA

recA_variant_matrix <- NULL
for(i in cutoff[5]:(cutoff[6]-1)){
  current <- list.string.diff(toString(cd_630_recA$sequence), toString(gene_key[i,4])) %>% as.matrix()
  current <- cbind(current, as.matrix(rep(gene_key$ID[i], nrow(current))))
  recA_variant_matrix <- rbind(recA_variant_matrix, as.matrix(current)) %>% as.data.frame()
}

recA_variant_matrix$position <- as.vector(as.numeric(as.character(recA_variant_matrix$position))) + 1539910
colnames(recA_variant_matrix)[4] <- "gene_ID"


#sodA
## issue, number 340 has a different length 
sodA_variant_matrix <- NULL
for(i in cutoff[6]:339){
  current <- list.string.diff(toString(cd_630_sodA$sequence), toString(gene_key[i,4])) %>% as.matrix()
  gene_id <- as.matrix(rep(gene_key$ID[i], nrow(current)))
  current <- cbind(current, gene_id)
  sodA_variant_matrix <- rbind(sodA_variant_matrix, as.matrix(current)) %>% as.data.frame()
}

##340 issue 
# the function used to identify the differences won't work unless the sequences are the same length. 
# for sequence 340 of sodA the 299th position is missing compared to cd630. 

#a
cd_630_short_340_a <- "gatgcacttgaaccttatatagataaagaaacaatgaaactgcatcatgataagcattatcaagcttatgttgataaattaaatgctgctcttgaaaaatatcctgagctttataattattctttatgtgaattattgcaaaatttagattctttacctaaagatattgctacaactgtaagaaataatgcaggtggagcttataatcataaattcttttttgatataatgacgccagaaaaaaccataccttctgaatctttaaaagaagctattgatagagactttggttcttttg"
gene_key_short_a <- "gatgcacttgaaccttatatagataaagaaacaatgaaactgcatcatgataagcattatcaagcttatgttgataaattaaatgctgctcttgaaaaatatcctgagctttataattattctttatgtgaattattgcaaaatttagattctttacctaaagatattgctacaactgtaagaaataatgcaggtggagcttataatcataaattcttttttgatataatgacgccagaaaaaaccataccttctgaatctttaaaagaagctattgatagagactttggttcttttg"
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
  print(gene_id)
  current <- cbind(current, gene_id)
  sodA_variant_matrix <- rbind(sodA_variant_matrix, as.matrix(current)) %>% as.data.frame()
}

### 385 issue 
#sequence 385 is shorted than cd630 
cd_630_shorted_385 <- "gatgcacttgaaccttatatagataaagaaacaatgaaactgcatcatgataagcattatcaagcttatgttgataaattaaatgctgctcttgaaaaatatcctgagctttataattattctttatgtgaattattgcaaaatttagattctttacctaaagatattgctacaactgtaagaaataatgcaggtggagcttataatcataaattcttttttgatataatgacgccagaaaaaaccataccttctgaatctttaaaagaagctattgatagagactttggttcttttgaaaaatttaagcaagagttccaaaaatctgctttagatgtctttggttctggttgggcttggcttgtagctactaaagatgggaaattatctattatgactactccaaatcaggatagccctgtaagtaaaaacctaactcctataatagg"
current <- list.string.diff(toString(cd_630_shorted_385), toString(gene_key[385,4])) %>% as.matrix()
gene_id <- as.matrix(rep(gene_key$ID[385], nrow(current)))
current <- cbind(current, gene_id)
sodA_variant_matrix <- rbind(sodA_variant_matrix, as.matrix(current)) %>% as.data.frame()

### gap issue 
issue_385 <- data.frame(position = (450), cd630 = "A", mlst ="-", V4 = "66")

sodA_variant_matrix <- rbind(sodA_variant_matrix, as.matrix(issue_385)) %>% as.data.frame()


for(i in 386:(cutoff[7]-1)){
  current <- list.string.diff(toString(cd_630_sodA$sequence), toString(gene_key[i,4])) %>% as.matrix()
  current <- cbind(current, as.matrix(rep(gene_key$ID[i], nrow(current))))
  sodA_variant_matrix <- rbind(sodA_variant_matrix, as.matrix(current)) %>% as.data.frame()
}
sodA_variant_matrix$position <- as.vector(as.numeric(as.character(sodA_variant_matrix$position))) + 1889811
colnames(sodA_variant_matrix)[4] <- "gene_ID"

#fill in gaps 
sodA_gene_key <- gene_key %>% filter(gene == "sodA") 
sodA_gene_key$length = sodA_gene_key$sequence %>% str_length()
sodA_gene_key$length = as.numeric(sodA_gene_key$length)
sodA_gaps <- filter(sodA_gene_key, length == 449)

#tpi
tpi_variant_matrix <- NULL
for(i in cutoff[7]:nrow(gene_key)){
  current <- list.string.diff(toString(cd_630_tpi$sequence), toString(gene_key[i,4])) %>% as.matrix()
  current <- cbind(current, as.matrix(rep(gene_key$ID[i], nrow(current))))
  tpi_variant_matrix <- rbind(tpi_variant_matrix, as.matrix(current)) %>% as.data.frame()
}

tpi_variant_matrix$position <- as.vector(as.numeric(as.character(tpi_variant_matrix$position))) + 3706953
colnames(tpi_variant_matrix)[4] <- "gene_ID"

