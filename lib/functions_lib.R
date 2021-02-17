library(vcfR)
library(tidyverse)
library(seqinr)


# Functions 
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


# "atpA", gene_key, positions_df, genomic_pos_df
## DEBUG
get_gene_id <- function(gene_name, g_key, gene_pos_df, genomic_position_df, gene_spec_vcf_df, cd630_seq)
{
  variant_matrix <- NULL 
    
      # now we need to change the sequence because there are differences between the CD630 version 
      # and our sample version
      # Use: gene_spec_vcf_df
      
      #changing the sequence
      temp <- cd630_seq$sequence
      if(nrow(gene_spec_vcf_df) > 0){
        for(i in 1:nrow(gene_spec_vcf_df)){
          local_str_pos <- gene_spec_vcf_df$POS[i] - genomic_pos_df$gene_pos[genomic_pos_df$gene_id==gene_name]
          substr(temp, local_str_pos, local_str_pos) <- gene_spec_vcf_df$ALT[i]
        }
      }
      # find sequence in gene_key for gene of interest
      temp_id <- gene_key %>% filter(gene == gene_name, sequence == temp) %>% pull(ID)
      if (temp_id %>% is_empty()){
        temp_id <- "Sequence match not found"
      }
      return(temp_id)
}

#write function that takes in 6/7 mlst gene ids and returns a mlst identification (ignore sodA)
# for each mlst it should say either fulfilled, partial, not 

get_mlst_id <- function(mlst_prof, adk_id, atpA_id, dxr_id, glyA_id, recA_id, tpi_id)
{
 temp <- mlst_profiles %>% group_by(ST) %>% mutate(match_count = sum(c(adk == adk_id,atpA == atpA_id, dxr == dxr_id, glyA == glyA_id, recA == recA_id, tpi== tpi_id)) ) %>% filter(match_count > 0)
  return(temp)
  
  # what is mlst_clade?
  
}

subset_vcf <- function(vcf_path){
  
  sample_vcf <- read.vcfR(vcf_path)
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
  
  results <- list('sample_adk' = sample_adk,
                   'sample_atpA' = sample_atpA,
                   'sample_dxr' = sample_dxr,
                   'sample_glyA' = sample_glyA, 
                   'sample_recA'=sample_recA, 
                  'sample_recA' = sample_recA, 
                  'sample_sodA' = sample_sodA, 
                  'sample_tpi' = sample_tpi)
  return(results)
}

vcf_to_mlst <- function(vcf_path, g_key, gene_pos_df, cd630_seq, mlst_prof){
  temp <- subset_vcf(vcf_path)
  adk <- get_gene_id("adk", g_key, gene_pos_df, temp$sample_adk, cd_630_adk$sequence)
  print(head(adk))
  atpA <- get_gene_id("atpA", g_key, gene_pos_df, temp$sample_atpA, cd630_seq)
  dxr <- get_gene_id("dxr", g_key, gene_pos_df, temp$sample_dxr, cd630_seq)
  glyA <- get_gene_id("glyA", g_key, gene_pos_df, temp$sample_glyA, cd630_seq)
  recA <- get_gene_id("recA", g_key, gene_pos_df, temp$sample_recA, cd630_seq)
  #sodA<- get_gene_id("sodA", g_key, gene_pos_df, temp$sample_sodA, cd630_seq)
  tpi <- get_gene_id("tpi", g_key, gene_pos_df, temp$sample_tpi, cd630_seq)
  
  mlst_df <- get_mlst_id(mlst_prof, adk_id = adk, atpA_id = atpA, dxr_id = dxr,glyA_id =  glyA,  recA_id = recA,tpi_id = tpi)
  return(mlst_df)
}

