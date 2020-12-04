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
get_gene_id <- function(gene_name, g_key, gene_pos_df, genomic_location, gene_spec_vcf_df, cd630_seq)
{
  variant_matrix <- NULL 
    
      # now we need to change the sequence becuase there are differences between the CD630 version 
      # and our sample version
      # Use: gene_spec_vcf_df
      
      #changing the sequence
      temp <- cd630_seq$sequence %>% toString()
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