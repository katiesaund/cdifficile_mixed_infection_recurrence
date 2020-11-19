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
make_variant_matrix <- function(gene_name, g_key, gene_pos_df, genomic_location, gene_spec_vcf_df, cd630_seq)
{
  variant_matrix <- NULL 
  # for i in the length of each mlst gene
  for(i in gene_pos_df$gene_pos[gene_pos_df$genes==gene_name] %>% unlist() %>% unname()){
    #get sequence for cd630 gene 
    #### THIS IS THE ISSUE ###   temp <- eval(as.name(paste("cd_630_",gene_name, sep = "")))$sequence
    temp <- cd630_seq$sequence
    
    if (nrow(gene_spec_vcf_df) > 0) {
      
      # now we need to change the sequence becuase there are differences between the CD630 version and our sample version
      # Use: gene_spec_vcf_df
      
      # get string difference between cd630 gene and test sequence 
      current <- list.string.diff(toString(temp), toString(g_key[i, 4])) %>% as.matrix()
      # add gene key information to the difference in current
      current <- cbind(current, as.matrix(rep(g_key$ID[i], nrow(current))))
      # add current to variant matrix
      variant_matrix <- rbind(variant_matrix, as.matrix(current)) %>% as.data.frame()
    } else {
      # same as CD630
    }
  }
  #print(as.vector(as.numeric(as.character(variant_matrix$position))) + genomic_location %>% filter(genes == gene_name) %>% pull(gene_pos))
  variant_matrix$position <- as.vector(as.numeric(as.character(variant_matrix$position))) + genomic_location %>% filter(genes == gene_name) %>% pull(gene_pos)
  colnames(variant_matrix)[4] <- "gene_ID"
  return(variant_matrix)
}
