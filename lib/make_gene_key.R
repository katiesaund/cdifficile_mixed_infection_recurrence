# 2019-11-08
# Bailey Garb 

# This script was finalized on 2021-03-06
# script to make gene key for MLST
# A gene key is defined as list of sequences and their corresponding gene id that maps to an MLST identifcation. 

# Cd630 MLST is 54 

library(seqinr)
library(tidyverse)

adk <- read.fasta("data/adk.fas")
atpA <- read.fasta("data/atpA.fas")
dxr <- read.fasta("data/dxr.fas")
glyA <- read.fasta("data/glyA.fas")
recA <- read.fasta("data/recA.fas")
sodA <- read.fasta("data/sodA.fas")
tpi <- read.fasta("data/tpi.fas")

#adk 

gene_key_adk <- data.frame("ID" =  1:length(adk), gene = "adk")
seq_adk <- data.frame(sequence = rep(0, length(adk)))
for (i in 1:length(adk)){
  seq_adk[i,] <- (paste(c(adk[[i]]), collapse = ""))
  
}
gene_key_adk <- cbind(gene_key_adk, seq_adk)

#atpA

gene_key_atpA <- data.frame("ID" =  1:length(atpA), gene = "atpA")
seq_atpA <- data.frame(sequence = rep(0, length(atpA)))
for (i in 1:length(atpA)){
  seq_atpA[i,] <- (paste(c(atpA[[i]]), collapse = ""))
  
}
gene_key_atpA <- cbind(gene_key_atpA, seq_atpA)

#dxr 

gene_key_dxr <- data.frame("ID" =  1:length(dxr), gene = "dxr")
seq_dxr <- data.frame(sequence = rep(0, length(dxr)))
for (i in 1:length(dxr)){
  seq_dxr[i,] <- (paste(c(dxr[[i]]), collapse = ""))
  
}
gene_key_dxr <- cbind(gene_key_dxr, seq_dxr)

#glyA

gene_key_glyA <- data.frame("ID" =  1:length(glyA), gene = "glyA")
seq_glyA <- data.frame(sequence = rep(0, length(glyA)))
for (i in 1:length(glyA)){
  seq_glyA[i,] <- (paste(c(glyA[[i]]), collapse = ""))
  
}
gene_key_glyA <- cbind(gene_key_glyA, seq_glyA)


#recA
gene_key_recA <- data.frame("ID" =  1:length(recA), gene = "recA")
seq_recA <- data.frame(sequence = rep(0, length(recA)))
for (i in 1:length(recA)){
  seq_recA[i,] <- (paste(c(recA[[i]]), collapse = ""))
  
}
gene_key_recA <- cbind(gene_key_recA, seq_recA)


#sodA
gene_key_sodA <- data.frame("ID" =  1:length(sodA), gene = "sodA")
seq_sodA <- data.frame(sequence = rep(0, length(sodA)))
for (i in 1:length(sodA)){
  seq_sodA[i,] <- (paste(c(sodA[[i]]), collapse = ""))
  
}
gene_key_sodA <- cbind(gene_key_sodA, seq_sodA)

#tpi
gene_key_tpi <- data.frame("ID" =  1:length(tpi), gene = "tpi")
seq_tpi <- data.frame(sequence = rep(0, length(tpi)))
for (i in 1:length(tpi)){
  seq_tpi[i,] <- (paste(c(tpi[[i]]), collapse = ""))
  
}
gene_key_tpi <- cbind(gene_key_tpi, seq_tpi)

#combine all gene keys

gene_key <- rbind(gene_key_adk, gene_key_atpA, gene_key_dxr, gene_key_glyA, gene_key_recA, gene_key_sodA, gene_key_tpi)

write.csv(gene_key, file = "data/gene_key.csv")
