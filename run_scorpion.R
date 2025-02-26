# Load SCORPION.
if(!require("SCORPION")){
  install.packages("SCORPION")
}
library("SCORPION")
if(!require("Seurat")){
  install.packages("Seurat")
}
library("Seurat")
if(!require(reshape2)){
  install.packages(reshape2)
}
library(reshape2)

# Get parameters.
args <- commandArgs(trailingOnly = TRUE)
print(paste0(args[2], ".csv"))
motifFile <- NULL
ppiFile <- NULL

# Load PANDA motif and PPI priors.
motifs <- read.table(motifFile,sep = "\t")
ppi <- read.table(ppiFile,sep = "\t")
colnames(ppi) <- c("protein1", "protein2", "combined_score")

# Load file.
directory <- args[1]
expression <- Seurat::Read10X(paste0(directory, "/"))

# Prepare for SCORPION.
scorpionInput <- list(tf = motifs, ppi = ppi, gex = expression)

# Run SCORPION.
scorpionOutput <- scorpion(tfMotifs = scorpionInput$tf,
                           gexMatrix = scorpionInput$gex,
                           ppiNet = scorpionInput$ppi)
meltedScorpionOutput <-melt(as.matrix(scorpionOutput$regNet))
colnames(meltedScorpionOutput) <- c("tf", "gene", "score")
write.csv(meltedScorpionOutput, paste0(args[2], ".csv"))
