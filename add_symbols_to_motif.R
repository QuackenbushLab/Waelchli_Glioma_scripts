# Assign gene names to motifs. This only needs to be done once.
if(!require("biomaRt")){
  install.packages("biomaRt")
}
library("biomaRt")
motifFile <- NULL
motifOutFile <- NULL

motifs <- read.table(motifFile,sep = "\t")
colnames(motifs) <- c("tf", "target", "mor")
motifs <- motifs[which(motifs$mor == 1),]
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
symbolMapping <- getBM(filters= "ensembl_gene_id", attributes= c("hgnc_symbol", "ensembl_gene_id"),
                       values= unique(motifs$target),mart= mart)
symbolMappingDedup <- do.call(rbind, lapply(unique(symbolMapping$ensembl_gene_id), function(gene){
  firstMatch <- which(symbolMapping$ensembl_gene_id == gene)[1]
  uniqueMapping <- data.frame(hgnc_symbol = symbolMapping[firstMatch, "hgnc_symbol"],
                              ensembl_gene_id = gene)
  return(uniqueMapping)
}))
rownames(symbolMappingDedup) <- symbolMappingDedup$ensembl_gene_id
motifsTargetName <- symbolMappingDedup[motifs$target, "hgnc_symbol"]
motifs$target <- motifsTargetName
write.table(motifs, motifOutFile, sep = "\t", quote = FALSE)