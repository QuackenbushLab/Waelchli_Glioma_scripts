setwd("/home/ubuntu/Waelchli_code/netZooR/netZooR")
library("roxygen2")
roxygen2::roxygenize()

# Set up directories.
sourceDirectory <- NULL
resultDirectory <- NULL
allGenes <- c(c("SOX4", "CD93", "COL4A2", "CALCRL", "HSP90B1", "COL4A1", "AKAP13", "HSPG"),
              c("ICAM1", "ICAM2", "ICAM3", "VCAM1", "SELP", "SELE", "MADCAM1",  "CD34", "GLYCAM1"),
              c("LGALS1", "LGALS9", "CD274", "FASLG", "HAVCR2", "PDCD1LG2", "IDO1", "STAB1"))
nullFile <- NULL
dir.create(resultDirectory)

# Read null distribution.
pandas<-sample(readRDS(nullFile), 1000)

# Separate out the recurrent and de novo GBMs.
# Function to read networks for a group and only return positive or thresholded values.
ReadNetworks<- function(group){
  
  # Read the files.
  files <- lapply(group, function(file){
    network <- read.csv(paste0(sourceDirectory, file), row.names = 1)
    rownames(network) <- paste(network$tf, network$gene, sep = "__")
    cat(".")
    return(network)
  })
  
  # Find the union of all edges.
  allRownames <- unique(unlist(lapply(files, function(net){
    return(rownames(net))
  })))
  
  # Create a new data frame.
  commonNets <- data.frame(tf = unlist(lapply(allRownames, function(row){return(strsplit(row, "__")[[1]][1])})),
                           gene = unlist(lapply(allRownames, function(row){return(strsplit(row, "__")[[1]][2])})))
  rownames(commonNets) <- allRownames
  for(i in 1:length(files)){
    commonNets[,i+2] <- NA
  }
  
  # Add the scores.
  for(i in 1:length(files)){
    commonNets[rownames(files[[i]]),i+2] <- files[[i]]$score
    str(commonNets)
  }
  
  return(commonNets)
}
denovoGBM <- ReadNetworks(c("GSM8101605_GBM_patient6_EC.csv",
                           "GSM8101549_GBM_28_08_2018.csv",
                           "GSM8101550_GBM_07122018.csv",
                           "GSM8101551_GBM_15-10-2019_ECs.csv",
                           "GSM8101552_GBM_EC_19_11_19.csv",
                           "GSM8101553_GBM_03-03-2020_sorted_3pr_v3.csv"))
saveRDS(denovoGBM, paste0(sourceDirectory, "denovoGBM.RDS"))
recurrentGBM <- ReadNetworks(c("GSM8101604_GBM_patient7_EC.csv",
                            "GSM8101608_GBM_patient8_EC.csv"))
saveRDS(recurrentGBM, paste0(sourceDirectory, "recurrentGBM.RDS"))

# Run BLOBFISH with precomputed p-values.
# Fetal CNS
fetalCNS <- readRDS(paste0(sourceDirectory, "fetalCNSNets.RDS"))
str(setdiff(allGenes, unique(fetalCNS$gene)))
fetalCNSSubnet <- RunBLOBFISH(geneSet = intersect(allGenes, fetalCNS$gene), 
                               pValueFile = paste0(resultDirectory, "fetalCNSPval.RDS"), loadPValues = TRUE,
                               networks = fetalCNS, alpha = 0.005, hopConstraint = 2, nullDistribution = pandas,
                               verbose = TRUE)
write.csv(fetalCNSSubnet, paste0(resultDirectory, "fetalCNS.csv"))

# Fetal periphery
fetalPeriphery <- readRDS(paste0(sourceDirectory, "fetalPeriphery.RDS"))
str(setdiff(allGenes, unique(fetalPeriphery$gene)))
fetalPeripherySubnet <- RunBLOBFISH(geneSet = intersect(allGenes, fetalPeriphery$gene), 
                                    pValueFile = paste0(resultDirectory, "fetalPeripheryPval.RDS"),
                                    loadPValues = FALSE,
                              networks = fetalPeriphery, alpha = 0.005, hopConstraint = 2, nullDistribution = pandas,
                              verbose = TRUE)
write.csv(fetalPeripherySubnet, paste0(resultDirectory, "fetalPeriphery.csv"))

# Adult controls
adultControls <- readRDS(paste0(sourceDirectory, "adultControls.RDS"))
str(setdiff(allGenes, unique(adultControls$gene)))
adultControlSubnet <- RunBLOBFISH(geneSet = intersect(allGenes, adultControls$gene), 
                                  pValueFile = paste0(resultDirectory, "adultControlsPval.RDS"),loadPValues = FALSE,
                                    networks = adultControls, alpha = 0.005, hopConstraint = 2, 
                                  nullDistribution = pandas,
                                    verbose = TRUE)
write.csv(adultControlSubnet, paste0(resultDirectory, "adultControl.csv"))

# AVM
avm <- readRDS(paste0(sourceDirectory, "AVM.RDS"))
str(setdiff(allGenes, unique(avm$gene)))
avmSubnet <- RunBLOBFISH(geneSet = intersect(allGenes, avm$gene), 
                         pValueFile = paste0(resultDirectory, "avmPval.RDS"),loadPValues = FALSE,
                                  networks = avm, alpha = 0.005, hopConstraint = 2, nullDistribution = pandas,
                                  verbose = TRUE)
write.csv(avmSubnet, paste0(resultDirectory, "avm.csv"))

# GBM
gbm <- readRDS(paste0(sourceDirectory, "GBM.RDS"))
str(setdiff(allGenes, unique(gbm$gene)))
gbmSubnet <- RunBLOBFISH(geneSet = intersect(allGenes, gbm$gene), 
                         pValueFile = paste0(resultDirectory, "gbmPval.RDS"),loadPValues = FALSE,
                         networks = gbm, alpha = 0.005, hopConstraint = 2, nullDistribution = pandas,
                         verbose = TRUE)
write.csv(gbmSubnet, paste0(resultDirectory, "gbm.csv"))

str(setdiff(allGenes, unique(lgg$gene)))

# LGG
lgg <- readRDS(paste0(sourceDirectory, "LGG.RDS"))
lggSubnet <- RunBLOBFISH(geneSet = intersect(allGenes, lgg$gene), 
                         pValueFile = paste0(resultDirectory, "lggPval.RDS"),loadPValues = FALSE,
                         networks = lgg, alpha = 0.005, hopConstraint = 2, nullDistribution = pandas,
                         verbose = TRUE)
write.csv(lggSubnet, paste0(resultDirectory, "lgg.csv"))

# MET
lungMetastasis <- readRDS(paste0(sourceDirectory, "MET.RDS"))
str(setdiff(allGenes, unique(lungMetastasis$gene)))
lungMetastasisSubnet <- RunBLOBFISH(geneSet = intersect(allGenes, lungMetastasis$gene), 
                                    loadPValues = TRUE,pValueFile = paste0(resultDirectory, "METPval.RDS"),
                         networks = lungMetastasis, alpha = 0.005, hopConstraint = 2, nullDistribution = pandas,
                         verbose = TRUE)
write.csv(lungMetastasisSubnet, paste0(resultDirectory, "lungMetastasis.csv"))

# MEN
meningioma <- readRDS(paste0(sourceDirectory, "MEN.RDS"))
str(setdiff(allGenes, unique(meningioma$gene)))
meningiomaSubnet <- RunBLOBFISH(geneSet = intersect(allGenes, meningioma$gene), 
                                pValueFile = paste0(resultDirectory, "MENPval.RDS"),loadPValues = FALSE,
                                    networks = meningioma, alpha = 0.005, hopConstraint = 2, nullDistribution = pandas,
                                    verbose = TRUE)
write.csv(meningiomaSubnet, paste0(resultDirectory, "meningioma.csv"))

# De Novo GBM
meningioma <- readRDS(paste0(sourceDirectory, "denovoGBM.RDS"))
str(setdiff(allGenes, unique(meningioma$gene)))
meningiomaSubnet <- RunBLOBFISH(geneSet = intersect(allGenes, meningioma$gene), 
                                pValueFile = paste0(resultDirectory, "denovoGBMPval.RDS"),loadPValues = FALSE,
                                networks = meningioma, alpha = 0.005, hopConstraint = 2, nullDistribution = pandas,
                                verbose = TRUE)
write.csv(meningiomaSubnet, paste0(resultDirectory, "denovoGBM.csv"))

# Recurrent GBM
meningioma <- readRDS(paste0(sourceDirectory, "recurrentGBM.RDS"))
str(setdiff(allGenes, unique(meningioma$gene)))
meningiomaSubnet <- RunBLOBFISH(geneSet = intersect(allGenes, meningioma$gene), 
                                pValueFile = paste0(resultDirectory, "recurrentGBMPval.RDS"),loadPValues = FALSE,
                                networks = meningioma, alpha = 0.005, hopConstraint = 2, nullDistribution = pandas,
                                verbose = TRUE)
write.csv(meningiomaSubnet, paste0(resultDirectory, "recurrentGBM.csv"))