# Find overlapping genes.
if(!require("ComplexHeatmap")){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("ComplexHeatmap")
}
library("ComplexHeatmap")
if(!require("grid")){
  install.packages("grid")
}
library("grid")
library("netZooR")

# Read in the networks.
networkDir <- NULL
combMatDir <- NULL
controlNet <- read.csv(paste0(networkDir, "adultControl.csv"), row.names = 1)
fetalCNS <- read.csv(paste0(networkDir, "fetalCNS.csv"), row.names = 1)
gbm <- read.csv(paste0(networkDir, "gbm.csv"), row.names = 1)
lgg <- read.csv(paste0(networkDir, "lgg.csv"), row.names = 1)
met <- read.csv(paste0(networkDir, "lungMetastasis.csv"), row.names = 1)
controlNet$score <- 1
fetalCNS$score <- 1
gbm$score <- 1
lgg$score <- 1
met$score <- 1

# Compute targeting with cell adhesion.
CrossNetworkOnly <- function(network, list1, list2){
  networkInList <- network[which(network$gene %in% c(list1, list2)),]
  networkAcrossListsList <- lapply(unique(networkInList$tf), function(tf){
    retval <- data.frame(tf = c(), gene = c())
    targetsInList1 <- intersect(networkInList[which(networkInList$tf == tf), "gene"], list1)
    targetsInList2 <- intersect(networkInList[which(networkInList$tf == tf), "gene"], list2)
    if(length(targetsInList1) > 0 && length(targetsInList2) > 0){
      retval <- networkInList[which(networkInList$tf == tf),]
    }
    return(retval)
  })
  networkAcrossLists <- do.call(rbind, networkAcrossListsList)
  return(networkAcrossLists)
}

# Plot the networks.
oncoFetalEndothelialMolecules <- c("SOX4", "CD93", "COL4A2", "CALCRL", "HSP90B1", "COL4A1", "AKAP13", "HSPG2")
oncoFetalEndothelialCellAdhesionMolecules <- c("ICAM1", "ICAM2", "ICAM3", "VCAM1", "SELP", "SELE", "MADCAM1",  "CD34", "GLYCAM1")
oncoFetalEndothelialImmunosuppressiveMolecules <- c("CLDN5", "JAM2", "TJP1", "LGALS3", "TNFSF10", "INSR", "PECAM1", "STAB1", "VWF",
                                                    "B2M")
allGenes <- c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules,
              oncoFetalEndothelialImmunosuppressiveMolecules)
geneColorMapping <- data.frame(gene = allGenes, color = c(rep(rgb(200 / 255,0 / 255,0 / 255,255 / 255), length(oncoFetalEndothelialMolecules)),
                                                          rep(rgb(0 / 255,175 / 255,255 / 255,255 / 255), length(oncoFetalEndothelialCellAdhesionMolecules)),
                                                          rep(rgb(150 / 255,200 / 255,255 / 255,255 / 255), length(oncoFetalEndothelialImmunosuppressiveMolecules))))
PlotNetworkLocal <- function(network, genesOfInterest,
                             tfColor = "palegreen4", nodeSize = 1, emphasizedVertexLabelSize = 0.25, emphasizedVertexSize = 5,
                             edgeWidth = 0.5, vertexLabels = NA, vertexLabelSize = 0.1, nodesToEmphasize = "SOX4",
                             vertexLabelOffset = 0.5, layoutBipartite = TRUE, geneColorMapping = NULL, makeThicker = NULL,
                             makeEdgesBlack = NULL){
  # Set the node attributes.
  uniqueNodes <- unique(c(network$tf, network$gene))
  str(uniqueNodes)
  nodeAttrs <- data.frame(node = uniqueNodes,
                          color = rep("gray", length(uniqueNodes)),
                          size = rep(nodeSize, length(uniqueNodes)),
                          label.color = "black", label.cex = vertexLabelSize,
                          label.dist = vertexLabelOffset)
  rownames(nodeAttrs) <- uniqueNodes
  str(nodeAttrs)
  
  # Add TF colors.
  nodeAttrs[which(uniqueNodes %in% network$tf), "color"] <- tfColor
  nodeAttrs[which(uniqueNodes %in% network$tf), "frame.color"] <- tfColor
  str(nodeAttrs)
  
  # Add gene colors and change the sizes of the labels.
  rownames(geneColorMapping) <- geneColorMapping$gene
  geneColorMapping <- geneColorMapping[intersect(rownames(geneColorMapping), uniqueNodes),]
  nodeAttrs[which(uniqueNodes %in% network$gene), "size"] <- emphasizedVertexSize
  nodeAttrs[which(uniqueNodes %in% network$gene), "label.cex"] <- emphasizedVertexLabelSize
  if(!is.null(geneColorMapping)){
    rgbValGene <- col2rgb(geneColorMapping$color)
    nodeAttrs[rownames(geneColorMapping), "color"] <- unlist(lapply(1:ncol(rgbValGene), function(i){
      return(rgb(red=min(rgbValGene[1, i], 255), 
                 green=min(rgbValGene[2, i], 255),
                 blue=min(rgbValGene[3, i], 225),
                 alpha=100, maxColorValue=255))
    }))
    print(nodeAttrs[rownames(geneColorMapping), "color"])
    nodeAttrs[rownames(geneColorMapping), "frame.color"] <- geneColorMapping$color
    print(nodeAttrs[rownames(geneColorMapping), "frame.color"])
  }
  str(geneColorMapping)
  
  # Modify nodes to emphasize.
  nodeAttrs[which(uniqueNodes %in% nodesToEmphasize), "label.cex"] <- emphasizedVertexLabelSize
  nodeAttrs[, "label.font"] <- 2
  rgbValTF <- col2rgb(tfColor)
  nodeAttrs[which(uniqueNodes %in% nodesToEmphasize), "color"] <- rgb(red=min(rgbValTF[1,] * 2, 255), 
                                                                      green=min(rgbValTF[2,] * 2, 255),
                                                                      blue=min(rgbValTF[3,] * 2, 225),
                                                                      alpha=125, maxColorValue=255)
  nodeAttrs[which(uniqueNodes %in% nodesToEmphasize), "size"] <- emphasizedVertexSize
  nodeAttrs[, "frame.width"] <- 1
  nodeAttrs[which(uniqueNodes %in% makeThicker), "frame.width"] <- 3
  nodeAttrs[which(uniqueNodes %in% makeThicker), "size"] <- emphasizedVertexSize * 2
  nodeAttrs[which(uniqueNodes %in% makeThicker), "label.cex"] <- emphasizedVertexLabelSize * 1.5
  str(nodeAttrs)
  
  # Add edge attributes.
  if(!is.null(geneColorMapping)){
    for(gene in rownames(geneColorMapping)){
      network[which(network$gene == gene), "color"] <- geneColorMapping[gene, "color"]
    }
  }
  network$width <- edgeWidth
  network[which(network$gene %in% makeEdgesBlack), "width"] <- edgeWidth * 2
  network[which(network$gene %in% makeEdgesBlack), "color"] <- "black"
  str(network)
  
  # Create a graph object.
  graph <- igraph::graph_from_data_frame(network, vertices = nodeAttrs, directed = FALSE)
  V(graph)$type <- V(graph)$name %in% network$tf
  str(graph)
  
  # Plot.
  labels <- V(graph)$name
  whichEmpty <- which(labels %in% setdiff(labels, vertexLabels))
  labels[whichEmpty] <- rep(NA, length(whichEmpty))
  
  if(layoutBipartite == TRUE){
    LO <- layout_as_bipartite(graph)
    LO <- LO[,c(2,1)]
    igraph::plot.igraph(graph, layout = LO, vertex.label = labels)
  }else{
    str(V(graph))
    LO <- layout_with_kk(graph)
    igraph::plot.igraph(graph, layout = LO, vertex.label = labels)
  }
}


# Now, filter the networks to only include SOX4 edges.
fetalCNS_SOX4 <- fetalCNS[which(fetalCNS$tf == "SOX4"),]
controlNet_SOX4 <- controlNet[which(controlNet$tf == "SOX4"),]
lgg_SOX4 <- lgg[which(lgg$tf == "SOX4"),]
gbm_SOX4 <- gbm[which(gbm$tf == "SOX4"),]
met_SOX4 <- met[which(met$tf == "SOX4"),]
dev.off()
par(mfrow = c(3,5))
par(mar=c(0,0,0,0)+.1)
PlotNetworkLocal(network = fetalCNS_SOX4[which(fetalCNS_SOX4$gene %in% oncoFetalEndothelialMolecules),],
                 genesOfInterest = c(fetalCNS_SOX4[which(fetalCNS_SOX4$gene %in% oncoFetalEndothelialMolecules), "tf"], 
                                     fetalCNS_SOX4[which(fetalCNS_SOX4$gene %in% oncoFetalEndothelialMolecules), "gene"]), 
                 vertexLabelOffset = 0,
                 vertexLabels = c(fetalCNS_SOX4$tf, fetalCNS_SOX4$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.5, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = controlNet_SOX4[which(controlNet_SOX4$gene %in% oncoFetalEndothelialMolecules),],
                 genesOfInterest = c(controlNet_SOX4$tf, controlNet_SOX4$gene), 
                 vertexLabelOffset = 0,
                 vertexLabels = c(controlNet_SOX4$tf, controlNet_SOX4$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.5, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = lgg_SOX4[which(lgg_SOX4$gene %in% oncoFetalEndothelialMolecules),],
                 genesOfInterest = c(lgg_SOX4$tf, lgg_SOX4$gene), 
                 vertexLabelOffset = 0,
                 vertexLabels = c(lgg_SOX4$tf, lgg_SOX4$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.5, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = gbm_SOX4[which(gbm_SOX4$gene %in% oncoFetalEndothelialMolecules),],
                 genesOfInterest = c(gbm_SOX4$tf, gbm_SOX4$gene), 
                 vertexLabelOffset = 0,
                 vertexLabels = c(gbm_SOX4$tf, gbm_SOX4$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.5, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = met_SOX4[which(met_SOX4$gene %in% oncoFetalEndothelialMolecules),],
                 genesOfInterest = c(met_SOX4$tf, met_SOX4$gene), 
                 vertexLabelOffset = 0,
                 vertexLabels = c(met_SOX4$tf, met_SOX4$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.5, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = CrossNetworkOnly(fetalCNS_SOX4, oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules),
                 genesOfInterest = c(fetalCNS_SOX4[, "tf"], fetalCNS_SOX4[, "gene"]), 
                 vertexLabelOffset = 0,
                 vertexLabels = c(fetalCNS_SOX4$tf, fetalCNS_SOX4$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.5, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = CrossNetworkOnly(controlNet_SOX4, oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules),
                 genesOfInterest = c(controlNet_SOX4$tf, controlNet_SOX4$gene), 
                 vertexLabelOffset = 0,
                 vertexLabels = c(controlNet_SOX4$tf, controlNet_SOX4$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.5, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = CrossNetworkOnly(lgg_SOX4, oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules),
                 genesOfInterest = c(lgg_SOX4$tf, lgg_SOX4$gene), 
                 vertexLabelOffset = 0,
                 vertexLabels = c(lgg_SOX4$tf, lgg_SOX4$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.5, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = CrossNetworkOnly(gbm_SOX4, oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules),
                 genesOfInterest = c(gbm_SOX4$tf, gbm_SOX4$gene), 
                 vertexLabelOffset = 0,
                 vertexLabels = c(gbm_SOX4$tf, gbm_SOX4$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.5, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = CrossNetworkOnly(met_SOX4, oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules),
                 genesOfInterest = c(met_SOX4$tf, met_SOX4$gene), 
                 vertexLabelOffset = 0,
                 vertexLabels = c(met_SOX4$tf, met_SOX4$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.5, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = CrossNetworkOnly(fetalCNS_SOX4, oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules),
                 genesOfInterest = c(fetalCNS_SOX4[, "tf"], fetalCNS_SOX4[, "gene"]), 
                 vertexLabelOffset = 0,
                 vertexLabels = c(fetalCNS_SOX4$tf, fetalCNS_SOX4$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.5, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = CrossNetworkOnly(controlNet_SOX4, oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules),
                 genesOfInterest = c(controlNet_SOX4$tf, controlNet_SOX4$gene), 
                 vertexLabelOffset = 0,
                 vertexLabels = c(controlNet_SOX4$tf, controlNet_SOX4$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.5, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = CrossNetworkOnly(lgg_SOX4, oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules),
                 genesOfInterest = c(lgg_SOX4$tf, lgg_SOX4$gene), 
                 vertexLabelOffset = 0,
                 vertexLabels = c(lgg_SOX4$tf, lgg_SOX4$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.5, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = CrossNetworkOnly(gbm_SOX4, oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules),
                 genesOfInterest = c(gbm_SOX4$tf, gbm_SOX4$gene), 
                 vertexLabelOffset = 0,
                 vertexLabels = c(gbm_SOX4$tf, gbm_SOX4$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.5, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = CrossNetworkOnly(met_SOX4, oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules),
                 genesOfInterest = c(met_SOX4$tf, met_SOX4$gene), 
                 vertexLabelOffset = 0,
                 vertexLabels = c(met_SOX4$tf, met_SOX4$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.5, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)

# # Combine the networks.
# edgeList <- unique(c(rownames(gbm), rownames(lgg), rownames(controlNet), rownames(fetalCNS),
#                      rownames(met)))
# edgeListTF <- unlist(lapply(edgeList, function(edge){return(strsplit(edge, "__")[[1]][1])}))
# edgeListGene <- unlist(lapply(edgeList, function(edge){return(strsplit(edge, "__")[[1]][2])}))
# edgeListDF <- data.frame(tf = edgeListTF, gene = edgeListGene, score = 0)
# rownames(edgeListDF) <- edgeList
# controlNet <- rbind(controlNet, edgeListDF[setdiff(rownames(edgeListDF), rownames(controlNet)),])
# fetalCNS <- rbind(fetalCNS, edgeListDF[setdiff(rownames(edgeListDF), rownames(fetalCNS)),])
# gbm <- rbind(gbm, edgeListDF[setdiff(rownames(edgeListDF), rownames(gbm)),])
# lgg <- rbind(lgg, edgeListDF[setdiff(rownames(edgeListDF), rownames(lgg)),])
# met <- rbind(met, edgeListDF[setdiff(rownames(edgeListDF), rownames(met)),])

# Check the overlap in transcription factors.
HubTFsPerc <- function(network, nodesToEmphasizeCutoffPercentile){
  degreeCentrality <- table(network$tf)
  xthQuantile <- quantile(degreeCentrality, probs = nodesToEmphasizeCutoffPercentile)
  print(xthQuantile)
  return(names(degreeCentrality)[which(degreeCentrality >= xthQuantile)])
}

# Check only for the oncofetal genes.
upFetal <- HubTFsPerc(fetalCNS[which(fetalCNS$gene %in% oncoFetalEndothelialMolecules),], 0.99)
upLGG <- HubTFsPerc(lgg[which(lgg$gene %in% oncoFetalEndothelialMolecules),], 0.99)
upGBM <- HubTFsPerc(gbm[which(gbm$gene %in% oncoFetalEndothelialMolecules),], 0.99)
upMET <- HubTFsPerc(met[which(met$gene %in% oncoFetalEndothelialMolecules),], 0.99)
upControl <- HubTFsPerc(controlNet[which(controlNet$gene %in% oncoFetalEndothelialMolecules),], 0.99)

tfUpList = list(FetalBrain = upFetal, LGG = upLGG, GBM = upGBM,
                Metastasis = upMET, AdultControls = upControl)
combMat = ComplexHeatmap::make_comb_mat(tfUpList)
ss <- ComplexHeatmap::set_size(combMat)
cs <- ComplexHeatmap::comb_size(combMat)
cd <- ComplexHeatmap::comb_degree(combMat)
ht <- ComplexHeatmap::UpSet(combMat, 
                            set_order = c("FetalBrain", "AdultControls", "LGG", 
                                          "GBM", "Metastasis"),
                            comb_order = c(3, 1, 4, 2, 6, 10,
                                           14, 8, 5, 9,
                                           7,
                                           12, 11,
                                           13),
                            top_annotation = ComplexHeatmap::HeatmapAnnotation(
                              "Shared Hub TFs\nAmong Core Onco\n-Fetal Signature" = ComplexHeatmap::anno_barplot(cs, 
                                                                              ylim = c(0, max(cs)*1.1),
                                                                              border = FALSE, 
                                                                              gp = gpar(fill = "black"), 
                                                                              height = grid::unit(4, "cm"),
                              ), 
                              annotation_name_side = "left", 
                              annotation_name_rot = 90
                            ),
                            right_annotation = ComplexHeatmap::rowAnnotation(
                              "Hub TFs" = ComplexHeatmap::anno_barplot(ss, 
                                                                       baseline = 0,
                                                                       border = FALSE, 
                                                                       gp = gpar(fill = "black"), 
                                                                       width = grid::unit(4, "cm")
                              )
                            ),
                            left_annotation = NULL,
                            show_row_names = TRUE
)
ht
od = ComplexHeatmap::column_order(ht)
ComplexHeatmap::decorate_annotation("Shared Hub TFs\nAmong Core Onco\n-Fetal Signature", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 10, col = "black"), rot = 45)
})
write.csv(attr(combMat, "data"), paste0(combMatDir, "/tfHubMatOncofetal.csv"))

# Check for oncofetal + cell adhesion. We note that the 99th percentile is 6 for 
# GBM, Control, and LGG
upFetalCA <- intersect(upFetal, HubTFsPerc(fetalCNS[which(fetalCNS$gene %in% oncoFetalEndothelialCellAdhesionMolecules),], 0.99))
upLGGCA <- intersect(upLGG, HubTFsPerc(lgg[which(lgg$gene %in% oncoFetalEndothelialCellAdhesionMolecules),], 0.99))
upGBMCA <- intersect(upGBM, HubTFsPerc(gbm[which(gbm$gene %in% oncoFetalEndothelialCellAdhesionMolecules),], 0.99))
upMETCA <- intersect(upMET, HubTFsPerc(met[which(met$gene %in% oncoFetalEndothelialCellAdhesionMolecules),], 0.99))
upControlCA <- intersect(upControl, HubTFsPerc(controlNet[which(controlNet$gene %in% oncoFetalEndothelialCellAdhesionMolecules),], 0.99))

HubTFs <- function(network, hardCutoff){
  degreeCentrality <- table(network$tf)
  return(names(degreeCentrality)[which(degreeCentrality >= hardCutoff)])
}
upFetalCA <- intersect(upFetal, HubTFs(fetalCNS[which(fetalCNS$gene %in% oncoFetalEndothelialCellAdhesionMolecules),], 6))
upLGGCA <- intersect(upLGG, HubTFs(lgg[which(lgg$gene %in% oncoFetalEndothelialCellAdhesionMolecules),], 6))
upGBMCA <- intersect(upGBM, HubTFs(gbm[which(gbm$gene %in% oncoFetalEndothelialCellAdhesionMolecules),], 6))
upMETCA <- intersect(upMET, HubTFs(met[which(met$gene %in% oncoFetalEndothelialCellAdhesionMolecules),], 6))
upControlCA <- intersect(upControl, HubTFs(controlNet[which(controlNet$gene %in% oncoFetalEndothelialCellAdhesionMolecules),], 6))

tfUpList = list(FetalBrain = upFetalCA, LGG = upLGGCA, GBM = upGBMCA,
                Metastasis = upMETCA, AdultControls = upControlCA)
combMat = ComplexHeatmap::make_comb_mat(tfUpList)
ss <- ComplexHeatmap::set_size(combMat)
cs <- ComplexHeatmap::comb_size(combMat)
cd <- ComplexHeatmap::comb_degree(combMat)
ht <- ComplexHeatmap::UpSet(combMat, 
                            set_order = c("FetalBrain", "AdultControls", "LGG", 
                                          "GBM", "Metastasis"),
                            comb_order = c(10, 1, 3, 5, 2, 4, 7, 6, 9, 8, 11, 12),
                            top_annotation = ComplexHeatmap::HeatmapAnnotation(
                              "Shared Hub TFs\nConnecting Core\nOnco-Fetal and Cell\nAdhesion Signatures" = ComplexHeatmap::anno_barplot(cs, 
                                                                              ylim = c(0, max(cs)*1.1),
                                                                              border = FALSE, 
                                                                              gp = gpar(fill = "black"), 
                                                                              height = grid::unit(4, "cm"),
                              ), 
                              annotation_name_side = "left", 
                              annotation_name_rot = 90
                            ),
                            right_annotation = ComplexHeatmap::rowAnnotation(
                              "Hub TFs" = ComplexHeatmap::anno_barplot(ss, 
                                                                       baseline = 0,
                                                                       border = FALSE, 
                                                                       gp = gpar(fill = "black"), 
                                                                       width = grid::unit(4, "cm")
                              )
                            ),
                            left_annotation = NULL,
                            show_row_names = TRUE
)
ht
od = ComplexHeatmap::column_order(ht)
ComplexHeatmap::decorate_annotation("Shared Hub TFs\nConnecting Core\nOnco-Fetal and Cell\nAdhesion Signatures", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 10, col = "black"), rot = 45)
})
write.csv(attr(combMat, "data"), paste0(combMatDir, "/tfHubMatCellAdhesion.csv"))

# Check for oncofetal + immunosuppressive.
upFetalIm <- intersect(upFetal, HubTFsPerc(fetalCNS[which(fetalCNS$gene %in% oncoFetalEndothelialImmunosuppressiveMolecules),], 0.99))
upLGGIm <- intersect(upLGG, HubTFsPerc(lgg[which(lgg$gene %in% oncoFetalEndothelialImmunosuppressiveMolecules),], 0.99))
upGBMIm <- intersect(upGBM, HubTFsPerc(gbm[which(gbm$gene %in% oncoFetalEndothelialImmunosuppressiveMolecules),], 0.99))
upMETIm <- intersect(upMET, HubTFsPerc(met[which(met$gene %in% oncoFetalEndothelialImmunosuppressiveMolecules),], 0.99))
upControlIm <- intersect(upControl, HubTFsPerc(controlNet[which(controlNet$gene %in% oncoFetalEndothelialImmunosuppressiveMolecules),], 0.99))

upFetalIm <- intersect(upFetal, HubTFs(fetalCNS[which(fetalCNS$gene %in% oncoFetalEndothelialImmunosuppressiveMolecules),], 8))
upLGGIm <- intersect(upLGG, HubTFs(lgg[which(lgg$gene %in% oncoFetalEndothelialImmunosuppressiveMolecules),], 8))
upGBMIm <- intersect(upGBM, HubTFs(gbm[which(gbm$gene %in% oncoFetalEndothelialImmunosuppressiveMolecules),], 8))
upMETIm <- intersect(upMET, HubTFs(met[which(met$gene %in% oncoFetalEndothelialImmunosuppressiveMolecules),], 8))
upControlIm <- intersect(upControl, HubTFs(controlNet[which(controlNet$gene %in% oncoFetalEndothelialImmunosuppressiveMolecules),], 8))

tfUpList = list(FetalBrain = upFetalIm, LGG = upLGGIm, GBM = upGBMIm,
                Metastasis = upMETIm, AdultControls = upControlIm)
combMat = ComplexHeatmap::make_comb_mat(tfUpList)
ss <- ComplexHeatmap::set_size(combMat)
cs <- ComplexHeatmap::comb_size(combMat)
cd <- ComplexHeatmap::comb_degree(combMat)
ht <- ComplexHeatmap::UpSet(combMat, 
                            set_order = c("FetalBrain", "AdultControls", "LGG", 
                                          "GBM", "Metastasis"),
                            comb_order = c(9, 1, 5,
                                           10, 6, 2, 3, 7,
                                           11, 4,
                                           12, 8,
                                           13),
                            top_annotation = ComplexHeatmap::HeatmapAnnotation(
                              "Shared Hub TFs\nConnecting Core Onco\n-Fetal and\nImmunosuppressive\nSignatures" = ComplexHeatmap::anno_barplot(cs, 
                                                                                                                                         ylim = c(0, max(cs)*1.1),
                                                                                                                                         border = FALSE, 
                                                                                                                                         gp = gpar(fill = "black"), 
                                                                                                                                         height = grid::unit(4, "cm"),
                              ), 
                              annotation_name_side = "left", 
                              annotation_name_rot = 90
                            ),
                            right_annotation = ComplexHeatmap::rowAnnotation(
                              "Hub TFs" = ComplexHeatmap::anno_barplot(ss, 
                                                                       baseline = 0,
                                                                       border = FALSE, 
                                                                       gp = gpar(fill = "black"), 
                                                                       width = grid::unit(4, "cm")
                              )
                            ),
                            left_annotation = NULL,
                            show_row_names = TRUE
)
ht
od = ComplexHeatmap::column_order(ht)
ComplexHeatmap::decorate_annotation("Shared Hub TFs\nConnecting Core Onco\n-Fetal and\nImmunosuppressive\nSignatures", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 10, col = "black"), rot = 45)
})
write.csv(attr(combMat, "data"), paste0(combMatDir, "/tfHubMatImmunosuppression.csv"))

# Plot the subset.
nodesOfImportanceCellAdhesion <- c(Reduce(intersect, list(upFetalCA, upLGGCA, upGBMCA, upMETCA, upControlCA)),
                                   setdiff(Reduce(intersect, list(upFetalCA, upGBMCA, upMETCA)), union(upLGGCA, upControlCA)),
                                   setdiff(Reduce(intersect, list(upFetalCA, upGBMCA, upMETCA, upLGGCA)), upControlCA))
nodesOfImportanceImmunosuppression <- setdiff(Reduce(intersect, list(upFetalIm, upGBMIm)), c(upLGGIm, upControlIm, upMETIm))
nodesOfImportance <- c(nodesOfImportanceCellAdhesion, nodesOfImportanceImmunosuppression)
nodesOfImportanceBoth <- intersect(nodesOfImportanceCellAdhesion, nodesOfImportanceImmunosuppression)

# Emphasize nodes.
dev.off()
par(mfrow = c(5,5))
par(mar=c(0,0,0,0)+.1)
PlotNetworkLocal(network = fetalCNS[intersect(which(fetalCNS$tf %in% c(nodesOfImportance, "SOX4")),
                                              which(fetalCNS$gene %in% oncoFetalEndothelialMolecules)),],
                 genesOfInterest = oncoFetalEndothelialMolecules,
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportance, "SOX4"),
                 vertexLabels = c(fetalCNS$tf, fetalCNS$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "SOX4", "AKAP13"), makeEdgesBlack = c("HSPG2", "AKAP13"))
PlotNetworkLocal(network = controlNet[intersect(which(controlNet$tf %in% c(nodesOfImportance, "SOX4")),
                                                which(controlNet$gene %in% oncoFetalEndothelialMolecules)),],
                 genesOfInterest = oncoFetalEndothelialMolecules,
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportance, "SOX4"),
                 vertexLabels = c(controlNet$tf, controlNet$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "SOX4", "AKAP13"), makeEdgesBlack = c("HSPG2", "AKAP13"))
PlotNetworkLocal(network = lgg[intersect(which(lgg$tf %in% c(nodesOfImportance, "SOX4")),
                                         which(lgg$gene %in% oncoFetalEndothelialMolecules)),],
                 genesOfInterest = oncoFetalEndothelialMolecules,
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportance, "SOX4"),
                 vertexLabels = c(lgg$tf, lgg$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "SOX4", "AKAP13"), makeEdgesBlack = c("HSPG2", "AKAP13"))
PlotNetworkLocal(network = gbm[intersect(which(gbm$tf %in% c(nodesOfImportance, "SOX4")),
                                         which(gbm$gene %in% oncoFetalEndothelialMolecules)),],
                 genesOfInterest = oncoFetalEndothelialMolecules,
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportance, "SOX4"),
                 vertexLabels = c(gbm$tf, gbm$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "SOX4", "AKAP13"), makeEdgesBlack = c("HSPG2", "AKAP13"))
PlotNetworkLocal(network = met[intersect(which(met$tf %in% c(nodesOfImportance, "SOX4")),
                                         which(met$gene %in% oncoFetalEndothelialMolecules)),],
                 genesOfInterest = oncoFetalEndothelialMolecules,
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportance, "SOX4"),
                 vertexLabels = c(met$tf, met$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "SOX4", "AKAP13"), makeEdgesBlack = c("HSPG2", "AKAP13"))
PlotNetworkLocal(network = fetalCNS[intersect(which(fetalCNS$tf %in% c(nodesOfImportanceCellAdhesion, "SOX4")),
                                              which(fetalCNS$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(fetalCNS$tf, fetalCNS$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "SOX4", "AKAP13"), makeEdgesBlack = c("HSPG2", "AKAP13"))
PlotNetworkLocal(network = controlNet[intersect(which(controlNet$tf %in% c(nodesOfImportanceCellAdhesion, "SOX4")),
                                                which(controlNet$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(controlNet$tf, controlNet$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "SOX4", "AKAP13"), makeEdgesBlack = c("HSPG2", "AKAP13"))
PlotNetworkLocal(network = lgg[intersect(which(lgg$tf %in% c(nodesOfImportanceCellAdhesion, "SOX4")),
                                         which(lgg$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(lgg$tf, lgg$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "SOX4", "AKAP13"), makeEdgesBlack = c("HSPG2", "AKAP13"))
PlotNetworkLocal(network = gbm[intersect(which(gbm$tf %in% c(nodesOfImportanceCellAdhesion, "SOX4")),
                                         which(gbm$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(gbm$tf, gbm$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "SOX4", "AKAP13"), makeEdgesBlack = c("HSPG2", "AKAP13"))
PlotNetworkLocal(network = met[intersect(which(met$tf %in% c(nodesOfImportanceCellAdhesion, "SOX4")),
                                         which(met$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(met$tf, met$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "SOX4", "AKAP13"), makeEdgesBlack = c("HSPG2", "AKAP13"))
PlotNetworkLocal(network = fetalCNS[intersect(which(fetalCNS$tf %in% c(nodesOfImportanceImmunosuppression, "SOX4")),
                                              which(fetalCNS$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, "SOX4"),
                 vertexLabels = c(fetalCNS$tf, fetalCNS$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "STAB1", "SOX4", "AKAP13"), 
                 makeEdgesBlack = c("HSPG2", "STAB1", "AKAP13"))
PlotNetworkLocal(network = controlNet[intersect(which(controlNet$tf %in% c(nodesOfImportanceImmunosuppression, "SOX4")),
                                                which(controlNet$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, "SOX4"),
                 vertexLabels = c(controlNet$tf, controlNet$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "STAB1", "SOX4", "AKAP13"), 
                 makeEdgesBlack = c("HSPG2", "STAB1", "AKAP13"))
PlotNetworkLocal(network = lgg[intersect(which(lgg$tf %in% c(nodesOfImportanceImmunosuppression, "SOX4")),
                                         which(lgg$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, "SOX4"),
                 vertexLabels = c(lgg$tf, lgg$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "STAB1", "SOX4", "AKAP13"), 
                 makeEdgesBlack = c("HSPG2", "STAB1", "AKAP13"))
PlotNetworkLocal(network = gbm[intersect(which(gbm$tf %in% c(nodesOfImportanceImmunosuppression, "SOX4")),
                                         which(gbm$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, "SOX4"),
                 vertexLabels = c(gbm$tf, gbm$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "STAB1", "SOX4", "AKAP13"), 
                 makeEdgesBlack = c("HSPG2", "STAB1", "AKAP13"))
PlotNetworkLocal(network = met[intersect(which(met$tf %in% c(nodesOfImportanceImmunosuppression, "SOX4")),
                                         which(met$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, "SOX4"),
                 vertexLabels = c(met$tf, met$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "STAB1", "SOX4", "AKAP13"), 
                 makeEdgesBlack = c("HSPG2", "STAB1", "AKAP13"))
PlotNetworkLocal(network = fetalCNS[intersect(which(fetalCNS$tf %in% c(nodesOfImportanceBoth, "SOX4")),
                                              which(fetalCNS$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                                         oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(fetalCNS$tf, fetalCNS$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "STAB1", "SOX4", "AKAP13"), 
                 makeEdgesBlack = c("HSPG2", "STAB1", "AKAP13"))
PlotNetworkLocal(network = controlNet[intersect(which(controlNet$tf %in% c(nodesOfImportanceBoth, "SOX4")),
                                                which(controlNet$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                                             oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion,"SOX4"),
                 vertexLabels = c(controlNet$tf, controlNet$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "STAB1", "SOX4", "AKAP13"), 
                 makeEdgesBlack = c("HSPG2", "STAB1", "AKAP13"))
PlotNetworkLocal(network = lgg[intersect(which(lgg$tf %in% c(nodesOfImportanceBoth, "SOX4")),
                                         which(lgg$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                               oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(lgg$tf, lgg$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "STAB1", "SOX4", "AKAP13"), 
                 makeEdgesBlack = c("HSPG2", "STAB1", "AKAP13"))
PlotNetworkLocal(network = gbm[intersect(which(gbm$tf %in% c(nodesOfImportanceBoth, "SOX4")),
                                         which(gbm$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                               oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(gbm$tf, gbm$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "STAB1", "SOX4", "AKAP13"), 
                 makeEdgesBlack = c("HSPG2", "STAB1", "AKAP13"))
PlotNetworkLocal(network = met[intersect(which(met$tf %in% c(nodesOfImportanceBoth, "SOX4")),
                                         which(met$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                               oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(met$tf, met$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "STAB1", "SOX4", "AKAP13"), 
                 makeEdgesBlack = c("HSPG2", "STAB1", "AKAP13"))
PlotNetworkLocal(network = fetalCNS[intersect(which(fetalCNS$tf %in% c(nodesOfImportance, "SOX4")),
                                              which(fetalCNS$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                                         oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(fetalCNS$tf, fetalCNS$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "STAB1", "SOX4", "AKAP13"), 
                 makeEdgesBlack = c("HSPG2", "STAB1", "AKAP13"))
PlotNetworkLocal(network = controlNet[intersect(which(controlNet$tf %in% c(nodesOfImportance, "SOX4")),
                                                which(controlNet$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                                             oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion,"SOX4"),
                 vertexLabels = c(controlNet$tf, controlNet$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "STAB1", "SOX4", "AKAP13"), 
                 makeEdgesBlack = c("HSPG2", "STAB1", "AKAP13"))
PlotNetworkLocal(network = lgg[intersect(which(lgg$tf %in% c(nodesOfImportance, "SOX4")),
                                         which(lgg$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                               oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(lgg$tf, lgg$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "STAB1", "SOX4", "AKAP13"), 
                 makeEdgesBlack = c("HSPG2", "STAB1", "AKAP13"))
PlotNetworkLocal(network = gbm[intersect(which(gbm$tf %in% c(nodesOfImportance, "SOX4")),
                                         which(gbm$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                               oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(gbm$tf, gbm$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "STAB1", "SOX4", "AKAP13"), 
                 makeEdgesBlack = c("HSPG2", "STAB1", "AKAP13"))
PlotNetworkLocal(network = met[intersect(which(met$tf %in% c(nodesOfImportance, "SOX4")),
                                         which(met$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                               oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(met$tf, met$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001, makeThicker = c("HSPG2", "STAB1", "SOX4", "AKAP13"), 
                 makeEdgesBlack = c("HSPG2", "STAB1", "AKAP13"))

# Emphasize nodes.
dev.off()
par(mfrow = c(5,5))
par(mar=c(0,0,0,0)+.1)
PlotNetworkLocal(network = fetalCNS[intersect(which(fetalCNS$tf %in% c(nodesOfImportance, "SOX4")),
                                              which(fetalCNS$gene %in% oncoFetalEndothelialMolecules)),],
                 genesOfInterest = oncoFetalEndothelialMolecules,
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportance, "SOX4"),
                 vertexLabels = c(fetalCNS$tf, fetalCNS$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = controlNet[intersect(which(controlNet$tf %in% c(nodesOfImportance, "SOX4")),
                                                which(controlNet$gene %in% oncoFetalEndothelialMolecules)),],
                 genesOfInterest = oncoFetalEndothelialMolecules,
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportance, "SOX4"),
                 vertexLabels = c(controlNet$tf, controlNet$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = lgg[intersect(which(lgg$tf %in% c(nodesOfImportance, "SOX4")),
                                         which(lgg$gene %in% oncoFetalEndothelialMolecules)),],
                 genesOfInterest = oncoFetalEndothelialMolecules,
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportance, "SOX4"),
                 vertexLabels = c(lgg$tf, lgg$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = gbm[intersect(which(gbm$tf %in% c(nodesOfImportance, "SOX4")),
                                         which(gbm$gene %in% oncoFetalEndothelialMolecules)),],
                 genesOfInterest = oncoFetalEndothelialMolecules,
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportance, "SOX4"),
                 vertexLabels = c(gbm$tf, gbm$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = met[intersect(which(met$tf %in% c(nodesOfImportance, "SOX4")),
                                         which(met$gene %in% oncoFetalEndothelialMolecules)),],
                 genesOfInterest = oncoFetalEndothelialMolecules,
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportance, "SOX4"),
                 vertexLabels = c(met$tf, met$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = fetalCNS[intersect(which(fetalCNS$tf %in% c(nodesOfImportanceCellAdhesion, "SOX4")),
                                              which(fetalCNS$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(fetalCNS$tf, fetalCNS$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = controlNet[intersect(which(controlNet$tf %in% c(nodesOfImportanceCellAdhesion, "SOX4")),
                                                which(controlNet$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(controlNet$tf, controlNet$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = lgg[intersect(which(lgg$tf %in% c(nodesOfImportanceCellAdhesion, "SOX4")),
                                         which(lgg$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(lgg$tf, lgg$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = gbm[intersect(which(gbm$tf %in% c(nodesOfImportanceCellAdhesion, "SOX4")),
                                         which(gbm$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(gbm$tf, gbm$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = met[intersect(which(met$tf %in% c(nodesOfImportanceCellAdhesion, "SOX4")),
                                         which(met$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(met$tf, met$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001)
PlotNetworkLocal(network = fetalCNS[intersect(which(fetalCNS$tf %in% c(nodesOfImportanceImmunosuppression, "SOX4")),
                                              which(fetalCNS$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, "SOX4"),
                 vertexLabels = c(fetalCNS$tf, fetalCNS$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001 
                 )
PlotNetworkLocal(network = controlNet[intersect(which(controlNet$tf %in% c(nodesOfImportanceImmunosuppression, "SOX4")),
                                                which(controlNet$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, "SOX4"),
                 vertexLabels = c(controlNet$tf, controlNet$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001 
                 )
PlotNetworkLocal(network = lgg[intersect(which(lgg$tf %in% c(nodesOfImportanceImmunosuppression, "SOX4")),
                                         which(lgg$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, "SOX4"),
                 vertexLabels = c(lgg$tf, lgg$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001 
                 )
PlotNetworkLocal(network = gbm[intersect(which(gbm$tf %in% c(nodesOfImportanceImmunosuppression, "SOX4")),
                                         which(gbm$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, "SOX4"),
                 vertexLabels = c(gbm$tf, gbm$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001 
                 )
PlotNetworkLocal(network = met[intersect(which(met$tf %in% c(nodesOfImportanceImmunosuppression, "SOX4")),
                                         which(met$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, "SOX4"),
                 vertexLabels = c(met$tf, met$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001 
                 )
PlotNetworkLocal(network = fetalCNS[intersect(which(fetalCNS$tf %in% c(nodesOfImportanceBoth, "SOX4")),
                                              which(fetalCNS$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                                         oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(fetalCNS$tf, fetalCNS$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001 
                 )
PlotNetworkLocal(network = controlNet[intersect(which(controlNet$tf %in% c(nodesOfImportanceBoth, "SOX4")),
                                                which(controlNet$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                                             oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion,"SOX4"),
                 vertexLabels = c(controlNet$tf, controlNet$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001 
                 )
PlotNetworkLocal(network = lgg[intersect(which(lgg$tf %in% c(nodesOfImportanceBoth, "SOX4")),
                                         which(lgg$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                               oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(lgg$tf, lgg$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001 
                 )
PlotNetworkLocal(network = gbm[intersect(which(gbm$tf %in% c(nodesOfImportanceBoth, "SOX4")),
                                         which(gbm$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                               oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(gbm$tf, gbm$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001 
                 )
PlotNetworkLocal(network = met[intersect(which(met$tf %in% c(nodesOfImportanceBoth, "SOX4")),
                                         which(met$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                               oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(met$tf, met$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001 
                 )
PlotNetworkLocal(network = fetalCNS[intersect(which(fetalCNS$tf %in% c(nodesOfImportance, "SOX4")),
                                              which(fetalCNS$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                                         oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(fetalCNS$tf, fetalCNS$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001 
                 )
PlotNetworkLocal(network = controlNet[intersect(which(controlNet$tf %in% c(nodesOfImportance, "SOX4")),
                                                which(controlNet$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                                             oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion,"SOX4"),
                 vertexLabels = c(controlNet$tf, controlNet$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001 
                 )
PlotNetworkLocal(network = lgg[intersect(which(lgg$tf %in% c(nodesOfImportance, "SOX4")),
                                         which(lgg$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                               oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(lgg$tf, lgg$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001 
                 )
PlotNetworkLocal(network = gbm[intersect(which(gbm$tf %in% c(nodesOfImportance, "SOX4")),
                                         which(gbm$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                               oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(gbm$tf, gbm$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001 
                 )
PlotNetworkLocal(network = met[intersect(which(met$tf %in% c(nodesOfImportance, "SOX4")),
                                         which(met$gene %in% c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                                               oncoFetalEndothelialCellAdhesionMolecules))),],
                 genesOfInterest = c(oncoFetalEndothelialMolecules, oncoFetalEndothelialImmunosuppressiveMolecules,
                                     oncoFetalEndothelialCellAdhesionMolecules),
                 vertexLabelOffset = 0, nodesToEmphasize = c(nodesOfImportanceImmunosuppression, 
                                                             nodesOfImportanceCellAdhesion, "SOX4"),
                 vertexLabels = c(met$tf, met$gene), geneColorMapping = geneColorMapping,
                 layoutBipartite = FALSE,edgeWidth = 0.3, emphasizedVertexSize = 20, emphasizedVertexLabelSize = 0.5,
                 vertexLabelSize = 0.0001 
                 )