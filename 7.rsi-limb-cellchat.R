library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratObject)
options(stringsAsFactors = FALSE)

x<-readRDS("rsi_3-split-250.rds")

#compare f &h
Idents(x)<-x$limb
f<-subset(x, idents="F")
Idents(f)<-f$cell.type3
h<-subset(x, idents="H")
Idents(h)<-h$cell.type3
levels(f)<-c("PNMP", "EMP", "AMP", "PDMP", "MMP", "ZMP", "CP", "Chondrocyte", "OP", "Osteoblast",
             "PE", "Epithelia", "Muscle", "Endothelia", "Neuron", "BM")
levels(h)<-c("PNMP", "EMP", "AMP", "PDMP", "MMP", "ZMP", "CP", "Chondrocyte", "OP", "Osteoblast",
             "PE", "Epithelia", "Muscle", "Endothelia", "Neuron", "BM")

data.input.f <- GetAssayData(f, assay = "RNA", slot = "data") 
labels.f <- Idents(f)
meta.f <- data.frame(group = labels.f, row.names = names(labels.f)) 
cellchat.f <- createCellChat(object = data.input.f, meta = meta.f, group.by = "group")
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
cellchat.f@DB <- CellChatDB
cellchat.f <- subsetData(cellchat.f) 
cellchat.f <- identifyOverExpressedGenes(cellchat.f)
cellchat.f <- identifyOverExpressedInteractions(cellchat.f)
ellchat.f <- computeCommunProb(cellchat.f)
cellchat.f <- filterCommunication(cellchat.f, min.cells = 10)
cellchat.f <- computeCommunProbPathway(cellchat.f)
cellchat.f <- aggregateNet(cellchat.f)
cellchat.f <- netAnalysis_computeCentrality(cellchat.f, slot.name = "netP")
saveRDS(cellchat.f, "cellchat-f.rds")

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.f, pattern = "outgoing",slot.name = "netP", signaling = pathways.show)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.f, pattern = "incoming", signaling = pathways.show)
ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.f, pattern = "outgoing",slot.name = "netP")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.f, pattern = "incoming", signaling = pathways.show)
ht1 + ht2

##calculate contribution score from PDMP
mat <- cellchat.f@net$weight
par(mfrow = c(4,4), xpd=TRUE)
groupSize <- as.numeric(table(cellchat.f@idents))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])
}
netVisual_bubble(cellchat.f, sources.use = 4, targets.use = c(1:10), remove.isolate = FALSE)

data.input.h <- GetAssayData(h, assay = "RNA", slot = "data")
labels.h <- Idents(h)
meta.h <- data.frame(group = labels.h, row.names = names(labels.h)) 
cellchat.h <- createCellChat(object = data.input.h, meta = meta.h, group.by = "group")
cellchat.h@DB <- CellChatDB
cellchat.h <- subsetData(cellchat.h) 
cellchat.h <- identifyOverExpressedGenes(cellchat.h)
cellchat.h <- identifyOverExpressedInteractions(cellchat.h)
cellchat.h <- computeCommunProb(cellchat.h)
cellchat.h <- filterCommunication(cellchat.h, min.cells = 10)
cellchat.h <- computeCommunProbPathway(cellchat.h)
cellchat.h <- aggregateNet(cellchat.h)
cellchat.h <- netAnalysis_computeCentrality(cellchat.h, slot.name = "netP")
saveRDS(cellchat.h, "cellchat-h.rds")

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.h, pattern = "outgoing", slot.name = "netP")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.h, pattern = "incoming")
ht1 + ht2

##m: same pipeline as F and H

cellchat.f <- readRDS("cellchat-f.rds")
cellchat.h <- readRDS("cellchat-h.rds")
cellchat.m <- readRDS("cellchat-m.rds")

cellchat.f <- updateCellChat(cellchat.f)
cellchat.h <- updateCellChat(cellchat.h)
cellchat.m <- updateCellChat(cellchat.m)

object.list <- list(F = cellchat.f, H = cellchat.h,M = cellchat.m)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

gg1 <- compareInteractions(cellchat, show.legend = F)
gg2 <- compareInteractions(cellchat, show.legend = F, measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c(1,3))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], 
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}



gg1 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, 
               show.raw = T,
               #cutoff.pvalue = 0.01, 
               comparison = c(1,2))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, 
               show.raw = T,
               #cutoff.pvalue = 0.01, 
               comparison = c(1,3))
gg1 + gg2

source("I:/works/whatever/manual/codes/r/cellchat-functional.R")
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

##DEG between F and H (same pipeline as between F and M)
pos.dataset = "F"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, 
                                       features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, 
                                       thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "F",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "H",ligand.logFC = -0.1, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairF.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairF.use.up, 
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
pairF.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairF.use.down,
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

# Chord diagram for all cluster DEG between F and H
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[1]], slot.name = 'net', 
                     net = net.up, lab.cex = 0.8, small.gap = 3.5, 
                     #sources.use = c(4), targets.use = c(1:16),
                     #reduce = 0.003,
                     #directional = 1,
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
netVisual_chord_gene(object.list[[2]], slot.name = 'net', 
                     net = net.down, lab.cex = 0.8, small.gap = 3.5, 
                     #sources.use = c(4), targets.use = c(1:16),
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))


##Signaling IN and OUT ALL
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, 
                                        slot.name = "netP",
                                        title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, 
                                        slot.name = "netP",
                                        title = names(object.list)[i+1], width = 5, height = 6)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, 
                                        slot.name = "netP",
                                        title = names(object.list)[i+2], width = 5, height = 6)
draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(object.list)[i+2], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))


