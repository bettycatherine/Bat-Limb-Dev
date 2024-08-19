library(dplyr)
library(Seurat)
library(patchwork)
library(tibble)
library(Rtsne)
library(sctransform)
library(ggplot2)

x<-readRDS("rsi_3-split-250.rds")
z<-readRDS("rsi_3-10x-250-f.rds")
v<-readRDS("rsi_3-10x-250-h.rds")
DefaultAssay(x)<-"integrated"
DefaultAssay(z)<-"SCT"
DefaultAssay(v)<-"SCT"
levels(x)<-c("PNMP", "EMP", "AMP", "PDMP", "MMP", "ZMP", "CP", "Chondrocyte", "OP", "Osteoblast",
             "PE", "Epithelia", "Muscle", "Endothelia", "Neuron", "BM")

f.anchors <- FindTransferAnchors(reference = x, query = z, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = f.anchors, refdata = x$cell.type3, dims = 1:30)
z <- AddMetaData(z, metadata = predictions)
z$prediction.match <- z$predicted.id
z$cell.type3 <- z$predicted.id
table(z$prediction.match)
table(z$predicted.id)
x <- RunUMAP(x, dims = 1:30, reduction = "pca", return.model = TRUE)
z <- MapQuery(anchorset = f.anchors, reference = x, query = z,
                           refdata = list(celltype = "cell.type3"), 
              reference.reduction = "pca", reduction.model = "umap")
p1 <- DimPlot(x,label = TRUE, label.size = 5,
              repel = TRUE)  + ggtitle("SPLIT-ALL")
p2 <- DimPlot(z, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 5, repel = TRUE) + ggtitle("F21-TL")
p1 + p2

saveRDS(z,"rsi_3-10x-250-f.rds")

h.anchors <- FindTransferAnchors(reference = x, query = v, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = h.anchors, refdata = x$cell.type3, dims = 1:30)
v <- AddMetaData(v, metadata = predictions)
v$prediction.match <- v$predicted.id
v$cell.type3<- v$predicted.id
table(v$prediction.match)
table(v$predicted.id)
x <- RunUMAP(x, dims = 1:30, reduction = "pca", return.model = TRUE)
v <- MapQuery(anchorset = h.anchors, reference = x, query = v,
              refdata = list(celltype = "cell.type3"), 
              reference.reduction = "pca", reduction.model = "umap")
p3 <- DimPlot(x, reduction = "umap", group.by = "cell.type3", label = TRUE, label.size = 5,
              repel = TRUE) + NoLegend() + ggtitle("SPLIT-ALL")
p4 <- DimPlot(v, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 5, repel = TRUE) + ggtitle("H21-TL")
p3+p4
saveRDS(v,"rsi_3-10x-250-h.rds")

m<-readRDS("mm-split-250.rds")
DefaultAssay(m)<-"integrated"
m.anchors <- FindTransferAnchors(reference = x, query = m, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = m.anchors, refdata = x$cell.type3, dims = 1:30)
m <- AddMetaData(m, metadata = predictions)
m$prediction.match <- m$predicted.id
m$cell.type3<-m$predicted.id
table(m$prediction.match)
table(m$predicted.id)
x <- RunUMAP(x, dims = 1:30, reduction = "pca", return.model = TRUE)
m <- MapQuery(anchorset = m.anchors, reference = x, query = m,
              refdata = list(celltype = "cell.type3"), 
              reference.reduction = "pca", reduction.model = "umap")
p5 <- DimPlot(x, reduction = "umap", group.by = "cell.type3",  ncol = 1,label = TRUE, label.size = 5,
              repel = TRUE) + NoLegend() + ggtitle("RSI-ALL")
p6 <- DimPlot(m, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 5, repel = TRUE) + ggtitle("MM-TL")
p5+p6

saveRDS(m, "mm-split-250.rds")


