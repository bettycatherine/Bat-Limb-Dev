library(ggplot2)
library(reshape2)
library(Seurat)
library(SeuratObject)

x<-readRDS("rsi_3-split-250.rds")
Idents(x)<-x$limb
f<-subset(x, idents = "F")
h<-subset(x, idents = "H")
Idents(f)<-f$cell.type3
Idents(h)<-h$cell.type3

DefaultAssay(f)<-"RNA"
DefaultAssay(h)<-"RNA"

gene2cor.f<-intersect(VariableFeatures(f), rownames(f.ave))
gene2cor.h<-intersect(VariableFeatures(h), rownames(h.ave))
gene2cor.fh<-intersect(VariableFeatures(f), rownames(h.ave))

cor.f<-cor(f.ave[gene2cor.f,], method="spearman")
cor.h<-cor(h.ave[gene2cor.h,], method="spearman")

cor.fh<-cor(f.ave[gene2cor.fh,],h.ave[gene2cor.fh,], method="spearman")

pheatmap::pheatmap(mat =cor.f, display_numbers = TRUE, color = hcl.colors(20,"Blue-Red3"), fontsize = 13, 
                   cluster_rows = F, cluster_cols = F)
pheatmap::pheatmap(mat =cor.h, display_numbers = TRUE, color = hcl.colors(20,"Blue-Red3"), fontsize = 13, 
                   cluster_rows = F, cluster_cols = F)
pheatmap::pheatmap(mat =cor.fh, display_numbers = TRUE, color = hcl.colors(20,"Blue-Red3"), fontsize = 13, 
                   cluster_rows = F, cluster_cols = F)

