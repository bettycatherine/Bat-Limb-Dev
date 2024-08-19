library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(tidyverse)
library(hrbrthemes)
library(circlize)
library(kableExtra)
options(knitr.table.format = "html")
library(viridis)
library(igraph)
library(ggraph)
library(colormap)

##all loom files were processed by pyscenic

library(SCopeLoomR)
loom <- open_loom("rsi_3-split-250-merged-scenic_integrated-output.loom", mode = "r")

# Read information from loom file:
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(loom)
values<-names(regulonsAucThresholds)
values<-as.numeric(values)
regulonsAucThresholds<-unname(regulonsAucThresholds)
#density plot od threshold, and used as threshold standard for comparing forelimb and hindlimb TF activity 
ggplot(as.data.frame(values), aes(x=values)) + geom_density()
names(values)<-regulonsAucThresholds
embeddings <- get_embeddings(loom)
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
clusterings <- get_clusterings_with_name(loom)
cellInfo <- data.frame(cellInfo)
##change
celltype<-read.csv("metadata-all.csv",row.names = 1)
cellInfo$celltype<-celltype
cellTypeColumn <- "Class"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "celltype"
close_loom(loom)
motifEnrichmentFile <- file.path("reg-all.csv")
motifEnrichment <- data.table::fread(motifEnrichmentFile, header=T, skip=1)[-3,]
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")

regulonsAUC
#Average Regulon Activity by cluster
regulonsAUC <- regulonsAUC[onlyNonDuplicatedExtended(rownames(regulonsAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                                     function(cells) rowMeans(getAUC(regulonsAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA, show_rownames = F)

auc_all <- regulonActivity_byCellType_Scaled

#####subset_celltypes
regulonActivity_byCellType_subset <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                                     function(cells) getAUC(regulonsAUC)[,cells])
for (i in 1:length(regulonActivity_byCellType_subset)) {
  regulonActivity_byCellType_subset[[i]] <- t(scale(t(regulonActivity_byCellType_subset[[i]]), center = T, scale=T))
}
pheatmap::pheatmap(regulonActivity_byCellType_subset$F_chondrocyte,
                   show_rownames = F,
                   show_colnames = F)

####subset each cell population, chondrocyte as an exmple
auc_ch<-auc_all[,2:3]
set.seed(123)
obj<-pheatmap::pheatmap(auc_ch, #fontsize_row=3, 
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA, kmeans_k = 10
                   )

data <- paste0(auc_ch, obj$kmeans$cluster)

data_long <- data %>%
  gather(key = 'key', value = 'value', -cluster)
data_long %>% 
  ggplot(aes(factor(cluster),value, fill=key)) + 
  geom_boxplot()+
  labs(title="auc of chondrocyte")

write.csv(data, "scenic-auc-chondrocyte.csv")


