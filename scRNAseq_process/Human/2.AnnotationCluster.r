library(dplyr)
library(Seurat)
library(stringr)
library(cowplot)
library(viridis)
library(ggpubr)
library(SingleR)
library(scater)
library(ExperimentHub)
library(AnnotationHub)
library(AnnotationDbi)
library(celldex)
library(ggsci)
library(stringr)
library(RColorBrewer)
library(scRepertoire)
library(scales)
dtVar <- Sys.Date() 
dtVar <- as.Date(dtVar, tz="UTC")

source('aPD1_Ery_scripts/scRNAseq_process/Human/stat.r')
workpath <- "aPD1_Ery_scripts/scRNAseq_process/Human/"
setwd(workpath)

rdsfile <- str_c("./P1_3_t-SNE_30PCA_0.6Resolution/P1_3_t-SNE_30PCA_0.6Resolution.rds")
immune.combined <- readRDS(rdsfile)
immune.combined$seurat_clusters <- Idents(immune.combined)
p1 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE, label.size = 5) + NoLegend() + theme_bw()
p1 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE, label.size = 5, split.by='sample_label') + NoLegend() + theme_bw()
ggsave(str_c('WholeCells.splitbysample_label.',dtVar,'.pdf'), p1, width=11, height=5)
DefaultAssay(immune.combined) <- "RNA"
major_cluster <- c('CD3E','CD3G','CD4','CD8A','CD8B','TRDC','NKG7','KLRB1','CD79A','MS4A1','MKI67','TYMS','FCGR3A','CD14','CD68','S100A8','S100A9','LYZ','CD163','IGKC','LILRA4','TCF4','IRF4','MZB1','XBP1','CD1C','FCER1A','CST3','CD86','PPBP','HBB','HBA2','HBA1','NLRP3','PLTP','IL1A','IL1B','CD34')
# major_cluster <- c('CD3E','CD3G','CD4','CD8A','CD8B','TRDC','NKG7','KLRB1','FCGR3A','CD79A','MS4A1','MKI67','TYMS','CD68','S100A8','S100A9','CD14','LYZ','CD163','ITGAX','CD1C','FCER1A','LILRA4','JCHAIN','IGKC','PPBP','HBB','HBA2','HBA1')
Idents(immune.combined) <- factor(immune.combined$seurat_clusters, levels=c(3,4,2,6,9,18,16,13,7,5,14,17,10,0,1,20,11,8,19,21,22,12,15))
# p1 <- DotPlot(immune.combined, features = major_cluster, cols = c("blue", "red", "grey"), dot.scale = 2) + RotatedAxis() 
p1 <- DotPlot(immune.combined, features = major_cluster, col.min=-1, cols=c('white','red','red4'), dot.scale = 3) + RotatedAxis() +theme_bw()+
        theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5)) 
ggsave(str_c('WholeCells.marker_dotplot.',dtVar,'.pdf'), p1, width=9.5, height=5)

immune.combined$celltype = dplyr::case_when(
      immune.combined$seurat_clusters %in% c(3,4) ~ "CD4 T Cell",
      immune.combined$seurat_clusters %in% c(2,6,9,18,16) ~ "CD8 T Cell",
      immune.combined$seurat_clusters %in% c(13) ~ "NKT Cell",
      immune.combined$seurat_clusters %in% c(7) ~ "NK Cell",
      immune.combined$seurat_clusters %in% c(5,14) ~ "B Cell",
      immune.combined$seurat_clusters %in% c(17) ~ "Proliferation Cell",
      immune.combined$seurat_clusters %in% c(10,0,20,11,8,19) ~ "Monocyte", # CD68+ FCGR3A+ 
      immune.combined$seurat_clusters %in% c(21,22) ~ "pDC",
      immune.combined$seurat_clusters %in% c(12) ~ "cDC",
      immune.combined$seurat_clusters %in% c(15) ~ "Platelet",
      immune.combined$seurat_clusters %in% c(1) ~ "Unknown")
total_celltype_levels = c('CD4 T Cell','CD8 T Cell','NKT Cell','NK Cell','B Cell','Proliferation Cell','Monocyte','Platelet','pDC','cDC',"Unknown")
total_celltype_colors = c(pal_nejm(alpha = 0.2)(8), pal_npg("nrc", alpha = 0.2)(10))
immune.combined$celltype <- factor(immune.combined$celltype, levels=total_celltype_levels)
Idents(immune.combined) <- immune.combined$celltype 
p1 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE, repel=TRUE, label.size = 2.5, cols=total_celltype_colors) + NoLegend() + theme_bw()
ggsave(str_c('WholeCells',dtVar,'.pdf'), p1, width=6, height=3.5)
p1 <- DimPlot(immune.combined, reduction = "tsne", label = FALSE, label.size = 2.5, repel=FALSE, cols=total_celltype_colors)+theme_bw() +
      theme(legend.position="bottom", axis.text=element_blank())+ NoLegend()+ylab('')+xlab('')
ggsave(str_c('WholeCells.',dtVar,'.no.pdf'), p1, width=3, height=3)
p1 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE, repel=TRUE, label.size = 2.5, cols=c(pal_nejm(alpha = 0.6)(8)[3:4], pal_npg("nrc", alpha = 0.6)(10)[3:10])) + NoLegend() + theme_bw()
ggsave('WholeCells_legend.pdf', p1, width=5, height=3.5)
Idents(immune.combined) <- factor(immune.combined$celltype, levels=c('CD4 T Cell','CD8 T Cell','NKT Cell','NK Cell','B Cell','Proliferation Cell','Monocyte','pDC','cDC','Platelet',"Unknown"))
p1 <- DotPlot(immune.combined, features = major_cluster, col.min=-1, cols=c('white','red','red4'), dot.scale = 3) + RotatedAxis() +theme_bw()+
        theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5)) 
ggsave(str_c('WholeCells_withanno.marker_dotplot.',dtVar,'.pdf'), p1, width=9, height=4.2)
### cellcycle
DefaultAssay(immune.combined) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
immune.combined <- CellCycleScoring(immune.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
p1 <- DimPlot(immune.combined, reduction = "tsne",cols=brewer.pal(9,"Paired"), pt.size=0.01)+theme_bw()#, split.by="Phase")
ggsave('WholeCells_cellcycle.pdf', p1, width=4, height=3.5)
### 保存结果
Idents(immune.combined) <- 'celltype'
saveRDS(immune.combined, file = str_c("./P1_3_t-SNE_30PCA_0.6Resolution.AnnoManual.rds"))



immune.combined <- readRDS(str_c("./P1_3_t-SNE_30PCA_0.6Resolution/P1_3_t-SNE_30PCA_0.6Resolution.AnnoManual.rds"))
Idents(immune.combined) <- 'celltype'
cellnumber_summary(immune.combined, names(table(immune.combined$sample_label)), total_celltype_levels, total_celltype_colors, 'immune.combined')