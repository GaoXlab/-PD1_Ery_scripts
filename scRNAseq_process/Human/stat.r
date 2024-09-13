library(dplyr)
library(Seurat)
library(stringr)
library(cowplot)
library(viridis)
library(ggpubr)
library(scater)
library(metap)
library(paletteer)
library(ggsci)
library(Hmisc)
library(reshape2)
library(plyr)
library(rlist)
library(tidyverse)
library(ggplot2)

data_summary <- function(data, varname, groupnames){
      require(plyr)
      summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
          sd = sd(x[[col]], na.rm=TRUE))
      }
      data_sum<-ddply(data, groupnames, .fun=summary_func,
                      varname)
      data_sum <- rename(data_sum, c("mean" = varname))
      return(data_sum)
    }

cellnumber_summary <- function(rds, sample_labels, celltypes, colors, name){
    sample_groups = sample_labels
    library(reshape2)
    Idents(rds) <- factor(rds$celltype, levels=celltypes)
    cell_num_inclusters <- c()
    cellnumber <- c()
    cellratio_in_sample <- c()
    sample_label_list <- names(table(rds$sample_label))
    for(clusters_i in celltypes){
    	cellnum_in_cluster <- WhichCells(rds, idents = clusters_i)
    	sample_label_cell_num <- c()
    	cellratio_in_sample_label <- c()
    	for(sample_label_i in sample_label_list){
            sample_label_cell_num <- c(sample_label_cell_num, length(grep(str_c('^', sample_label_i, '_'), cellnum_in_cluster)))
    		cellratio_in_sample_label <- c(cellratio_in_sample_label, length(grep(str_c('^', sample_label_i, '_'), cellnum_in_cluster))/table(rds$sample_label)[sample_label_i])
      }
    	cellnumber <- rbind(cellnumber, sample_label_cell_num)
    	cellratio_in_sample <- rbind(cellratio_in_sample, cellratio_in_sample_label)
    	cell_num_inclusters <- c(cell_num_inclusters, length(cellnum_in_cluster))
    }
    cells_statistics <- as.data.frame(cbind(cellnumber, cellnumber/cell_num_inclusters, cellratio_in_sample))
    colnames(cells_statistics) <- c(str_c('CellNum',sample_label_list), str_c('CellRatio_in_Clusters',sample_label_list), str_c('CellRatio_in_Samples',sample_label_list))
    rownames(cells_statistics) <- celltypes
    write.table(cells_statistics, str_c(name, 'cells_statistics.log'), sep='\t')
    # cells_statistics <- read.table('BoneMarrow_cellsratio.csv', sep=',')
    df3 <- melt(cbind(rownames(cells_statistics), cells_statistics[grep('CellRatio_in_Samples', colnames(cells_statistics))]))
    df3$label <- gsub("CellRatio_in_Samples", "", as.character(df3$variable))
    colnames(df3) <- c('celltype', 'variable', 'value', 'label')
    df3$label <- factor(df3$label, levels=sample_labels)
    df3$celltype <- factor(df3$celltype, levels=celltypes)
    df3$value100 <- df3$value*100
    df3$group <- gsub('_1$|_2$|_3$|_4$','', df3$label)
    df3$group <- factor(df3$group, sample_groups)
    p1 <- ggplot(data=df3, aes(x=label, y=value100, fill=celltype)) + geom_bar(stat="identity")+theme_classic()+facet_grid(.~celltype)+
          theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.5))+ylab('Percentage (%)')+
          scale_fill_manual(values=colors) + theme(legend.position='bottom')
    p2 <- ggplot(data=df3, aes(x=label, y=value100, fill=celltype)) + geom_bar(stat="identity")+theme_classic()+ylab('Percentage (%)')+
          scale_fill_manual(values=colors) + theme(legend.position='right')
    p3 <- ggplot(data=df3, aes(x=group, y=value100, fill=celltype)) + geom_bar(stat="identity")+theme_classic()+ylab('Percentage (%)')+
          scale_fill_manual(values=colors) + theme(legend.position='right')
    g <- plot_grid(p1, p2, p3, rel_widths=c(1,1,0.4), ncol=3)
    ggsave(str_c("./", name, "cellsratio_InAnnotationClusters.barplot2.pdf"), g, width=16, height=5)
}
