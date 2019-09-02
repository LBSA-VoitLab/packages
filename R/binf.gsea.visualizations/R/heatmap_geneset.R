#' gsea_bar_plot Function
#'
#' This function allows you to heatmaps for genes specific to genesets
#' @param gsea_res_folder folder location of GSEA result
#' @keywords barplot gsea
#' @export
#' @examples
#' gsea_bar_plot()

library("reader")
library("readxl")
library("xlsx")
library("ggplot2")
library("gplots")

plot_gene_set_heatmap <- function(geneSets,gene_set_groups,expression_data){
  genes <- c()
  ls <- data.frame(matrix(ncol = dim(expression_data)[2], nrow = 0))
  colnames(ls) <- colnames(expression_data)
  for (i in c(1:length(geneSets))){
    gns <- intersect(rownames(expression_data),gsg[[geneSets[i]]])
    addGns <- setdiff(gns,genes)
    genes <- append(genes,addGns)
    thisData <- expression_data[addGns,]
    ls <- rbind(ls,thisData)
  }
  col.order <- hclust(dist(t(ls)))$order
  dat_new <- ls[, col.order]
  dat_new$genes <- rownames(dat_new)
  
  df_molten <- melt(dat_new)
  df_molten_add <- as.data.frame(tstrsplit(df_molten$variable, "_",fixed =T))
  df_molten <- cbind(df_molten,df_molten_add)
  colnames(df_molten) <- c("genes","id","count","Exp","type","sample","TP")
  df_molten$genes <- factor(df_molten$genes, levels = genes)
  
  p<- ggplot(data = df_molten,
             aes(x = id, y = genes, fill = count)) + 
    geom_raster() +
    xlab("Genes") + 
    scale_fill_distiller(palette = "RdYlBu", trans = "log10") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(angle = 0, hjust = 1)) + 
    ggtitle("ggplot heatmap")
  
  return(p)
}

plot_gene_set_heatmap.2 <- function(geneSets,gene_set_groups,expression_data){
  library(RColorBrewer) 
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  
  geneGrCol <- colorRampPalette(brewer.pal(9, "Spectral"))(length(geneSets))
  
  
  genes <- c()
  geneGrCollist <- c()
  
  ls <- data.frame(matrix(ncol = dim(expression_data)[2], nrow = 0))
  colnames(ls) <- colnames(expression_data)
  for (i in c(1:length(geneSets))){
    gns <- intersect(rownames(expression_data),gsg[[geneSets[i]]])
    addGns <- setdiff(gns,genes)
    addGnsCol <- rep(geneGrCol[i],length(addGns))
    genes <- append(genes,addGns)
    geneGrCollist <- append(geneGrCollist,addGnsCol)
    thisData <- expression_data[addGns,]
    ls <- rbind(ls,thisData)
  }
  dat2 <- ls
  dat2_samples <- as.data.frame(tstrsplit(colnames(dat2), "_",fixed =T))
  colnames(dat2_samples) <- c("Exp","type","sample","TP")
  dat2_samples$Status <- rep("Post",dim(dat2_samples)[1])
  dat2_samples$Status[dat2_samples$TP %in% c("T01","T02","T2B")] <- rep("Pre",length(dat2_samples$Status[dat2_samples$TP %in% c("T01","T02","T2B")]))
  
  expGrCol <- colorRampPalette(brewer.pal(9, "Spectral"))(length(levels(factor(paste(dat2_samples$Status,dat2_samples$TP,sep = "_")))))
  
  dev.off() 
  p <- heatmap.2(as.matrix(dat2),
                 trace="none",
                 scale="row",
                 RowSideColors = geneGrCollist,
                 ColSideColors = expGrCol[factor(paste(dat2_samples$Exp,dat2_samples$Status,sep = "_"))],
                 col="terrain.colors",
                 dendrogram ="column",
                 margins = c(10, 5),
                 Rowv = F) #col=hmcol

  
  return(p)
}

