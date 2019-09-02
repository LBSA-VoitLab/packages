#' plot_multi_GSEA_heatmap_modular_group Function
#'
#' This function allows you to heatmaps for genes specific to genesets
#' @param folder location of GSEA result
#' @keywords barplot gsea
#' @export
#' @examples
#' @docType data
#' @usage data(modular_data)
#' getAllGenes()

data(modular_data)

modular_data <- readRDS("Z:/projects-bmed/projects/packages/R/binf.modular/data/modular_data.RData")
#Function to get all genes from modular analysis
getAllGenes <- function(){

}


#gsea heatmap for modular genesets
plot_multi_GSEA_heatmap_modular <- function(gsea_pathnames,name_list,fdr = NULL,pVal = NULL){
  if(is.null(fdr)) fdr <- 0.25
  if(is.null(pVal)) pVal <- 0.05

  df.2 = data.table(stringsAsFactors = FALSE)
  names = c()
  for (gsea_results in gsea_pathnames){
    posList <- list.files(path=gsea_results,pattern="gsea_report_for_.*_pos_.*.xls")
    negList <- list.files(path=gsea_results,pattern="gsea_report_for_.*_neg_.*.xls")
    posData <- read.csv2(paste0(gsea_results,"/",posList), header=TRUE,sep="\t")
    negData <- read.csv2(paste0(gsea_results,"/",negList), header=TRUE,sep="\t")
    full <- rbind(posData, negData)
    names <- union(names,as.character(full[(as.numeric(as.character(full$NOM.p.val)) < pVal) & (as.numeric(as.character(full$FDR.q.val)) < fdr),"NAME"]))

    dt <- data.table(Name = full$NAME,NES = full$NES)
    if (match(gsea_results,gsea_pathnames) == 1){
      df.2 <- dt
    }
    else {
      df.2 <- merge(df.2,dt,by = "Name",all = T)
    }
  }
  df <- as.matrix(df.2[df.2$Name %in% names,-1])
  df[is.na(df)] <- 0
  row.names(df) <- df.2$Name[df.2$Name %in% names]
  colnames(df) <- name_list

  df.1 <- df[!grepl(".*:TBD",rownames(df)),]
  return(plot_heatmap(df.1))
}

#gsea heatmap for modular with full list of modules grouped by functional
plot_multi_GSEA_heatmap_modular_group <- function(gsea_pathnames,name_list){

  df.2 = data.table(stringsAsFactors = FALSE)
  names = c()
  for (gsea_results in gsea_pathnames){
    posList <- list.files(path=gsea_results,pattern="gsea_report_for_.*_pos_.*.xls")
    negList <- list.files(path=gsea_results,pattern="gsea_report_for_.*_neg_.*.xls")
    posData <- read.csv2(paste0(gsea_results,"/",posList), header=TRUE,sep="\t")
    negData <- read.csv2(paste0(gsea_results,"/",negList), header=TRUE,sep="\t")
    full <- rbind(posData, negData)
    names <- union(names,as.character(full[,"NAME"]))

    dt <- data.table(Name = full$NAME,NES = full$NES)
    if (match(gsea_results,gsea_pathnames) == 1){
      df.2 <- dt
    }
    else {
      df.2 <- merge(df.2,dt,by = "Name",all = T)
    }
  }
  gr <- modular_data[,c(1,3)]
  colnames(gr)[2] <- "Group"
  gr$Name <- toupper(paste(gr$Module.ID,gr$Group,sep=":"))
  df.3 <- merge(df.2,gr,by = "Name")
  df.4 <- df.3[df.3$Group != 'TBD',]
  df.4$Group = as.factor(df.4$Group)
  df.5 <- df.4[order(df.4$Group),]

  cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
  warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
  cols = c(rev(cool), rev(warm))
  data_col <- colorRampPalette(cols)(255)

  gr_cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(nlevels(df.5$Group))
  gr_cols <- colorRamps::primary.colors(30)

  df <- as.matrix(df.5[df.5$Name %in% names,c(2,3)])
  class(df) <- "numeric"
  row.names(df) <- df.5$Name
  colnames(df) <- name_list
  df[is.na(df)] <- 0

  #  dev.off()
  #  p <- heatmap.2(df,mar=c(10,20),trace="none",
  #            col = data_col,RowSideColors = gr_cols[df.5$Module.functional.association.title],
  #            Rowv = F,Colv = F)

  return(heatmap.2(df,mar=c(10,20),trace="none",dendrogram = 'none',
                   col = data_col,RowSideColors = gr_cols[df.5$Group],
                   Rowv = F,Colv = F))
}

#gsea heatmap for modular with full list of modules clustered by rows
plot_multi_GSEA_heatmap_modular_cluster <- function(gsea_pathnames,name_list){

  df.2 = data.table(stringsAsFactors = FALSE)
  names = c()
  for (gsea_results in gsea_pathnames){
    posList <- list.files(path=gsea_results,pattern="gsea_report_for_.*_pos_.*.xls")
    negList <- list.files(path=gsea_results,pattern="gsea_report_for_.*_neg_.*.xls")
    posData <- read.csv2(paste0(gsea_results,"/",posList), header=TRUE,sep="\t")
    negData <- read.csv2(paste0(gsea_results,"/",negList), header=TRUE,sep="\t")
    full <- rbind(posData, negData)
    names <- union(names,as.character(full[,"NAME"]))

    dt <- data.table(Name = full$NAME,NES = full$NES)
    if (match(gsea_results,gsea_pathnames) == 1){
      df.2 <- dt
    }
    else {
      df.2 <- merge(df.2,dt,by = "Name",all = T)
    }
  }
  gr <- modular_data[,c(1,3)]
  colnames(gr)[2] <- "Group"
  gr$Name <- toupper(paste(gr$Module.ID,gr$Group,sep=":"))
  df.3 <- merge(df.2,gr,by = "Name")
  df.4 <- df.3[df.3$Group != 'TBD',]
  df.4$Group = as.factor(df.4$Group)
  df.5 <- df.4[order(df.4$Group),]

  #  cool = YlGnBu(50)
  #  warm = YlOrRd(50)
  #  tot = PRGn(100)
  cool = rainbow(50, start=rgb2hsv(col2rgb('green'))[1], end=rgb2hsv(col2rgb('blue'))[1])
  warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('green'))[1])
  cols = c(rev(cool), rev(warm))
  #  cols = tot
  data_col <- colorRampPalette(cols)(255)

  #  gr_cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(nlevels(df.5$Group))
  gr_cols <- colorRamps::primary.colors(30)

  df <- as.matrix(df.5[df.5$Name %in% names,c(2,3)])
  class(df) <- "numeric"
  row.names(df) <- df.5$Name
  colnames(df) <- name_list
  df[is.na(df)] <- 0

  #  dev.off()
  #  p <- heatmap.2(df,mar=c(10,20),trace="none",
  #                 col = data_col,RowSideColors = gr_cols[df.5$Group],
  #                 Colv = F)

  return(heatmap.2(df,mar=c(10,20),trace="none",
                   col = data_col,RowSideColors = gr_cols[df.5$Group],
                   Colv = F))
}

plot_heatmap <- function(df){
  require(reshape2)
  require(ggplot2)

  melted_df <- melt(df)
  melted_df$value=as.numeric(levels(melted_df$value))[melted_df$value]

  p <- ggplot(data = melted_df, aes(x=Var2, y=Var1, fill=value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-2.5,2.5), space = "Lab",
                         name="NES")
  return(p)
}
