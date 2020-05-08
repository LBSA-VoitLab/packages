#' plot_multi_GSEA_heatmap Function
#'
#' This function allows you to heatmaps for genes specific to genesets
#' @param gsea_res_folder folder location of GSEA result
#' @keywords barplot gsea
#' @export
#' @examples
#' plot_multi_GSEA_heatmap()

plot_multi_GSEA_heatmap <- function(gsea_pathnames,name_list,fdr = NULL,pVal = NULL){
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

  return(plot_heatmap(df))
}

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

#only plot for modules mentioned
plot_multi_GSEA_heatmap_modular.2 <- function(gsea_pathnames,name_list,module_list,fdr = NULL,pVal = NULL){
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
  
  df.1 <- df[module_list,]
  return(plot_heatmap(df.1))
}

#only plot for modules that either clears the criterion
plot_multi_GSEA_heatmap_modular.3 <- function(gsea_pathnames,name_list,fdr = NULL,pVal = NULL, nes = NULL){
  if(is.null(fdr)) fdr <- 0.25
  if(is.null(pVal)) pVal <- 0.05
  if(is.null(pVal)) nes <- 0
  
  df.2 = data.table(stringsAsFactors = FALSE)
  names = c()
  for (gsea_results in gsea_pathnames){
    posList <- list.files(path=gsea_results,pattern="gsea_report_for_.*_pos_.*.xls")
    negList <- list.files(path=gsea_results,pattern="gsea_report_for_.*_neg_.*.xls")
    posData <- read.csv2(paste0(gsea_results,"/",posList), header=TRUE,sep="\t")
    negData <- read.csv2(paste0(gsea_results,"/",negList), header=TRUE,sep="\t")
    full <- rbind(posData, negData)

    dt <- data.table(Name = full$NAME,NES = full$NES, pVal = full$NOM.p.val, fdr = full$FDR.q.val)
    if (match(gsea_results,gsea_pathnames) == 1){
      df.2 <- dt
    }
    else {
      df.2 <- merge(df.2,dt,by = "Name",all = T)
    }
  }
  
  df <- df.2[(abs(as.numeric(as.character(df.2$NES.x))) > nes | 
               abs(as.numeric(as.character(df.2$NES.y))) > nes) &
               (as.numeric(as.character(df.2$fdr.x)) < fdr |
               as.numeric(as.character(df.2$fdr.y)) < fdr) &
               (as.numeric(as.character(df.2$pVal.x)) < pVal |
               as.numeric(as.character(df.2$pVal.y)) < pVal),c(1,2,5)]
  
  df[is.na(df)] <- 0
  df <- df[order(df$Name),]
  df.1 <- as.matrix(df[,-1])
  row.names(df.1) <- df$Name
  colnames(df.1) <- name_list
#  p <- plot_heatmap(df.1)
  return(df.1)
}


plot_heatmap <- function(df){
  require(reshape2)
  require(ggplot2)
  
  melted_df <- melt(df)
  melted_df$value=as.numeric(levels(melted_df$value))[melted_df$value]
  
  p <- ggplot(data = melted_df, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile(color = "white") +
    scale_fill_gradient2(low="red1", high="dodgerblue2", mid = gray(0.9), 
                         midpoint = 0, limit = c(-2.5,2.5), space = "Lab", 
                         name="NES")
  return(p)
}
