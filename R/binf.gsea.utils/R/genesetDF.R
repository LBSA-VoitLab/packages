#' get_geneset_table
#'
#' This function allows to create a list all genesets with enrichment
#' @param gsea_res_folderget 
#' @keywords GSEA folder
#' @export
#' @examples
#' get_gsea_df()


get_gsea_df <- function(gsea_res_folder){
  posList <- list.files(path=gsea_res_folder,pattern="gsea_report_for_.*_pos_.*.xls")
  negList <- list.files(path=gsea_res_folder,pattern="gsea_report_for_.*_neg_.*.xls")
  posData <- read.csv2(paste0(gsea_res_folder,"/",posList), header=TRUE,sep="\t")
  negData <- read.csv2(paste0(gsea_res_folder,"/",negList), header=TRUE,sep="\t")
  df <- rbind(posData,negData)
  df$NES <- as.numeric(as.character(df$NES))
  df <- df[!is.na(df$NES),]
  df <- df[order(-df$NES),]
  return(df[,c("NAME","NES")])
}


