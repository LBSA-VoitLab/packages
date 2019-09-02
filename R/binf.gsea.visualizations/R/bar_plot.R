#' gsea_bar_plot Function
#'
#' This function allows you to create bar plots directlty from gsea results
#' @param gsea_res_folder folder location of GSEA result
#' @keywords barplot gsea
#' @export
#' @examples
#' gsea_bar_plot()

library("reader")
library("readxl")
library("xlsx")
library("ggplot2")
library("dplyr")

gsea_bar_plot <- function(gsea_res_folder,posCount,negCount,FDRq,pV){

  s <- getSignificantGenesetsTable(gsea_res_folder,posCount,negCount,FDRq,pV)

  p <- ggplot2::ggplot(s, ggplot2::aes(x=reorder(NAME,NES),y=NES)) +
    ggplot2::geom_bar(stat = "identity",ggplot2::aes(fill = NOM.p.val)) +
    ggplot2::scale_y_continuous() +
    ggplot2::coord_flip()
  return(p)
}

getSignificantGenesetsTable <- function(gsea_res_folder,posCount,negCount,FDRq,pV){
  if(missing(posCount)) {
    posCount <- 20
  } else {
    posCount <- posCount
  }
  if(missing(negCount)) {
    negCount <- 20
  } else {
    negCount <- negCount
  }
  if(missing(FDRq)) {
    FDRq <- 0.25
  } else {
    FDRq <- FDRq
  }
  if(missing(pV)) {
    pV <- 0.01
  } else {
    pV <- pV
  }

  posList <- list.files(path=gsea_res_folder,pattern="gsea_report_for_.*_pos_.*.xls")
  negList <- list.files(path=gsea_res_folder,pattern="gsea_report_for_.*_neg_.*.xls")
  posData <- read.csv2(paste0(gsea_res_folder,"/",posList), header=TRUE,sep="\t")
  negData <- read.csv2(paste0(gsea_res_folder,"/",negList), header=TRUE,sep="\t")
  subPos <- posData[(as.numeric(as.character(posData$FDR.q.val)) < FDRq) & (as.numeric(as.character(posData$NOM.p.val)) < pV),]
  subNeg <- negData[(as.numeric(as.character(negData$FDR.q.val)) < FDRq) & (as.numeric(as.character(negData$NOM.p.val)) < pV),]
  sub <- rbind(dplyr::top_n(subPos,posCount,wt = NES),dplyr::top_n(subNeg,negCount,wt = NES))
  sub$NES=as.numeric(levels(sub$NES))[sub$NES]
  sub$NOM.p.val=as.numeric(levels(sub$NOM.p.val))[sub$NOM.p.val]

  return(sub[rev(order(as.numeric(as.character(sub$NES)))),])
}
