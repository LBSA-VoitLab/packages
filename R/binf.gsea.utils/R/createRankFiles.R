#' createRankFiles
#'
#' This function allows to create a mapping for 2 sets of gene names based on NCBI gene_info file
#' @param gene_set_groups 
#' @keywords GSEA gmt rank
#' @export
#' @examples
#' createRankFileFrom2DE()
#' createTableAndSelect()
#' createRankFileForGSEA()
#' createFilteredRankFileForGSEA()
#' createDistFilteredRankFileForGSEA()

createRankFileFrom2DE <- function(res1,res2,mode=NULL,filename=NULL,p=NULL,fc=NULL){
  if(is.null(p)) p <- 0.01
  if(is.null(fc)) fc <- 1
  if(is.null(mode)) mode <- "intersect"
  
  intra.6.filter <- createTableAndSelect(res1,p,fc)
  res1.filter <- intra.6.filter$sel
  intra.7.filter <- createTableAndSelect(res2,p,fc)
  res2.filter <- intra.7.filter$sel
  
  rank.1.6 <- data.frame(round(sign(res1.filter$fc) * (-log(res1.filter$pVal)), digits=7))
  colnames(rank.1.6) <- "Rank"
  rank.6 <- data.frame(rank.1.6[!is.na(rank.1.6$Rank),])
  row.names(rank.6) <- row.names(res1.filter)[!is.na(rank.1.6$Rank)]
  colnames(rank.6) <- "Rank"
  
  rank.2.7 <- data.frame(round(sign(res2.filter$fc) * (-log(res2.filter$pVal)), digits=7))
  colnames(rank.2.7) <- "Rank"
  rank.7 <- data.frame(rank.2.7[!is.na(rank.2.7$Rank),])
  row.names(rank.7) <- row.names(res2.filter)[!is.na(rank.2.7$Rank)]
  colnames(rank.7) <- "Rank"
  
  if (mode == "intersect"){
    de <- merge(rank.6, rank.7, by=0, all=F)
    de$Rank <- rowMeans(cbind(abs(de$Rank.x),abs(de$Rank.y))) * round(sign(de$Rank.x)) * round(sign(de$Rank.y))
    de.f <-  data.frame(de$Rank)
    colnames(de.f) <- "Rank"
    row.names(de.f) <- de$Row.names
  }
  
  if (mode == "1"){
    se <- setdiff(rownames(rank.6), rownames(rank.7))
    de.f <- data.frame(rank.6[rownames(rank.6) %in% se,])
    colnames(de.f) <- "Rank"
    rownames(de.f) <- rownames(rank.6)[rownames(rank.6) %in% se]
  }
  
  if (mode == "2"){
    se <- setdiff(rownames(rank.7), rownames(rank.6))
    de.f <- data.frame(rank.7[rownames(rank.7) %in% se,])
    colnames(de.f) <- "Rank"
    rownames(de.f) <- rownames(rank.7)[rownames(rank.7) %in% se]
  }
  
  fn <- paste0("Z:/projects-bmed/projects/E6v7/GSEA/rank_files/",filename)
  write.table(de.f,fn,sep="\t",col.names = F,quote = F)
}

#from deseq results to table and selected table
createTableAndSelect <- function(res,p = NULL,fc = NULL){
  if(is.null(p)) p <- 0.01
  if(is.null(fc)) fc <- 1
  
  tab.1 <- data.frame(fc = res$log2FoldChange,pVal = res$pvalue,row.names = row.names(res))
  tab <- tab.1[!is.na(tab.1$pVal),]
  sel <- tab[tab$pVal < p & abs(tab$fc) > fc,]
  return(list(tab = tab,sel = sel))
}

#Create a rank file for GSEA
createRankFileForGSEA <- function(res,name){
  rank.1 <- data.frame(round(sign(res$log2FoldChange) * (-log(res$pvalue)), digits=7))
  colnames(rank.1) <- "Rank"
  rank <- data.frame(rank.1[!is.na(rank.1$Rank ),])
  row.names(rank) <- row.names(res)[!is.na(rank.1$Rank )]
  colnames(rank) <- "Rank"
  fn <- paste0("Z:/projects-bmed/projects/E6v7/GSEA/rank_files/",name)
  write.table(rank,fn,sep="\t",col.names = F,quote = F)
}

#Create a rank file for GSEA with filters
createFilteredRankFileForGSEA <- function(res,name,p=NULL,fc=NULL){
  if(is.null(p)) p <- 0.01
  if(is.null(fc)) fc <- 1
  
  inter.67.filter <- createTableAndSelect(res,p,fc)
  res.filter <- inter.67.filter$sel
  rank.1 <- data.frame(round(sign(res.filter$fc) * (-log(res.filter$pVal)), digits=7))
  colnames(rank.1) <- "Rank"
  rank <- data.frame(rank.1[!is.na(rank.1$Rank ),])
  row.names(rank) <- row.names(res.filter)[!is.na(rank.1$Rank )]
  colnames(rank) <- "Rank"
  fn <- paste0("Z:/projects-bmed/projects/E6v7/GSEA/rank_files/",name)
  write.table(rank,fn,sep="\t",col.names = F,quote = F)
}

#Create a rank file for GSEA with filters and distance
createDistFilteredRankFileForGSEA <- function(res,name,p=NULL,fc=NULL,dist=NULL, d=NULL){
  if(is.null(p)) p <- 0.01
  if(is.null(fc)) fc <- 1
  if(is.null(dist)){
    createFilteredRankFileForGSEA(res,name,p,fc)
    stop("distance file not specified")
  }
  if(is.null(d)) d <- 0.98
  
  inter.67.filter <- createTableAndSelect(res,p,fc)
  res.filter <- inter.67.filter$sel
  ind <- dist$dist < d
  geneList <- dist$gene[ind]
  ind2 <- rownames(res.filter) %in% geneList
  res.filter.2 <- res.filter[ind2,]
  
  rank.1 <- data.frame(round(sign(res.filter.2$fc) * (-log(res.filter.2$pVal)), digits=7))
  colnames(rank.1) <- "Rank"
  rank <- data.frame(rank.1[!is.na(rank.1$Rank ),])
  row.names(rank) <- row.names(res.filter.2)[!is.na(rank.1$Rank )]
  colnames(rank) <- "Rank"
  
  fn <- paste0("Z:/projects-bmed/projects/E6v7/GSEA/rank_files/",name)
  write.table(rank,fn,sep="\t",col.names = F,quote = F)
}

#Create a rank file for GSEA with exagerations : more significant => higher 
createExageratedRankFileForGSEA <- function(res,name,p=NULL,fc=NULL){
  if(is.null(p)) p <- 0.01
  if(is.null(fc)) fc <- 1
  
  inter.67.filter <- createTableAndSelect(res,p,fc)
  res.filter <- inter.67.filter$tab
  
  res.filter$exageration = apply(res.filter, 1, function(x) {if (x["fc"] > 1 & x["pVal"] < 0.01){16}
    else if (x["fc"] > 0.5 & x["pVal"] < 0.01){8}
    else if (x["fc"] > 0.5 & x["pVal"] < 0.05){4}
    else if (x["fc"] > 0.5 & x["pVal"] < 0.1){2}
    else {1}})
  
  rank.1 <- data.frame(round(sign(res.filter$fc) * (-log(res.filter$pVal)) * res.filter$exageration, digits=7))
  colnames(rank.1) <- "Rank"
  rank <- data.frame(rank.1[!is.na(rank.1$Rank ),])
  row.names(rank) <- row.names(res.filter)[!is.na(rank.1$Rank )]
  colnames(rank) <- "Rank"
  
  
  fn <- paste0("Z:/projects-bmed/projects/E6v7/GSEA/rank_files/",name)
  write.table(rank,fn,sep="\t",col.names = F,quote = F)
}
