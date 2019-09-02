#' gene_symbol2enterez Function
#'
#' This function allows you to create bar plots directlty from gsea results
#' @param gsea_res_folder folder location of GSEA result
#' @keywords barplot gsea
#' @export
#' @examples
#' gene_symbol2enterez()

library("reader")
#source('https://bioconductor.org/biocLite.R')
#biocLite('org.Hs.eg.db')
#biocLite('BSgenome.Mfascicularis.NCBI.5.0')
library('org.Hs.eg.db')
library('org.Mmu.eg.db')

ls1 <- gene_symbol2enterez(rownames(E06_final_df))
ls2 <- gene_symbol2enterez(rownames(E07_final_df))

gene_symbol2enterez <- function(geneList){
  a <-mapIds(org.Hs.eg.db, geneList, 'ENTREZID', 'SYMBOL')
  mapped_genes <- mappedkeys(a)
  id <- as.list(a[mapped_genes])
#  return(id)
  return(as.list(a))
}
