
#' getOutliers
#'
#' This function allows to perform deseq2 diff exp analysis 
#' @param gene_set_groups 
#' @keywords GSEA gmt
#' @export
#' @examples
#' getOutliers()



#from deseq results table (or sub table) to find outlier DEGs
getOutliers <- function(tab,value = "fc",mode = "neg"){
  if(is.null(mode)) mode <- "both"
  if(is.null(value)) value <- "fc"

  data_mean = mean(tab[["fc"]])
  data_std =  sd(tab[["fc"]])
  # identify outliers
  cut_off = data_std * 2
  lower = data_mean - cut_off 
  upper = data_mean + cut_off
  
  # identify outliers
  if (mode == "pos") outliers = tab[tab[["fc"]] > upper,]
  if (mode == "neg") outliers = tab[tab[["fc"]] < lower,]
  else outliers = tab[tab[["fc"]] < lower | tab[["fc"]] > upper,]

  return(outliers)
}
