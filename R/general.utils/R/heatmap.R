#' plot_heatmap Function
#'
#' This function allows you to heatmaps for genes specific to genesets
#' @param df dataframe with values corresponding to the color
#' @keywords heatmap
#' @export
#' @examples
#' plot_heatmap()


plot_heatmap <- function(df){
  df$rownames <- rownames(df)
  df.m <- melt(df)
#  df.m <- ddply(df.m, .(variable), transform,rescale = rescale(value))
  p <- ggplot(df.m, aes(variable, rownames)) + 
    geom_tile(aes(fill = value),colour = "white") + 
    scale_fill_gradient2(low = "red", mid = "white",high = "blue", midpoint = 0,na.value = "grey50")
  return(p)
}

#plot heatmap cluster rows
plot_heatmap.2 <- function(df){
  df.dendro <- as.dendrogram(hclust(d = dist(x = df)))
  df.order <- order.dendrogram(df.dendro)
  
  df$rownames <- rownames(df)
  df.m <- melt(df)
  
  # Order the levels according to dendogram order
  df.m$rownames <- factor(x = df.m$rownames,
                                 levels = rownames(df)[df.order], 
                                 ordered = TRUE)
  
  #  df.m <- ddply(df.m, .(variable), transform,rescale = rescale(value))
  p <- ggplot(df.m, aes(variable, rownames)) + 
    geom_tile(aes(fill = value),colour = "white") + 
    scale_fill_gradient2(low = "red", mid = "white",high = "blue", midpoint = 0,na.value = "grey50")
  return(p)
}