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
    scale_fill_gradient2(low = "red1", mid = "grey60",high = "dodgerblue2", midpoint = 0,na.value = "grey10")
  return(p)
}

#plot heatmap of DEGs with circles
plot_heatmap_genes_shape <- function(df){
#  df.dendro <- as.dendrogram(hclust(d = dist(x = df)))
#  df.order <- order.dendrogram(df.dendro)
  
  df$rownames <- rownames(df)
  df.m <- melt(df)
  
  # Order the levels according to dendogram order
  df.m$rownames <- fct_rev(factor(x = df.m$rownames,
                          levels = rownames(df), 
                          ordered = TRUE))

  
  #  df.m <- ddply(df.m, .(variable), transform,rescale = rescale(value))
  p <- ggplot(df.m, aes(variable, rownames)) + 
#    geom_text(size = 5) +
#    geom_tile(aes(fill = value),colour = "white") + 
    geom_point(aes(color = value), size =3)  +
    scale_color_gradient2(low = "dodgerblue2", mid = "grey90",high = "red1", 
                          midpoint = 0, limits=c(-1.5, 1.5), oob=squish, na.value = "grey40") #+
    #theme_bw()
  return(p)
}

#plot heatmap of DEGs with circles and different sizes
plot_heatmap_genes_shape_size <- function(df){
  #  df.dendro <- as.dendrogram(hclust(d = dist(x = df)))
  #  df.order <- order.dendrogram(df.dendro)
  
  df$rownames <- rownames(df)
  df.m <- melt(df)
  
  # Order the levels according to dendogram order
  df.m$rownames <- fct_rev(factor(x = df.m$rownames,
                                  levels = rownames(df), 
                                  ordered = TRUE))
  
  df.m$abs <- abs(df.m$value)
  #  df.m <- ddply(df.m, .(variable), transform,rescale = rescale(value))
  p <- ggplot(df.m, aes(variable, rownames)) + 
    geom_point(aes(color = value, size =abs))  +
    geom_point(data = df.m[is.na(df.m$abs),], color = "grey40", size = 3) +
    scale_color_gradient2(low = "dodgerblue2", mid = "grey90",high = "red1", 
                          midpoint = 0, limits=c(-1.5, 1.5), oob=squish, na.value = "grey40") 
  #theme_bw()
  return(p)
}
