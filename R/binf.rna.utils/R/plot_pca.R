#' plot_pca_deseq2
#'
#' This function allows to perform deseq2 diff exp analysis
#' @param gene_set_groups
#' @keywords GSEA gmt
#' @export
#' @examples
#' plot_pca_deseq2()

plot_pca_deseq2 <- function(rld){
  data <- plotPCA(rld, intgroup=c("Subject", "Status"), returnData=TRUE)
  percentVar <- round(100 * attr(data, "percentVar"))
  p <- ggplot(data, aes(PC1, PC2, color=Status)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))

  ### geom_label_repel
  q <- p +
    geom_label_repel(aes(label = name),
                     box.padding   = 0.15,
                     point.padding = 0.3,
                     segment.color = 'grey50') +
    theme_classic()
  return(q)
}

#pca with labels
plot_pca_df <- function(df,samples){
  pc <- prcomp(t(df))
  df_out <- as.data.frame(pc$x)
  df_out$group <- samples$Status
  df_out$group2 <- samples$Exp
  df_out$name <- samples$name

  percentage <- round(pc$sdev / sum(pc$sdev) * 100, 2)
  percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

  p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group,shape=group2))
  p<-p+geom_point(size= 3) + xlab(percentage[1]) + ylab(percentage[2])
  ### geom_label_repel
  p <- p +
    geom_label_repel(aes(label = name),
                     box.padding   = 0.15,
                     point.padding = 0.3,
                     segment.color = 'grey50') +
    theme_classic()
  return(p)
}

#pca without labels
plot_pca_df.2 <- function(df,samples){
  pc <- prcomp(t(df))
  df_out <- as.data.frame(pc$x)
  df_out$group <- samples$Status
  df_out$group2 <- samples$Exp
  df_out$group3 <- paste0(samples$Exp,samples$Status)
  df_out$name <- samples$name

  percentage <- round(pc$sdev / sum(pc$sdev) * 100, 2)
  percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

  p <- ggplot(df_out,aes(x=PC1,y=PC2,color=group3,shape=group)) +
    scale_shape_manual(values=c(16, 2, 17, 2))+
    scale_color_manual(values=c('#8367e6', '#6793e6', '#d767e6', '#D9AC6F', '#d98f6f','#d9c46f'))
  p <- p + geom_point(size= 3) + xlab(percentage[1]) + ylab(percentage[2])
  ### geom_label_repel
  p <- p +
    theme_classic()
  return(p)
}
