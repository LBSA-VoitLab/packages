#' rankedParallelPlotAuto Function
#'
#' This function allows you to heatmaps for genes specific to genesets
#' @param gsea_res_folder folder location of GSEA result
#' @keywords barplot gsea
#' @export
#' @examples
#' rankedParallelPlotAuto()

rankedParallelPlotAuto <- function(folder1, folder2, labels, top_n){
  df_1 <- get_gsea_df(folder1)
  df_2 <- get_gsea_df(folder2)
  df.1 <- df_1[1:top_n,]
  df.2 <- df_2[1:top_n,]
  df.new <- data.frame(pathway = df.1$NAME,rank=1:top_n,label = rep(labels[1],top_n))
  df.new <- rbind(df.new,data.frame(pathway = df.2$NAME,rank=1:top_n,label = rep(labels[1],top_n)))
  
  p <- rankedParallelPlot(df.new, labels)
  return(p)
}

#Default coloring
rankedParallelPlot <- function(df.new, labels){
  df.new$col <- sapply(df.new$pathway,function(x) {
    d=diff(df.new[df.new$pathway %in% x,"rank"])
    if (length(d) == 0) return(0)
    else return(sign(d)*(abs(d)+1)/25)})
  df.new$sz <- sapply(df.new$pathway,function(x) {
    d=diff(df.new[df.new$pathway %in% x,"rank"])
    if (length(d) == 0) return(0)
    else return((abs(d)+1)/25)})
  

  p<- ggplot(data = df.new,aes(x=label,y=-rank,group=pathway))+
    geom_line(aes(color = col, alpha = 1,size = sz)) +
    #  scale_color_brewer(palette="RdBu") +
    scale_colour_gradient2(low="red", high="blue", mid = gray(0.8)) +
    geom_point(aes(color = col, alpha = 1), size = 4) +
    geom_text(data = df.new %>% filter(label == labels[1]),
              aes(label = pathway) , hjust = 1, size = 3, nudge_x = -0.03) +
    geom_text(data = df.new %>% filter(label == labels[2]),
              aes(label = pathway) , hjust = 0, size = 3, nudge_x = 0.03) +
    theme(legend.position="none")
  return(p)
}

#Custom coloring column
rankedParallelPlot.2 <- function(df.new, labels, col_column){
  df.new$col <- sapply(df.new$pathway,function(x) {
    d=diff(df.new[df.new$pathway %in% x,"rank"])
    if (length(d) == 0) return(0)
    else return(sign(d)*(abs(d)+1)/25)})
  
  df.new$sz <- sapply(df.new$pathway,function(x) {
    d=diff(df.new[df.new$pathway %in% x,"rank"])
    if (length(d) == 0) return(0)
    else return((abs(d)+1)/25)})

  
  p<- ggplot(data = df.new,aes(x=label,y=-rank,group=pathway))+
    geom_line(aes(color = ancestor2, alpha = 1,size = sz)) +
    scale_color_manual(values = carto.pal(pal1 = "multi.pal", n1 = 20)) +
    #  scale_color_brewer(palette="RdBu") +
    geom_point(aes(, color = ancestor2, alpha = 1), size = 4) +
    scale_color_manual(values = carto.pal(pal1 = "multi.pal", n1 = 20)) +
#    scale_colour_gradient2(low="red", high="blue", mid = gray(0.8)) +
    geom_text(data = df.new %>% filter(label == labels[1]),
              aes(label = pathway, color = ancestor2) , hjust = 1, size = 3, nudge_x = -0.03) +
    geom_text(data = df.new %>% filter(label == labels[2]),
              aes(label = pathway, , color = ancestor2) , hjust = 0, size = 3, nudge_x = 0.03) +
    theme()
  
  p<-p + guides(color = guide_legend(reverse = TRUE))
  return(p)
}

#Default coloring with manual color pallete
rankedParallelPlot.3 <- function(df.new, labels){
  df.new$col <- sapply(df.new$pathway,function(x) {
    d=diff(df.new[df.new$pathway %in% x,"rank"])
    if (length(d) == 0) return(0)
    else return(sign(d)*(abs(d)+1)/25)})
  df.new$sz <- sapply(df.new$pathway,function(x) {
    d=diff(df.new[df.new$pathway %in% x,"rank"])
    if (length(d) == 0) return(0)
    else return((abs(d)+1)/25)})
  
  sub1 <- filter(df.new,label == labels[1])
  sub1$ord <- order(sub1$rank)
  sub2 <- filter(df.new,label == labels[2])
  sub2$ord <- order(sub2$rank)
  
  df.new <- rbind(sub1,sub2)
  
  col1 <- ggsci::pal_material("red")
  col2 <- ggsci::pal_material("blue")

#  col_man <- c(col1(10)[1:10],col2(10)[10:1])
  
  p<- ggplot(data = df.new,aes(x=label,y=-ord,group=pathway))+
    geom_line(aes(color = col, alpha = 0.5,size = sz)) +
    
    #  scale_color_brewer(palette="RdBu") +
    scale_colour_gradient2(low="red1", high="dodgerblue2", mid = gray(0.6)) +
#    scale_color_manual(values = col_man) +
#    scale_color_gsea() +
#    scale_color_material("red") +
    geom_point(aes(color = col, alpha = 0.5), size = 4) +
    geom_text(data = df.new %>% filter(label == labels[1]),
              aes(label = pathway) , hjust = 1, size = 5, nudge_x = -0.03) +
    geom_text(data = df.new %>% filter(label == labels[2]),
              aes(label = pathway) , hjust = 0, size = 5, nudge_x = 0.03) +
    theme(legend.position="none")
  
    return(p)
}


#Default coloring with manual color pallete and alluvium package for curves
rankedParallelPlot.4 <- function(df.new, labels){
  df.new$col <- sapply(df.new$pathway,function(x) {
    d=diff(df.new[df.new$pathway %in% x,"rank"])
    if (length(d) == 0) return(0)
    else return(sign(d)*(abs(d)+1)/25)})
  df.new$sz <- sapply(df.new$pathway,function(x) {
    d=diff(df.new[df.new$pathway %in% x,"rank"])
    if (length(d) == 0) return(0)
    else return((abs(d)+1)/25)})
  
  sub1 <- filter(df.new,label == labels[1])
  sub1$ord <- order(sub1$rank)
  sub2 <- filter(df.new,label == labels[2])
  sub2$ord <- order(sub2$rank)
  
  df.new <- rbind(sub1,sub2)
  
  col1 <- ggsci::pal_material("red")
  col2 <- ggsci::pal_material("blue")
  
  #  col_man <- c(col1(10)[1:10],col2(10)[10:1])
  
  p<- ggplot(data = df.new,aes(x=label,y=-ord,group=pathway))+
    geom_line(aes(color = col, alpha = 0.5,size = sz)) +
    
    #  scale_color_brewer(palette="RdBu") +
    scale_colour_gradient2(low="red1", high="dodgerblue2", mid = gray(0.6)) +
    #    scale_color_manual(values = col_man) +
    #    scale_color_gsea() +
    #    scale_color_material("red") +
    geom_point(aes(color = col, alpha = 0.5), size = 4) +
    geom_text(data = df.new %>% filter(label == labels[1]),
              aes(label = pathway) , hjust = 1, size = 5, nudge_x = -0.03) +
    geom_text(data = df.new %>% filter(label == labels[2]),
              aes(label = pathway) , hjust = 0, size = 5, nudge_x = 0.03) +
    theme(legend.position="none")
  
  return(p)
}
