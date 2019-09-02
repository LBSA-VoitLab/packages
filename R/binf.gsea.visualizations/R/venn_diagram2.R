#' plot_venn_diagram2 Function
#'
#' This function allows you to create venn diagram for 2 sets
#' @param list1 list 1
#' @param list2 list 2
#' @keywords venn diagram
#' @export
#' @examples
#' plot_venn_diagram2()

library(VennDiagram)

plot_venn_diagram2 <- function(list1,list2,cat){
  if(missing(cat)) {
    c <- c("1","2")
  } else {
    c <- cat
  }
  grid.newpage()
  p <-  draw.pairwise.venn(length(list1), length(list2), length(intersect(list1,list2)), category = c,fill = c("light blue", "pink"))
#  draw.pairwise.venn(5, 6, 4, category = c)
  return(p)
}
