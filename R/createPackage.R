#This script consists the basic steps required to build a package

install.packages("devtools")
install.packages("roxygen2")
library("devtools")


setwd("Z:/projects-bmed/projects/packages/R/")
create("binf.gsea.visualizations")
create("binf.gsea.utils")
create("binf.rna.utils")
create("binf.genes.utils")
create("general.utils")
create_package("binf.trrust")
create_package("binf.modular")
create_package("binf.cellMarker")

install("binf.gsea.utils")
install("binf.gsea.visualizations")
install("binf.rna.utils")
install("binf.genes.utils")
install("general.utils")
install("binf.modular")
setwd("Z:/projects-bmed/projects/packages/R/")
install("binf.trrust")
install("binf.cellMarker")