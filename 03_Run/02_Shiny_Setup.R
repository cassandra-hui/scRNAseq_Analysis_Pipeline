if (!any(rownames(installed.packages()) == "shiny")){
  install.packages("shiny")
}
library(shiny)

if (!any(rownames(installed.packages()) == "markdown")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("markdown")
}
library(markdown)

if (!any(rownames(installed.packages()) == "tidyr")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("tidyr")
}
library(tidyr)

download.file("https://github.com/ucdavis-bioinformatics/scRNA_shiny/archive/master.zip", "scRNA_shiny.zip")
#zipf <- file.choose(new=FALSE)
zipf <- "scRNA_shiny.zip"
outdir <- "scRNA_shiny"
unzip(zipf, exdir=outdir)


