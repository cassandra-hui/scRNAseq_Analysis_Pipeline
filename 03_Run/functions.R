## Functions

# Convenience functions
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, name, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

ReadObject <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}


########################
## 10X


plot_cellranger_cells <- function(ind){
  xbreaks = c(1,1e1,1e2,1e3,1e4,1e5,1e6)
  xlabels = c("1","10","100","1000","10k","100K","1M")
  ybreaks = c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000)
  ylabels = c("1","2","5","10","2","5","100","2","5","1000","2","5","10k","2","5","100K","2","5","1M")
  
  pl1 <- data.frame(index=seq.int(1,ncol(d10x.data[[ind]])),
                    nCount_RNA = sort(Matrix:::colSums(d10x.data[[ind]])+1,decreasing=T),
                    nFeature_RNA = sort(Matrix:::colSums(d10x.data[[ind]]>0)+1,decreasing=T)) %>%
    ggplot() +
    scale_color_manual(values=c("red2","blue4"), labels=c("Features", "UMI"), name=NULL) +
    ggtitle(paste("CellRanger filltered cells:",ids[ind],sep=" ")) + xlab("Barcodes") + ylab("counts (UMI or Features") +
    scale_x_continuous(trans = 'log2', breaks=xbreaks, labels = xlabels) +
    scale_y_continuous(trans = 'log2', breaks=ybreaks, labels = ylabels) +
    geom_line(aes(x=index, y=nCount_RNA, color = "UMI"), size=1.75) +
    geom_line(aes(x=index, y=nFeature_RNA, color = "Features"), size=1.25)
  
  return(pl1)
}

