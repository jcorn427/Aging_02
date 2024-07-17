#Author Leanne Whitmore 
library(ggplot2)


theme_Publicationdot <- function(base_size = 14, base_family = "arial") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size)
    + theme(
      plot.title = element_text(
        face = "bold",
        size = rel(1.2), hjust = 0.5
      ),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text = element_text(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.key.size = unit(0.2, "cm"),
      legend.margin = unit(0, "cm"),
      legend.title = element_text( size= 9, face = "italic"),
      legend.text = element_text( size= 7, angle = 90, hjust = 1, vjust = 1),
      plot.margin = unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    ))
}

dotplot4bulkheatmap <- function(clusters, clustermatrix, pvals, breaks=NULL, width=5.5, 
                                labelsize=8, height=6, scalesize=8, figurename,sizebreaks=4) {
  df <- data.frame(matrix(ncol = length(colnames(clustermatrix)), 
                          nrow = length(unique(clusters))))
  colnames(df) <- colnames(clustermatrix)
  rownames(df) <-unique(clusters)
  dfpval <- data.frame(matrix(ncol = length(colnames(clustermatrix)), 
                              nrow = length(unique(clusters))))
  colnames(dfpval) <- colnames(clustermatrix)
  rownames(dfpval) <-unique(clusters)
  colnames(pvals) <- colnames(clustermatrix)
  print(unique(clusters))
  for (cluster in unique(clusters)) {
    print(cluster)
    genes = which(clusters==cluster)
    print(length(genes))
    for (col in colnames(clustermatrix)) {
      tmp = median(clustermatrix[names(genes), col])
      df[as.character(cluster), col] <- tmp
      countpvals =sum(pvals[names(genes), col] < 0.05)
      dfpval[as.character(cluster), col] <- countpvals
    }
  }
  print(df)
  dfmelt <- reshape2::melt(as.matrix(df))
  colnames(dfmelt)[3] <- "LFC"
  dfpvalsmelt <- reshape2::melt(as.matrix(dfpval))
  colnames(dfpvalsmelt)[3] <- "pval.count"
  if (all.equal(dfmelt$Var2, dfpvalsmelt$Var2)==TRUE) {
    dftotal <- cbind(dfmelt, dfpvalsmelt[,3])
    colnames(dftotal)[4] <- "pval.count"
  } else {
    stop("WARNING: DATA NOT IN THE SAME ORDER")
  }
  myPalette <- colorRampPalette(c("blue", "skyblue", "white", "orange", "red"))(100)
  sc <- scale_fill_gradientn(
    colours = myPalette,
    values = scales::rescale(c(
      min(dftotal$LFC), 0,
      0, max(dftotal$LFC)
    )), 
  )
  dftotal$LFC <- as.numeric(dftotal$LFC)
  dftotal$pval.count <- as.numeric(dftotal$pval.count)
  
  pl <- ggplot(data = dftotal, aes(
    x = Var2,
    y = factor(as.character(Var1),levels=as.character(rev(unique(clusters)))),
    size = pval.count
  )) + geom_point(aes(fill=LFC),
                  colour="black",pch=21) + theme_Publicationdot() + labs(fill="Median\nLFC", size="Adj. Pval\nCount") +
    scale_size_binned(range = c(1, scalesize), n.breaks = sizebreaks) + geom_vline(xintercept = breaks+.5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = labelsize)) +
    sc
  ggsave(figurename, height = height, width = width, dpi = 300)
}
