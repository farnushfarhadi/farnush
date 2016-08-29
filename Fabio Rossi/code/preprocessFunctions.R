### standardizing


logTransform_standardizeMyTable <- function(table , idx)
{
  
  #table <- failHandler(table)
  names <- colnames(table)
  
  t<- table [,idx:dim(table)[2]]
  r = dim(t)[1];
  c = dim(t)[2];
  t %>% as.matrix() -> p
  p %>% as.numeric() -> p1
  matrix(p1 , nrow = r , ncol = c) -> t
  log (t + 1) -> t 
  t.means <- colMeans(t)
  t.stdevs <- apply(t, 2, sd)
  
  for (i in 1:dim (t)[2])
  {
    t[,i] <-  (t[,i] - t.means[i]) / t.stdevs[i]
  }
  cbind (table[,1:(idx-1 )] , t) -> tableFinal
  colnames(tableFinal) <- names
  return(tableFinal)
}

logTransform <- function(table , idx)
{
  #table <- failHandler(table)
  (names <- colnames(table))
  
  t<- table [,idx:dim(table)[2]]
  r = dim(t)[1];
  c = dim(t)[2];
  t %>% as.matrix() -> p
  p %>% as.numeric() -> p1
  matrix(p1 , nrow = r , ncol = c) -> t
  log (t + 1) -> t 
  if (idx ==2)
  {
    table[,1:(idx-1)] %>% as.character() -> firstCol
    cbind (firstCol , t) -> tableFinal
  } else{
    cbind (table[,1:(idx-1)] , t) -> tableFinal
  }
  
  colnames(tableFinal) <- names
  return(tableFinal)
}
failHandler <- function (table)
{
  names <- colnames(table)
  
  t<- table [, idx:dim(table)[2]]
  as.matrix(t) -> t
  which (t == "FAIL") -> idx
  for (i in 1:length(idx))
  {
    as.integer(idx[i] / dim(t)[1]) -> p
    idx[i] - p * dim (t)[1] -> r
    if (r == 0)
    {
      t [idx[i]/p,p] <- "0"
    } else { 
      t[r ,p+1] <- "0"
    }
  }
  cbind (table[,1:(idx-1)] , t) -> tableFinal
  colnames(tableFinal) <- names
  return (tableFinal)
}

ss_corHeatmap <- function (table , toWrite , names, idx)
{
  cor_matrix <- cor (table [ , idx:dim(table)[2]] , table [ , idx:dim(table)[2]])
  #cor_matrix <- cor (t(table [ , idx:dim(table)[2]] ), t(table [ , idx:dim(table)[2]]) )
  colnames(cor_matrix) = rownames(cor_matrix) = names
  diag(cor_matrix) <- NA
  cols<-c(rev(brewer.pal(9,"YlOrRd")), "#FFFFFF")
  #cols<-colorRampPalette(brewer.pal(9,"Greens"))
  #cols<-colorRampPalette(brewer.pal(9,"Greens"))
  par(cex.main=0.8)
  # heatmap.2(cor_matrix, Rowv=NA, Colv=NA, symm=T, scale = NULL , trace="none", dendrogram="none", 
  #           col=cols, cexCol=0.7, cexRow=0.55 , margins=c(8,8) , srtCol=45  
  #           , main = toWrite )
  heatmap.2(cor_matrix,  symm=T, scale = NULL , trace="none",  
            col=cols, cexCol=0.6, cexRow=0.5 , margins=c(5,5) , srtCol=90 , 
            main = toWrite )
  # table[ , - c(1,2)] %>% as.matrix() -> toHeat
  # heatmap.2(t(toHeat), scale = NULL , trace="none",  
  #           col=cols, cexCol=0.6, cexRow=0.5 , margins=c(5,5) , srtCol=90 , 
  #           main = toWrite )
  #heatmap.2(cor_matrix , col=cols, cexCol=0.7, cexRow=0.55 , margins=c(8,8) , srtCol=45 , scale = NULL)
  
  #heatmap.2(cor_matrix,  symm=T, scale = NULL , trace="none", 
  #         col=cols, cexCol=0.7, cexRow=0.55 , margins=c(8,8) , srtCol=45  
  #          , main = toWrite )
} 
read.table("../../elana/exmaple 2/group1.txt" , header = FALSE) %>% unlist() -> g 
read.table("../../elana/exmaple 2/group2.txt" , header = FALSE) %>% unlist() %>% unname() -> g1
g3 <- read.table("../../elana/example 4/group2.txt") %>% unlist() %>% unname()
gg <- c( g , g1 , g3 )
sapply(gg , function(x) {which (FAP_damaged_log_filtered_n$tracking_id == x)}) -> idx
FAP_damaged_log_filtered_n[unname(idx) , - c(1,2)] -> part
rownames(part) <- FAP_damaged_log_filtered_n$tracking_id[unname(idx) ] %>% as.character()
cols<-c(rev(brewer.pal(9,"YlOrRd")), "#FFFFFF")

heatmap.2(as.matrix(part) , scale = NULL ,  Rowv = FALSE , Colv = TRUE , dendrogram = "column" , trace="none",  
          col=cols, cexCol=0.6, cexRow=0.5 , margins=c(5,5) , srtCol=90 , 
          main = "example 2" )


heatmap.2(as.matrix(part) , scale = NULL ,  Rowv = TRUE , Colv = FALSE , dendrogram = "row" , trace="none",  
          col=cols, cexCol=0.6, cexRow=0.5 , margins=c(5,5) , srtCol=90 , 
          main = "example 2" )

heatmap.2(as.matrix(part) , scale = NULL ,  Rowv = TRUE , Colv = TRUE , dendrogram = "both" , trace="none",  
          col=cols, cexCol=0.6, cexRow=0.5 , margins=c(5,5) , srtCol=90 , 
          main = "example 2" )

plotExpressionDist <- function (allSamples , t1 , t2)
{
  as.matrix(allSamples) -> allSamples
  as.numeric(allSamples) -> allSamples
  # the original range is 0 to 1M. hist is nt good so we transform. 
  log (allSamples + 1) -> allSamples
  par (mfrow = c(3 , 1))
  truehist(allSamples , prob = FALSE , main = "gene expression distribution" , 
           xlab="gene expression" , ylab = "density")
  axis(side=1, at=c(0:10))
  
  allSamples[which (allSamples >t1)] -> genes02
  truehist(genes02, main = "genes with >0.3 (log scale) expression distribution" , 
           xlab="gene expression" , ylab = "density" , prob = FALSE)
  axis(side=1, at=c(0:10 , by = 1))
  
  allSamples[which (allSamples >t2)] -> genes07
  truehist(genes07, main = "genes with >0.7 (log scale) expression distribution" , 
           xlab="gene expression" , ylab = "density" , prob = FALSE)
  axis(side=1, at=c(0:10 , by = 1))
}

filterLowExpressedGenes <- function(table  , threshold , percentage)
{
  if (percentage == 1)
  {
    apply(table[,-c(1,2)], 1, function(x) {all(as.numeric(x)>threshold )}) -> res
    which (res %>% unname() == FALSE ) -> lowIdx
    table [ - lowIdx , ] -> table 
    return (table)
  }
}


QuantileNormalize <- function (table , idx)
{
  r = dim (table)[1]
  c= dim (table)[2] - (idx-1) 
  table [ , - c(1:(idx-1))] %>% as.matrix() %>% as.numeric() %>% matrix (nrow = r, ncol = c) -> table_m
  normalize.quantiles(table_m , copy = TRUE) -> table_normal
  (colnames (table_normal ) <- colnames(table [ , - c(1:(idx-1))]))
  if (idx ==2)
  {
    table[,1:(idx-1)] %>% as.character() -> firstCol
    cbind (firstCol , table_normal) -> tableFinal
  } else{
    cbind (table[,1:(idx-1)] , table_normal) -> tableFinal
  }
  (colnames (tableFinal ) <- colnames(table ))
  return (tableFinal)
}


dataQC <- function (table , threshold , percentage , title , names , idx)
{
  logTransform (table , idx) -> table_log
  filterLowExpressedGenes(table_log , threshold , percentage) -> table_log_filtered 
  QuantileNormalize (table_log_filtered , idx) -> table_log_filtered_n 
  ss_corHeatmap (table_log_filtered_n , title , names , idx)
  return(as.data.frame( table_log_filtered_n) )
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


