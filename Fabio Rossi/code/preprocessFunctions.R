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
  names <- colnames(table)
  
  t<- table [,idx:dim(table)[2]]
  r = dim(t)[1];
  c = dim(t)[2];
  t %>% as.matrix() -> p
  p %>% as.numeric() -> p1
  matrix(p1 , nrow = r , ncol = c) -> t
  log (t + 1) -> t 
  cbind (table[,1:(idx-1)] , t) -> tableFinal
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
  #heatmap.2(cor_matrix , col=cols, cexCol=0.7, cexRow=0.55 , margins=c(8,8) , srtCol=45 , scale = NULL)
  
  #heatmap.2(cor_matrix,  symm=T, scale = NULL , trace="none", 
  #         col=cols, cexCol=0.7, cexRow=0.55 , margins=c(8,8) , srtCol=45  
  #          , main = toWrite )
} 

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


QuantileNormalize <- function (table)
{
  r = dim (table)[1]
  c= dim (table)[2] -2 
  table [ , - c(1,2)] %>% as.matrix() %>% as.numeric() %>% matrix (nrow = r, ncol = c) -> table_m
  normalize.quantiles(table_m , copy = TRUE) -> table_normal
  colnames (table_normal ) <- colnames(table [ , - c(1,2)])
  return (cbind (table[ , c(1,2)] , table_normal))
}


dataQC <- function (table , threshold , percentage , title , names , idx)
{
  logTransform (table , idx) -> table_log
  filterLowExpressedGenes(table_log , threshold , percentage) -> table_log_filtered 
  QuantileNormalize (table_log_filtered) -> table_log_filtered_n 
  ss_corHeatmap (table_log_filtered_n , title , names , idx)
  return(table_log_filtered_n)
}

