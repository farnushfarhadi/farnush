corAverage <- function (table  , genes , z)
{
  n = length(genes)
  sapply(genes,function(x) {which (table$tracking_id == x)})  -> c_idx
  if (z)
  { ### NOTE: It does not matter if you normalize or not. cor is same 
    t = table [ , -c(1,2)] %>% as.matrix() 
    cluster.means <- apply(t, 1, mean)
    cluster.stdevs <- apply(t, 1, sd)
    
    for (i in 1:dim (t)[1])
    {
      t[i,] <-  (t[i,] - cluster.means[i]) / cluster.stdevs[i]
    }
    colnames(t) = colnames(table)[-c(1,2)]
    table = cbind (table [ , c(1,2)] , t)
  }

  cor (t(as.matrix(table[Reduce(union , c_idx) , -c(1,2)]) )) %>% round(4) -> cor_scores
  avg_cor <- (cor_scores[lower.tri(cor_scores)] %>% sum() ) / (n * (n-1) )
  #return (avg_cor)
  return (cor_scores)
}

doCorAverageForClusters <- function (clusters , table , z)
{
  cor_avg = c()
  for (i in 1:length(unique (clusters)) )
  {
    table$tracking_id [which (clusters==i)] %>% as.character() -> genes
    if (length(genes) > 1 )
    {
      cor_avg = c(cor_avg , corAverage(table , genes , z) )
    }
  }
  return(cor_avg)
}

doCorAverageForClusters (clusters = EC_WT_train1_kmean_Euc_norm[[1]] , table = EC_WT_train1 , z = 1) -> l
# EC WT - use 2-2 , 3-2 and 10-2 to find clusters 
EC_WT_train1 <- EC_WT_log_filtered_n[ , - c(5,6,8,9,15)]
EC_WT_train2 <- EC_WT_log_filtered_n[ , - c(4,6,7,9,14)]
doCorAverageForClusters (clusters = EC_WT_train1_kmean_Euc_norm[[1]] , table = EC_WT_train2) 



