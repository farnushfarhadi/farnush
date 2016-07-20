kmeanpp <- function(table , K )
{
  mean_names = c()
  means = c()
  idx_init = c()
  idx_init = c(idx_init , sample (1:dim(table)[1])[1] ); 
  mean_names = c(mean_names , table$tracking_id[idx_init] %>% as.character())
  means = rbind (means , table [ idx_init, - c(1,2) ] %>% as.numeric() )
  # ++ initialization
  
  for (i in 2:K)
  {
    dist = abs (cor (t(as.matrix(table[ , -c(1,2)]) ) , t(as.matrix(table[ idx_init, -c(1,2)])) )  )
    #rownames(dist ) = gene_names[-idx_init]
    minDist = apply(dist, 1, min) 
    idx_init = c(idx_init , which.max(minDist))
    means = rbind (means , table [ tail(idx_init,1), - c(1,2) ] %>% as.numeric() ) # choosing the mean that the minimum dist is maximized
    
  }
  rownames(means) = table$tracking_id[idx_init] %>% as.character()
  
  idx = 1 
  while (1)
  {
    print (paste ("iteration " , idx , sep=""))
    #assign each data point to closest mean
    dist = abs (cor (t(as.matrix(table[ , -c(1,2)]) )  , t(as.matrix(means))) )
    
    apply(dist, 1, function(x) {which.min(x)}) %>% unlist() -> clusters
    
    means_old = means
    # updating means
    for (i in 1:K)
    {
      t = as.matrix (table [which (clusters == i) , - c(1,2) ] )
      apply(t, 2, mean) -> m 
      means[i,] = m %>% unname()
    }
    
    if (max(max(abs(means_old - means))) < 1e-5)
      break;
  }
    
}