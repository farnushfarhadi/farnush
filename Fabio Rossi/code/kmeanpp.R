
library (pdist)
kmeanpp <- function(table , K , maxIter , cor , modified , standardized)
{
  if (standardized)
  {
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

  mean_names = c()
  means = c()
  idx_init = c()
  (idx_init = c(idx_init , sample (1:dim(table)[1])[1] ));
  #(mean_names = c(mean_names , table$tracking_id[idx_init] %>% as.matrix() %>% as.character()) )
  (means = rbind (means , table [ idx_init, - c(1,2) ] %>% as.matrix() %>% as.numeric() ) )
  # ++ initialization
  
  print ("Initializing means...")
  for (i in 2:K)
  {
    if (i%% 50 ==0) print (paste (i , "means initialized" , sep = " "))
    
    
    if (cor )
    {
      #dist = abs (cor (t(as.matrix(table[ , -c(1,2)]) )  , t(as.matrix(means))) ) %>% round(4)
      dist = 1- (cor (t(as.matrix(table[ , -c(1,2)]) )  , t(as.matrix(means))) ) %>% round(4)
      # minDist = apply(dist, 1, function(x) {if (sort(x)[length(x)] == 1) {return(sort(x)[length(x)-1] )} else {
      #   return(sort(x)[length(x)] )} }) %>% unname()
      minDist = apply(dist, 1, min)
      (idx_init = c(idx_init , which.max(minDist[which (minDist != 0 )])[1]))
    } else {
      # Euclidean distance
      dist = pdist (as.matrix(table[ , -c(1,2)]) , as.matrix(means) ) %>% as.matrix()
      # minDist = apply(dist, 1, function(x) {if (sort(x)[1] == 0) {return(sort(x)[2] )} else {
      #   return(sort(x)[1] )} }) %>% unname()
      minDist = apply(dist, 1, min)
      (idx_init = c(idx_init , which.max(minDist[which (minDist != 0 )])[1]))
    }

    (means = rbind (means , table [ tail(idx_init,1), - c(1,2) ] %>% as.matrix() %>% as.numeric() ) )# choosing the mean that the minimum dist is maximized
    
  }
  #arownames(means) = table$tracking_id[idx_init] %>% as.character()
  
  idx = 1 
  while (1)
  {
    print (paste ("iteration " , idx , sep=""))
    #assign each data point to closest mean
    if (cor )
    {
      dist = 1- (cor (t(as.matrix(table[ , -c(1,2)]) )  , t(as.matrix(means))) ) %>% round(4)
      apply(dist, 1, function(x) {which.min(x)}) %>% unname() -> clusters
    } else {
      # Euclidean distance
      dist = pdist (as.matrix(table[ , -c(1,2)]) , as.matrix(means) ) %>% as.matrix()
      apply(dist, 1, function(x) {which.min(x)}) %>% unname() -> clusters
    }
  
    means_old = means
    # updating means
    for (i in 1:K)
    {
      if (modified)
      {
        if (which (clusters == i) %>% length() >1 )
        { # a cluster can be empty! 
          t = as.matrix (table [which (clusters == i) , - c(1,2) ] )
          # normalized mean
          norm.const = sum (t^2) %>% sqrt()
          (apply(t, 2, sum) -> m)  
          means[i,] = (m/norm.const) %>% unname()
          if (which(is.na(means)) %>% length() > 1)
          {
            print (paste ("updating mean problem: " , i , sep = ""))
            break;
            m
          }
        }
      }else{
        if (which (clusters == i) %>% length() >1 )
        { # a cluster can be empty! 
          t = as.matrix (table [which (clusters == i) , - c(1,2) ] )
          (apply(t, 2, mean) -> m)  # seems to remain same in different time points
          means[i,] = m %>% unname()
          if (which(is.na(means)) %>% length() > 1)
          {
            print (paste ("updating mean problem: " , i , sep = ""))
            break;
            m
          }
        }
      }

    }
    if (idx == maxIter)
      break;
    if (cor)
    {
      print (max(max(abs(means_old - means))) )
      if (max(max(abs(means_old - means))) < 1e-4)
        break;
    }else {
      print (max(max(abs(means_old - means))) )
      if (max(max(abs(means_old - means))) < 1e-4)
        break;
    }

    
    idx = idx +1
  }
  
  # BIC 
  bic = 0 ;
  for (i in 1:K)
  {
    # \sum_n (m_k(n)-x_n)^2  where m_k(n) is the mean associated to the nth data point
    (which (clusters == i) -> c_idx)
    m = dim (means)[2]
#     if (cor)
#     {
#       dist = (1- (cor (t(as.matrix(table[c_idx , -c(1,2)]) )  , means[i,])) %>% round(4) ) %>% sum()
#     }else {
#       dist = (pdist (as.matrix(table[ c_idx, -c(1,2)]) , means[i,] ) %>% as.matrix())^2 %>% sum()
#     }
    dist = (pdist (as.matrix(table[ c_idx, -c(1,2)]) , means[i,] ) %>% as.matrix())^2 %>% sum()
    bic = bic + dist;
    
  }
  bic = bic + log (dim(table)[1]) * K * dim(means)[2]
  
  ret = list (clusters , bic);
  names(ret) = c("clusters" , "BIC")
  return (ret)
    
}

## - question: optimal K
sapply (c(10:600 ) , kmeanpp , table= FAP_damaged_noRep_log_filtered ,maxIter = 200 , cor = 0) -> all_for_k
bic = c()
for (i in seq(2,length(all_for_k),2))
{
  bic = c(bic , all_for_k[[i]])
}
c(10:600) -> k 
plot ( k , bic , xlab = "Number of clusters" , ylab = "BIC score" , main = "K-mean with Euclidean distance" , col = "deepskyblue3")
# 42 is the best number  

kmeanpp (FAP_damaged_noRep_log_filtered , K = 500 , maxIter = 300 , cor = 0) -> FAP_damaged_kmean
kmeanpp (FAP_damaged_noRep_log_filtered , K = 500 , maxIter = 200 , cor = 1) -> FAP_damaged_kmean_cor
kmeanpp (FAP_damaged_noRep_log_filtered , K = 250 , maxIter = 300 , cor = 1) -> FAP_damaged_kmean_cor_k250
save (FAP_damaged_kmean_cor_k250 , file = "../../code/clustering res/FAP_damaged_kmean_cor_k250.RData")
load ("../../code/clustering res/FAP_damaged_kmean_cor_k250.RData")
sapply (seq (10 , 600 , 5 ) , kmeanpp , table= FAP_damaged_noRep_log_filtered ,maxIter = 400 , cor = 1) -> all_for_k_euclidean
plotBIC <- function (t , k , title)
{
  bic = c()
  for (i in seq(2,length(t),2))
  {
    bic = c(bic , t[[i]])
  }
  print(length(bic))
  plot ( k , bic , xlab = "Number of clusters" , ylab = "BIC score" , main = title, col = "deepskyblue3")
  #axis(side=1, at=c(k , by = 50))
}
pdf ("../../code/clustering res/kmean_fap_damaged_euclidean.pdf")
"../../code/clustering res/kmean_fap_damaged_pearson2.pdf"
table = FAP_damaged_noRep_log_filtered
witeplots <- function (path , clusters , table , z)
{
  pdf (path)
  for (i in 1:length(clusters))
  {
    table$tracking_id [which (clusters==i)] %>% as.character() -> genes
    if (length(genes) > 1 )
    {
      print (plotCluster(table , genes , paste ("cluster" , i , " ") , z) )
    }
  }
  dev.off()
}
witeplots("../../code/clustering res/kmean_fap_damaged_pearson_k250.pdf" ,FAP_damaged_kmean_cor_k250[[1]], FAP_damaged_noRep_log_filtered)

#### for Elena 
i = 93
clusters = FAP_damaged_kmean_cor_k250[[1]]
table = FAP_damaged_noRep_log_filtered
table$tracking_id [which (clusters==i)] %>% as.character() -> genes
which (clusters==i) -> idx
dist = cor ( t (as.matrix(table[idx , -c(1,2)]) )  )
colnames(dist) = rownames(dist) = genes
apply(dist, 1, function(x){which (x<0)}) -> sign_cor
pekh = c()
for (i in 1:length(sign_cor))
{
  pekh = union (pekh , sign_cor[[i]] %>% unname() )
}
genes [sign_cor[[1]] %>% unname()] -> l 
plotCluster(table , l , paste ("group1"))  -> p1
genes [-sign_cor[[1]] %>% unname()] -> l2
plotCluster(table , l2 , paste ("group2")) -> p2
multiplot (p1 , p2)
#example 2 is cluster 174 with k = 250
#example 3 is cluster 144 with k = 250
#example 4 is cluster 93 with k = 250
write.table (l , file = "../../elana/example 4/group1.txt", row.names = FALSE, col.names = FALSE , quote = FALSE)
write.table (l2 , file = "../../elana/example 4/group2.txt", row.names = FALSE, col.names = FALSE , quote = FALSE)
###
sapply(c (1:500) , function(x){which (FAP_damaged_kmean == x) %>% length()}) %>% unname() -> num
sapply(c (1:500) , function(x){which (FAP_damaged_kmean_cor == x) %>% length()}) %>% unname() -> num_cor

### - replicates as CV 
# kmean with one replicates and validate clusters in the other replicates
# EC WT - use 2-1 , 3-1 and 10-1 to find clusters 
EC_WT_train1 <- EC_WT_log_filtered_n[ , - c(5,6,8,9,15)]
sapply (seq (10 , 500 , 10) , kmeanpp , table= EC_WT_train1 ,maxIter = 200 , cor = 0) -> EC_WT_train1_k
sapply (seq (10 , 150 ,10 ) , kmeanpp , table= EC_WT_train1 ,maxIter = 200 , cor = 1 , modified = 0 , standardized = 0) -> EC_WT_train1_cor_k
plotBIC(EC_WT_train1_cor_k ,seq (10 , 150 ,10 )  , title = "EC WT k-mean 1-cor")
sapply (c(20:60) , kmeanpp , table= EC_WT_train1 ,maxIter = 200 , cor = 0) -> EC_WT_train1_k1

# works 
sapply (seq(50,150,4) , kmeanpp , table= EC_WT_train1 ,maxIter = 200 , cor = 0, modified = 0 , standardized = 1) -> EC_WT_train1_k_normalized_Euc
plotBIC (EC_WT_train1_k_normalized_Euc , seq(50,150,4) , "K-mean - standardized - Euclidean dist - classic centroid")
kmeanpp (table = EC_WT_train1 , K = 74 , maxIter = 200 , cor = 0 , modified = 0 , standardized = 1) -> EC_WT_train1_kmean_Euc_norm
witeplots("../../code/clustering res/kmean_EC_WT_train1_k74_normalized_Euc.pdf" ,EC_WT_train1_kmean_Euc_norm[[1]], EC_WT_train1 , z =1)
#sapply (seq(60,80,1) , kmeanpp , table= EC_WT_train1 ,maxIter = 200 , cor = 0 , modified = 0) -> l
plotBIC ( l , c(60:80))
sapply (seq(2,100,5) , kmeanpp , table= EC_WT_train1 ,maxIter = 200 , cor = 0 , modified = 1) -> EC_WT_train1_k_normalized_cModified
plotBIC( EC_WT_train1_k , seq (10 , 500 , 10))
plotBIC( EC_WT_train1_k1 , c(20:60) ) # 29
plotBIC( EC_WT_train1_k_normalized , seq(20,250,5) )
plotBIC( EC_WT_train1_k_normalized_cModified , seq(2,100,5) , title = "EC WT k-mean-standardized-modified centroid-Euclidean dist")
kmeanpp(EC_WT_train1 , K = 29 , maxIter = 200 , cor = 0) -> EC_WT_train1_kmean_Euc
kmeanpp(EC_WT_train1 , K = 73 , maxIter = 300 , cor = 0 , modified = 0) -> EC_WT_train1_kmean_norm
witeplots("../../code/clustering res/kmean_EC_WT_train1_k29.pdf" ,EC_WT_train1_kmean_Euc[[1]], EC_WT_train1)
witeplots("../../code/clustering res/kmean_EC_WT_train1_k29_normalized.pdf" ,EC_WT_train1_kmean_Euc[[1]], EC_WT_train1 , z =1)
witeplots("../../code/clustering res/kmean_EC_WT_train1_k70_normalized.pdf" ,EC_WT_train1_kmean_norm[[1]], EC_WT_train1 , z =1)
witeplots("../../code/clustering res/kmean_EC_WT_train1_k73_normalized.pdf" ,EC_WT_train1_kmean_norm[[1]], EC_WT_train1 , z =1)
minmin <- function (vector)
{
  min = 1000;
  for (i in seq(2 , length(vector) , 2))
    if (vector[[i]] < min)
      min = vector[[i]]
    
  return (i)
}
