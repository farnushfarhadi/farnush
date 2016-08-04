#method (EC_damaged_log_filtered , 0.95 ) -> res95_notRep
load ("../../code/clustering res/method1/90EC_damaged_noRep.RData")
load("../../code/clustering res/method1/95EC_damaged_noRep.RData")
res95_notRep_num <- sapply(res95_notRep , function(x) {length (x)}) 
res95_notRep_num  %>% order() %>% tail (50) -> idx
res95_notRep_num [idx]
logTransform(EC_damaged) -> EC_damaged_log
logTransform(EC_WT) -> EC_damaged_log
plotCluster(EC_damaged_log , res95_notRep[[tail(idx ,1)]] , "standardized" , 1) -> pn
plotCluster(EC_damaged_log , res95_notRep[[tail(idx ,1)]] , "original" , 0) -> p
multiplot(pn , p , cols = 2)




plotPairedKOvsWT <- function (wt_table , ko_table , cluster_res , wt_days , ko_days , num , path, jj)
{
  cluster_res_num <- sapply(cluster_res , function(x) {length (x)}) 
  cluster_res_num  %>% order() %>% tail (num) -> idx
  
  pdf(file=path) 
  for (i in num:1)
  {
    print(i)
    par(mfrow = c(2,2))
    print (plotCluster( wt_table , cluster_res[[  idx[i] ]]  , paste("WT-cluster" , i , sep = "-") , wt_days , jj) )
    print (plotCluster( ko_table , cluster_res[[  idx[i] ]]  , paste( "KO-cluster" , i , sep = "-") , ko_days , jj) )
  }
  dev.off()
}

plotPairedKOvsWT_v2 <- function (wt_table , ko_table , cluster_res , wt_days , ko_days , num , path, jj)
{
  cluster_res_num <- sapply(cluster_res , function(x) {length (x)}) 
  cluster_res_num  %>% order() %>% tail (num) -> idx
  
  pdf(file=path) 
  for (i in num:1)
  {
    print(i)
    p1 <- plotCluster( wt_table , cluster_res[[  idx[i] ]]  , paste("WT-cluster" , i , sep = "-") , wt_days , jj) 
    p2 <- plotCluster( ko_table , cluster_res[[  idx[i] ]]  , paste( "KO-cluster" , i , sep = "-") , ko_days , jj) 
    multiplot(p1 , p2 , rows = 2)
  }
  dev.off()
}

## looking at WT and KO patterns 
ec_path = "../../code/clustering res/method1/EC/WT_KO_patterns.pdf"
plotPairedKOvsWT (EC_WT_log_filtered_n_1 , EC_damaged_log_filtered_n , EC_WT_m1_80 , EC_WT_days1 , EC_damaged_days , 200 ,ec_path ,3)
fap_path = "../../code/clustering res/method1/FAP/WT_KO_patterns.pdf"
plotPairedKOvsWT (FAP_WT_log_filtered_n_1 , FAP_damaged_log_filtered_n_1 , FAP_WT_m1_80 , FAP_WT_days1 , FAP_ko_days1 , 200 ,fap_path ,3)
muscle_path = "../../code/clustering res/method1/muscle/WT_KO_patterns.pdf"
plotPairedKOvsWT (muscleProgenitors_WT_log_filtered_n_1 , muscleProgenitors_damaged_log_filtered_n_1 , 
                  muscle_WT_m1_80 , muscle_WT_days1 , muscle_ko_days1 , 300 , muscle_path ,2)
