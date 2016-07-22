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


## looking at WT patterns 
## 