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
EC_WT_m1_80_size <- sapply(EC_WT_m1_80 , function(x) {length (x)}) 
EC_WT_m1_80_size  %>% order() %>% tail (100) -> idx
i = 98
ii <- idx[101-i]

days = c("D0" , "D1" , "D2"  , "D3" , "D5" , "D6" , "D7" , "D10" )
plotCluster(EC_damaged_log_filtered , EC_WT_m1_80[[ii]] , paste( "cluster", i  , sep = " "),days )
plotCluster(EC_WT_log_filtered_n_1 ,EC_WT_m1_80[[ii]] , paste( "cluster", i  , sep = " "),EC_WT_days1  )


load("../../code/clustering res/method1/EC/EC_WT90.RData")
EC_WT_m1_90_size <- sapply(EC_WT_m1_90 , function(x) {length (x)}) 
EC_WT_m1_90_size  %>% order() %>% tail (100) -> idx
i = 48
ii <- idx[101-i]

days = c("D0" , "D1" , "D2"  , "D3" , "D5" , "D6" , "D7" , "D10" )
plotCluster(EC_damaged_log_filtered , EC_WT_m1_90[[ii]] %>% unique(), paste( "cluster", i  , sep = " "),days )
plotCluster(EC_WT_log_filtered_n_1 ,EC_WT_m1_90[[ii]] , paste( "cluster", i  , sep = " "),EC_WT_days1  )

