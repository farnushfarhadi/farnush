#labels 
## 'Question' for asking the question
## 'Resault' for showing the answer of questions


## LOAD Data from here
load (file = "../saved in R /EC_WT.RData")
load (file = "../saved in R /EC_damaged.RData")
load (file = "../saved in R /FAP_WT.RData")
load (file = "../saved in R /FAP_damaged.RData")
load (file = "../saved in R /inflammatory_WT.RData")
load (file = "../saved in R /muscleProgenitors_WT.RData")
load (file = "../saved in R /muscleProgenitors_damaged.RData")

## we handle fail values, log transform and then standardize tables
logTransform_standardizeMyTable(EC_WT) -> EC_WT_processed
logTransform(EC_WT) -> EC_WT_log

logTransform_standardizeMyTable(EC_damaged) -> EC_damaged_processed
logTransform(EC_damaged) -> EC_damaged_log
logTransform_standardizeMyTable(FAP_WT) -> FAP_WT_processed
logTransform(FAP_damaged) -> FAP_damaged_log
logTransform_standardizeMyTable(FAP_damaged) -> FAP_damaged_processed
logTransform_standardizeMyTable(inflammatory_WT) -> inflammatory_WT_processed # HIDATA !!!
# Ubb gene has HIDATA value
logTransform_standardizeMyTable(muscleProgenitors_WT) -> muscleProgenitors_WT_processed
logTransform_standardizeMyTable(muscleProgenitors_damaged) -> muscleProgenitors_damaged_processed

- ## Question: sample-sample correlation heatmap 


ss_corHeatmap(EC_WT_processed , "EC_WT")
ss_corHeatmap(EC_damaged_processed , "EC_damaged")
ss_corHeatmap(FAP_WT_processed , "FAP_WT")
ss_corHeatmap (FAP_damaged_processed , "FAP_damaged")
ss_corHeatmap(muscleProgenitors_WT_processed , "muscle progenitors_WT")
ss_corHeatmap(muscleProgenitors_damaged_processed , "muscle progenitors_damaged")
## Result: 
  

- # Question: dist. of genes in one sample!? 
  
hist (EC_WT_log$Endothelium.0D.NTX[which (EC_WT_log$Endothelium.0D.NTX >0.5) ])
hist (EC_WT_log$Endothelium.0D.NTX)
  


#### 
# issues to talk 
# FAIL
# HIDATA
# duplicated genes 



#STEM
read.table("../../stem/EC.txt" , header = TRUE , sep= "\t") -> stem_dat
stem_dat[-lowData_idx , ] -> stem_dat
for (i in 1:dim(stem_dat)[2])
{
  as.matrix(stem_dat) -> stem_dat
  (which (stem_dat [ , i] == "FAIL")  -> idx)
  if (length(idx) >0)
    stem_dat [idx , i] <- ""
}
write.table (stem_dat , file = "~/Desktop/internship-ubc/stem/EC_damaged_processed.txt"  , sep = "\t" , quote =FALSE , row.names = FALSE)


# - Question: finding low expressed genes with the distribution. - all samples 
cbind (EC_WT[ , -c(1,2)] , EC_damaged[ , -c(1,2)] , FAP_WT[,-c(1,2)] , FAP_damaged[,-c(1,2)] , 
inflammatory_WT[,-1] , muscleProgenitors_WT[,-1] , muscleProgenitors_damaged[,-1]) -> allSamples


plotExpressionDist (allSamples , 0.3 , 0.7)
# Result: 1 as threshold in log scale 

# - Question: finding low expressed genes and normalizing each cell individually
plotExpressionDist (EC_damaged[ , - c(1,2)] , 0.3 , 0.6)
# - Question: filtering low expressed genes
# threshold is 0.7 log scale - removing genes less than threshold in at least one time points

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

# - Question: how do we treat biological replicates
# averaging over biological replicates
# EC_damaged is fine
# FAP_damaged day 3 has two replicates
averageReplicates <- function (table , )
FAP_damaged  -> FAP_damaged_noRep
(FAP_damaged_noRep[ , 6] %>% as.matrix() %>% as.numeric()+ FAP_damaged_noRep [,7] %>% as.matrix() %>% as.numeric()) /2 -> FAP_damaged_noRep[,6]
FAP_damaged_noRep[ , -7] -> FAP_damaged_noRep

# Muscle progenitors damaged day 3 has two replicates
muscleProgenitors_damaged -> muscleProgenitors_damaged_noRep 
(muscleProgenitors_damaged_noRep[ , 3] %>% as.matrix() %>% as.numeric()+ muscleProgenitors_damaged_noRep[,4] %>% as.matrix() %>% as.numeric()) /2 -> muscleProgenitors_damaged_noRep[,3]
muscleProgenitors_damaged_noRep[ , -4] -> muscleProgenitors_damaged_noRep



# - z transfomr for cor?
#- repeated genes 
# var filtering genes

# - question: quantile normalization 

#checking samples distribution 
# par (mfrow = c(2, 3))
# hist (log (FAP_damaged_m[,1] + 1 ) , col = "red")
# hist (log (FAP_damaged_m[,2] + 1 ) , col = "blue")
# hist (log (FAP_damaged_m[,3] + 1 ) , col = "green")
# hist (log (FAP_damaged_normal[,1] + 1 ) , col = "red")
# hist (log (FAP_damaged_normal[,2] + 1 ) , col = "blue")
# hist (log (FAP_damaged_normal[,3] + 1 ) , col = "green")



# - Question: preprocessing  data
logTransform (EC_WT) -> EC_WT_log
filterLowExpressedGenes(EC_WT_log , 0.7 , 1) -> EC_WT_log_filtered # 9396 16 genes
QuantileNormalize (EC_WT_log_filtered) -> EC_WT_log_filtered_n 
ss_corHeatmap (EC_WT_log_filtered_n , "EC WT-quantile normalized")
logTransform (EC_damaged) -> EC_damaged_log
filterLowExpressedGenes(EC_damaged_log , 0.7 , 1) -> EC_damaged_log_filtered # 9025 10 genes
QuantileNormalize (EC_damaged_log_filtered) -> EC_damaged_log_filtered_n 
ss_corHeatmap (EC_damaged_log_filtered_n , "EC KO-quantile normalized")

logTransform (FAP_WT) -> FAP_WT_log
filterLowExpressedGenes(FAP_WT_log , 0.7 , 1) -> FAP_WT_log_filtered # 8909 17 genes
QuantileNormalize (FAP_WT_log_filtered) -> FAP_WT_log_filtered_n 
ss_corHeatmap (FAP_WT_log_filtered_n , "FAP WT-quantile normalized")
logTransform (FAP_damaged) -> FAP_damaged_log
filterLowExpressedGenes(FAP_damaged_log , 0.7 , 1) -> FAP_damaged_log_filtered # 8341 11 genes
QuantileNormalize (FAP_damaged_log_filtered) -> FAP_damaged_log_filtered_n 
ss_corHeatmap (FAP_damaged_log_filtered_n , "FAP KO-quantile normalized")

logTransform (muscleProgenitors_WT) -> muscleProgenitors_WT_log
filterLowExpressedGenes(muscleProgenitors_WT_log , 0.7 , 1) -> muscleProgenitors_WT_log_filtered # 9014 8 genes
QuantileNormalize (muscleProgenitors_WT_log_filtered) -> muscleProgenitors_WT_log_filtered_n 
ss_corHeatmap (muscleProgenitors_WT_log_filtered_n , "muscleProgenitors WT-quantile normalized")
logTransform (muscleProgenitors_damaged) -> muscleProgenitors_damaged_log
filterLowExpressedGenes(muscleProgenitors_damaged_log , 0.7 , 1) -> muscleProgenitors_damaged_log_filtered # 9535 8 genes
QuantileNormalize (muscleProgenitors_damaged_log_filtered) -> muscleProgenitors_damaged_log_filtered_n 
ss_corHeatmap (muscleProgenitors_damaged_log_filtered_n , "muscleProgenitors KO-quantile normalized")

logTransform (inflammatory_WT) -> inflammatory_WT_log
filterLowExpressedGenes(inflammatory_WT_log , 0.7 , 1) -> inflammatory_WT_log_filtered # 2263 15 genes
QuantileNormalize (inflammatory_WT_log_filtered) -> inflammatory_WT_log_filtered_n 
ss_corHeatmap (inflammatory_WT_log_filtered_n , "inflammatory WT-quantile normalized")

logTransform (FAP_damaged_noRep) -> FAP_damaged_noRep_log
filterLowExpressedGenes(FAP_damaged_noRep_log , 0.7 , 1) -> FAP_damaged_noRep_log_filtered # 8394 genes
QuantileNormalize (FAP_damaged_noRep_log_filtered) -> l 

#heatmap of all samples 
cbind (EC_WT[, c(1,2)] , allSamples) -> allSamplestoHeat
logTransform(allSamplestoHeat) -> allSamplestoHeat_log
filterLowExpressedGenes(allSamplestoHeat , 0.7 , 1) -> allSamplestoHeat_filtered # 2263 15 genes
QuantileNormalize (allSamplestoHeat_filtered) -> allSamplestoHeat_filtered_n 
ss_corHeatmap (allSamplestoHeat_filtered_n , "all samples-quantile normalized")


logTransform (muscleProgenitors_damaged_noRep) -> muscleProgenitors_damaged_noRep_log


filterLowExpressedGenes(EC_damaged_log , 0.7 , 1) -> EC_damaged_log_filtered # 9025 genes
filterLowExpressedGenes(muscleProgenitors_damaged_noRep_log , 0.7 , 1) -> muscleProgenitors_damaged_noRep_log_filtered # 9623