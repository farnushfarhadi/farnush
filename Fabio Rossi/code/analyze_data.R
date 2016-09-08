#labels 
## 'Question' for asking the question
## 'Resault' for showing the answer of questions
# a change to test for comit

## LOAD Data from here

setwd("Fabio Rossi/data/saved in R/")
load (file = "EC_WT.RData")
load (file = "EC_damaged.RData")
load (file = "FAP_WT.RData")
load (file = "FAP_damaged.RData")
load (file = "inflammatory_WT.RData")
load (file = "muscleProgenitors_WT.RData")
load (file = "muscleProgenitors_damaged.RData")

## we handle fail values, log transform and then standardize tables
# logTransform_standardizeMyTable(EC_WT) -> EC_WT_processed
# logTransform(EC_WT) -> EC_WT_log
# 
# logTransform_standardizeMyTable(EC_damaged) -> EC_damaged_processed
# logTransform(EC_damaged) -> EC_damaged_log
# logTransform_standardizeMyTable(FAP_WT) -> FAP_WT_processed
# logTransform(FAP_damaged) -> FAP_damaged_log
# logTransform_standardizeMyTable(FAP_damaged) -> FAP_damaged_processed
# logTransform_standardizeMyTable(inflammatory_WT) -> inflammatory_WT_processed
# logTransform_standardizeMyTable(muscleProgenitors_WT) -> muscleProgenitors_WT_processed
# logTransform_standardizeMyTable(muscleProgenitors_damaged) -> muscleProgenitors_damaged_processed

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


# - Question: finding low expressed genes with the distribution. - all samples 
cbind (EC_WT[ , -c(1,2)] , EC_damaged[ , -c(1,2)] , FAP_WT[,-c(1,2)] , FAP_damaged[,-c(1,2)] , 
inflammatory_WT[,-1] , muscleProgenitors_WT[,-1] , muscleProgenitors_damaged[,-1]) -> allSamples


plotExpressionDist (allSamples , 0.3 , 0.7)
# Result: 1 as threshold in log scale 

# - Question: finding low expressed genes and normalizing each cell individually
plotExpressionDist (EC_damaged[ , - c(1,2)] , 0.3 , 0.6)
# - Question: filtering low expressed genes
# threshold is 0.7 log scale - removing genes less than threshold in at least one time points

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
EC_WT_days <- c("D0" , "D2-1" , "D2-2" , "D2-3" , "D3-1" , "D3-2" , "D3-3" , "D4" , "D5" , "D6" , "D7" , "D10-1" , "D10-2" , "D14" )
dataQC (EC_WT , 0.7 , 1 , "EC WT-quantile normalized" , EC_WT_days , 3) -> EC_WT_log_filtered_n # 9396 16 genes
EC_damaged_days <- c("D0" , "D1" , "D2" , "D3" , "D5" , "D6" , "D7" , "D10")
dataQC (EC_damaged , 0.7 , 1 , "EC KO-quantile normalized" , EC_damaged_days , 3) -> EC_damaged_log_filtered_n # 9025 10 genes
FAP_WT_days <- c("D0" , "D1" , "D2-1" , "D2-2" , "D2-3" , "D3-1" , "D3-2" , "D3-3" , "D4-1" , "D4-2" , "D5" , "D6" , "D7" , "D10" , "D14" )
dataQC (FAP_WT , 0.7 , 1 , "FAP WT-quantile normalized" , FAP_WT_days , 3) -> FAP_WT_log_filtered_n # 8909 17 genes
FAP_damaged_days <- c("D0" , "D1" , "D2" , "D3-1" , "D3-2" , "D4" , "D5" , "D6"  , "D10")
dataQC (FAP_damaged , 0.7 , 1 , "FAP KO-quantile normalized" , FAP_damaged_days ,3) -> FAP_damaged_log_filtered_n # 8341 10 genes
muscleProgenitors_WT_days <- c( "D1" , "D2" , "D3-1" , "D3-2" , "D5" ,  "D7"  , "D10")
dataQC (muscleProgenitors_WT , 0.7 , 1 , "muscleProgenitors WT-quantile normalized", muscleProgenitors_WT_days,2) -> muscleProgenitors_WT_log_filtered_n # 9014 8 genes
muscleProgenitors_damaged_days <- c( "D0" , "D3-1" , "D3-2" ,"D4" , "D5" ,"D6"  , "D10")
dataQC (muscleProgenitors_damaged , 0.7 , 1 , "muscleProgenitors KO-quantile normalized" ,muscleProgenitors_damaged_days, 2) -> muscleProgenitors_damaged_log_filtered_n # 9535 8 genes
inflammatory_WT_days <- c( "D1-1", "D1-2" ,"D2-1" , "D2-2" , "D2-3" , "D3-1" , "D3-2" , "D3-3" , "D4"  ,"D5" ,"D6","D7-1" ,"D7-2" , "D10")
dataQC (inflammatory_WT , 0.7 , 1 , "inflammatory WT-quantile normalized", inflammatory_WT_days,2) -> inflammatory_WT_log_filtered_n # 2263 15 genes



##################
# bulding tables for different cell types 
EC_WT_log_filtered_n_rep1 = EC_WT_log_filtered_n[ , c(1,2,3,4,7,10,11,12,13,14,16)]
EC_WT_days_rep1 <- c("D0" , "D2-1" , "D3-1"  , "D4" , "D5" , "D6" , "D7" , "D10-1"  , "D14" )
EC_KO_log_filtered_n_rep1 = EC_damaged_log_filtered_n
EC_KO_days_rep1 <- c("D0" , "D1" , "D2" , "D3" , "D5" , "D6" , "D7" , "D10")
FAP_WT_log_filtered_n_rep1 = FAP_WT_log_filtered_n[ , c(1,2,3,4,5,8,11,13,14,15,16,17)]
FAP_WT_days_rep1 <- c("D0" ,"D1" , "D2-1" , "D3-1"  , "D4-1" , "D5" , "D6" , "D7" , "D10"  , "D14" )
FAP_KO_log_filtered_n_rep1 = FAP_damaged_log_filtered_n[ , c(1,2,3,4,5,6,8,9,10,11)]
FAP_KO_days_rep1 <- c("D0" , "D1" , "D2" , "D3-1" , "D4" , "D5" , "D6"  , "D10")

logTransform (FAP_damaged_noRep) -> FAP_damaged_noRep_log
filterLowExpressedGenes(FAP_damaged_noRep_log , 0.7 , 1) -> FAP_damaged_noRep_log_filtered # 8394 genes
QuantileNormalize (FAP_damaged_noRep_log_filtered) -> l 

#heatmap of all samples 
cbind (EC_WT[, c(1,2)] , allSamples) -> allSamplestoHeat
logTransform(allSamplestoHeat , 3) -> allSamplestoHeat_log
filterLowExpressedGenes(allSamplestoHeat_log , 0.7 , 1) -> allSamplestoHeat_filtered # 2263 15 genes
QuantileNormalize (allSamplestoHeat_filtered) -> allSamplestoHeat_filtered_n 
allSamples_names = c(paste("EC_WT",EC_WT_days),paste("EC_KO",EC_damaged_days),
                     paste("FAP_WT",FAP_WT_days),paste("FAP_WT",FAP_damaged_days),
                     paste("ifml_WT",inflammatory_WT_days),
                     paste("mus_WT",muscleProgenitors_WT_days),paste("mus_KO",muscleProgenitors_damaged_days))
ss_corHeatmap (allSamplestoHeat_filtered_n , "all samples-quantile normalized" ,allSamples_names,3)

sample <- allSamplestoHeat_filtered_n[, - c(1,2)]
colnames(sample) <- allSamples_names
cols<-c(rev(brewer.pal(9,"YlOrRd")), "#FFFFFF")
par(cex.main=0.8)
heatmap.2(sample%>% as.matrix() %>% as.numeric() %>% matrix(nrow= dim(sample)[1] , ncol = dim(sample)[2] ),  scale = NULL , trace="none",  
          col=cols, cexCol=0.6, cexRow=0.5 , margins=c(5,5) , srtCol=90 , Colv = TRUE , Rowv = FALSE,
          main = "genes - all samples" )



logTransform (muscleProgenitors_damaged_noRep) -> muscleProgenitors_damaged_noRep_log


filterLowExpressedGenes(EC_damaged_log , 0.7 , 1) -> EC_damaged_log_filtered # 9025 genes
filterLowExpressedGenes(muscleProgenitors_damaged_noRep_log , 0.7 , 1) -> muscleProgenitors_damaged_noRep_log_filtered # 9623