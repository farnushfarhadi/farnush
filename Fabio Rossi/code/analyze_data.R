#labels 
## 'Question' for asking the question
## 'Resault' for showing the answer of questions

library (dplyr)
library(RColorBrewer)
library(ggplot2)
library(gplots)
library(MASS)
library(gdata)
library(reshape)
library (Rmisc)
##### DATA LOADING #####
#setwd ("~/Desktop/internship-ubc/data/txt version")
setwd ("~/Desktop/Fabio-summer 2016/data/txt version")

read.table("EC_wildType.txt" , header = TRUE , sep = "\t") -> EC_WT
read.table("EC_damaged.txt" , header = TRUE , sep = "\t") -> EC_damaged
 
read.table("EC_damaged.txt" , header = TRUE , sep = "\t" , stringsAsFactors = FALSE) -> EC_damaged2

read.table("FAP_wildType.txt" , header = TRUE , sep = "\t") -> FAP_WT
read.table("FAP_damaged.txt" , header = TRUE , sep = "\t") -> FAP_damaged

read.table("inflammatory_wildType.txt" , header = TRUE , sep = "\t") -> inflammatory_WT


read.table("muscleProgenitors_wildType.txt" , header = TRUE , sep = "\t") -> muscleProgenitors_WT
read.table("muscleProgenitors_damaged.txt" , header = TRUE , sep = "\t") -> muscleProgenitors_damaged

##### QC ####
- ## Question: same ids? what else? 
identical(EC_WT$tracking_id , EC_damaged$tracking_id)
identical(EC_WT$tracking_id , FAP_WT$tracking_id)
identical(EC_WT$tracking_id , FAP_damaged$tracking_id)
identical(EC_WT$tracking_id , inflammatory_WT$tracking_id)
identical(EC_WT$tracking_id , muscleProgenitors_WT$tracking_id)
identical(EC_WT$tracking_id , muscleProgenitors_damaged$tracking_id)

## Resault: YES. Ids are same . but there are some duplicated ids. Inflammatory and muscle does not have locus column

- ## Question: missing values? 
  
  # they are labeled by LOWDATA
 
apply (EC_WT , 2 , function(x){which(x == "LOWDATA")}) -> EC_WT_ld_idx
apply (EC_damaged , 2 , function(x){which(x == "LOWDATA")}) -> EC_damaged_ld_idx

apply (FAP_WT , 2 , function(x){which(x == "LOWDATA")}) -> FAP_WT_ld_idx
apply (FAP_damaged , 2 , function(x){which(x == "LOWDATA")}) -> FAP_damaged_ld_idx

apply (inflammatory_WT , 2 , function(x){which(x == "LOWDATA")}) -> inflammatory_WT_ld_idx

apply (muscleProgenitors_WT , 2 , function(x){which(x == "LOWDATA")}) -> muscleProgenitors_WT_ld_idx
apply (muscleProgenitors_damaged , 2 , function(x){which(x == "LOWDATA")}) -> muscleProgenitors_damaged_ld_idx
checkList <- function(list1 , list2)
{
  for (i in 3:length(list1))
  {
    for (j in 3:length(list2))
    {
      if (!identical(as.character(list1[i]) , as.character(list2[j])))
        print (paste (i , j , "FALSE!" , sep = " "))
    }
  }
}


checkList (EC_WT_ld_idx ,EC_damaged_ld_idx )
checkList (EC_WT_ld_idx ,FAP_WT_ld_idx )
checkList (EC_WT_ld_idx ,FAP_damaged_ld_idx )
checkList (EC_WT_ld_idx ,inflammatory_WT_ld_idx )
checkList (EC_WT_ld_idx ,muscleProgenitors_WT_ld_idx )
checkList (EC_WT_ld_idx ,muscleProgenitors_damaged_ld_idx )


## Resault: genes with LOWDATA (312 out of 23847) are same for all the samples (all replicates in different time points),
## so we can remove same genes from all samples
lowData_idx <- EC_WT_ld_idx[4] %>% unlist() %>% unname()

EC_WT[-lowData_idx , ] -> EC_WT
EC_damaged[-lowData_idx , ] -> EC_damaged
FAP_WT[-lowData_idx , ] -> FAP_WT
FAP_damaged[-lowData_idx , ] -> FAP_damaged
inflammatory_WT[-lowData_idx , ] -> inflammatory_WT
muscleProgenitors_WT[-lowData_idx , ] -> muscleProgenitors_WT
muscleProgenitors_damaged[-lowData_idx , ] -> muscleProgenitors_damaged
## dim : 23535 X samples_number 3 last cells does not have the chromosome column

# - Question: which genes are failed?
findFailGenes <- function (table)
{
  which (table == "FAIL") -> failGenesIdx
  arrayInd(failGenesIdx, dim(table)) -> failGenesIdx
  genesFailed <- table[failGenesIdx [,1],1] %>% as.character() %>% unique()
  return (genesFailed)
}

findFailGenes(EC_WT) -> EC_WT_f
findFailGenes(EC_damaged) -> EC_damaged_f
findFailGenes(FAP_WT) -> FAP_WT_f
findFailGenes(FAP_damaged) -> FAP_damaged_f
findFailGenes(inflammatory_WT) -> inflammatory_WT_f
findFailGenes(muscleProgenitors_WT) -> muscleProgenitors_WT_f
findFailGenes(muscleProgenitors_damaged) -> muscleProgenitors_damaged_f
geneList <- list (EC_WT_f ,EC_damaged_f , FAP_WT_f , FAP_damaged_f , inflammatory_WT_f , muscleProgenitors_WT_f , muscleProgenitors_damaged_f )
geneListFailed <- Reduce(union,geneList) 
length(geneListFailed) # 17 genes
# result

# - Task: filtering genes with FAIL label
sapply(geneListFailed , function(x) { which (EC_WT$tracking_id == x)}) %>% unname() -> failGenesIdx
arrayInd(which (inflammatory_WT == "HIDATA"), dim(inflammatory_WT)) -> hidata_idx 
idxToRemove <- c(failGenesIdx,hidata_idx[1] )
EC_WT[-idxToRemove , ] -> EC_WT
EC_damaged[-idxToRemove , ] -> EC_damaged
FAP_WT[-idxToRemove , ] -> FAP_WT
FAP_damaged[-idxToRemove , ] -> FAP_damaged
inflammatory_WT[-idxToRemove , ] -> inflammatory_WT
muscleProgenitors_WT[-idxToRemove , ] -> muscleProgenitors_WT
muscleProgenitors_damaged[-idxToRemove , ] -> muscleProgenitors_damaged
# dim 23518 X time_points
# Result



### standardizing

logTransform_standardizeMyTable <- function(table)
{
  
  #table <- failHandler(table)
  names <- colnames(table)
  
  t<- table [,3:dim(table)[2]]
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
  cbind (table[,1:2] , t) -> tableFinal
  colnames(tableFinal) <- names
  return(tableFinal)
}

logTransform <- function(table)
{
  #table <- failHandler(table)
  names <- colnames(table)
  
  t<- table [,3:dim(table)[2]]
  r = dim(t)[1];
  c = dim(t)[2];
  t %>% as.matrix() -> p
  p %>% as.numeric() -> p1
  matrix(p1 , nrow = r , ncol = c) -> t
  log (t + 1) -> t 
  cbind (table[,1:2] , t) -> tableFinal
  colnames(tableFinal) <- names
  return(tableFinal)
}
failHandler <- function (table)
{
  names <- colnames(table)
  
  t<- table [, 3:dim(table)[2]]
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
  cbind (table[,1:2] , t) -> tableFinal
  colnames(tableFinal) <- names
  return (tableFinal)
}


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

ss_corHeatmap <- function (table , toWrite)
{
  cor_matrix <- cor (table [ , 3:dim(table)[2]] , table [ , 3:dim(table)[2]])
  colnames(cor_matrix) = rownames(cor_matrix) = colnames(table)[3:length(colnames(table))]
  diag(cor_matrix) <- NA
  cols<-c(rev(brewer.pal(9,"YlOrRd")), "#FFFFFF")
  #cols<-colorRampPalette(brewer.pal(9,"Greens"))
  #cols<-colorRampPalette(brewer.pal(9,"Greens"))
  par(cex.main=0.8)
  heatmap.2(cor_matrix, Rowv=NA, Colv=NA, symm=T, scale = NULL , trace="none", dendrogram="none", 
            col=cols, cexCol=0.7, cexRow=0.55 , margins=c(8,8) , srtCol=45  
            , main = toWrite )
  
  heatmap.2(cor_matrix )
  
  #heatmap.2(cor_matrix,  symm=T, scale = NULL , trace="none", 
  #         col=cols, cexCol=0.7, cexRow=0.55 , margins=c(8,8) , srtCol=45  
  #          , main = toWrite )
} 
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
QuantileNormalize <- function (table)
{
  r = dim (table)[1]
  c= dim (table)[2] -2 
  table [ , - c(1,2)] %>% as.matrix() %>% as.numeric() %>% matrix (nrow = r, ncol = c) -> table_m
  normalize.quantiles(table_m , copy = TRUE) -> table_normal
  colnames (table_normal ) <- colnames(table [ , - c(1,2)])
  return (cbind (table[ , c(1,2)] , table_normal))
}


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