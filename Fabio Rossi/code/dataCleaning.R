library (dplyr)
library(RColorBrewer)
library(ggplot2)
library(gplots)
library(MASS)
library(gdata)
library(reshape)
library (Rmisc)
library(preprocessCore)
##### DATA LOADING #####
#setwd ("~/Desktop/internship-ubc/data/txt version")
setwd ("Fabio Rossi/data/txt version")

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

# saving data in R 
save (EC_WT , file = "../saved in R /EC_WT.RData")
save (EC_damaged , file = "../saved in R /EC_damaged.RData")
save (FAP_WT , file = "../saved in R /FAP_WT.RData")
save (FAP_damaged , file = "../saved in R /FAP_damaged.RData")
save (inflammatory_WT , file = "../saved in R /inflammatory_WT.RData")
save (muscleProgenitors_WT , file = "../saved in R /muscleProgenitors_WT.RData")
save (muscleProgenitors_damaged , file = "../saved in R /muscleProgenitors_damaged.RData")
