#source("https://bioconductor.org/biocLite.R")
#biocLite("KEGGREST")




# converting gene ID to entrez ID 
library (org.Hs.eg.db)
genes_symbol_entrez <- as.list(org.Hs.egALIAS2EG)
genes_symbol_entrez_n <- names (genes_symbol_entrez)
# filtering 
sapply(genes_symbol_entrez , function(x){x %>% length() > 1 }) -> moreThan1  # 118K
genes_symbol_entrez [! moreThan1] -> genes_symbol_entrez_unique # 114K

##### Process data for using kegg
#data.frame (symbol = names(genes_symbol_entrez_unique) %>% toupper(), id = genes_symbol_entrez_unique %>% unname()) -> database
cbind (names(genes_symbol_entrez_unique) %>% toupper() ,genes_symbol_entrez_unique %>% unname() ) -> database
colnames(database) = c("symbol" , "id")
genes_mydata = c()
genes_mydata = cbind(genes_mydata , EC_WT$tracking_id %>% as.character() %>% toupper())
colnames(genes_mydata) = "symbol"
#data.frame( symbol = EC_WT$tracking_id %>% as.character() %>% toupper())-> genes_mydata

merge (genes_mydata , database , by.x = "symbol" , by.y = "symbol" ) -> merged # 15622 genes out of 23516 genes

mergeCellWithID <- function (table , merged , jj)
{
  table$tracking_id = toupper(table$tracking_id)
  merge(x = table ,y= merged , by.x = "tracking_id" , by.y = "symbol" ) -> table_id
  table_id = table_id [ ,  c(1:(jj-1) , dim(table_id)[2] , jj:(dim(table_id)[2]-1))]
  return (table_id)
}
mergeCellWithID (EC_WT , merged , 3) -> EC_WT_id # 15688 genes for all next files as well
mergeCellWithID (EC_damaged , merged , 3) -> EC_damaged_id
mergeCellWithID (FAP_WT , merged , 3) -> FAP_WT_id
mergeCellWithID (FAP_damaged , merged ,3) -> FAP_damaged_id
mergeCellWithID (inflammatory_WT , merged ,2) -> inflammatory_WT_id
mergeCellWithID (muscleProgenitors_WT , merged , 2) -> muscleProgenitors_WT_id
mergeCellWithID (muscleProgenitors_damaged, merged , 2) -> muscleProgenitors_damaged_id


EC_WT_days <- c("D0" , "D2-1" , "D2-2" , "D2-3" , "D3-1" , "D3-2" , "D3-3" , "D4" , "D5" , "D6" , "D7" , "D10-1" , "D10-2" , "D14" )
dataQC (EC_WT_id , 0.7 , 1 , "EC WT-quantile normalized" , EC_WT_days , 4) -> EC_WT_id_log_filtered_n # 8230 16 genes
EC_damaged_days <- c("D0" , "D1" , "D2" , "D3" , "D5" , "D6" , "D7" , "D10")
dataQC (EC_damaged_id , 0.7 , 1 , "EC KO-quantile normalized" , EC_damaged_days , 4) -> EC_damaged_id_log_filtered_n # 7921 10 genes
FAP_WT_days <- c("D0" , "D1" , "D2-1" , "D2-2" , "D2-3" , "D3-1" , "D3-2" , "D3-3" , "D4-1" , "D4-2" , "D5" , "D6" , "D7" , "D10" , "D14" )
dataQC (FAP_WT_id , 0.7 , 1 , "FAP WT-quantile normalized" , FAP_WT_days , 4) -> FAP_WT_id_log_filtered_n # 7889 17 genes
FAP_damaged_days <- c("D0" , "D1" , "D2" , "D3-1" , "D3-2" , "D4" , "D5" , "D6"  , "D10")
dataQC (FAP_damaged_id , 0.7 , 1 , "FAP KO-quantile normalized" , FAP_damaged_days ,4) -> FAP_damaged_id_log_filtered_n # 7405 10 genes
muscleProgenitors_WT_days <- c( "D1" , "D2" , "D3-1" , "D3-2" , "D5" ,  "D7"  , "D10")
dataQC (muscleProgenitors_WT_id , 0.7 , 1 , "muscleProgenitors WT-quantile normalized", muscleProgenitors_WT_days,3) -> muscleProgenitors_WT_id_log_filtered_n # 7826 8 genes
muscleProgenitors_damaged_days <- c( "D0" , "D3-1" , "D3-2" ,"D4" , "D5" ,"D6"  , "D10")
dataQC (muscleProgenitors_damaged_id , 0.7 , 1 , "muscleProgenitors KO-quantile normalized" ,muscleProgenitors_damaged_days, 3) -> muscleProgenitors_damaged_id_log_filtered_n # 8151 8 genes
inflammatory_WT_days <- c( "D1-1", "D1-2" ,"D2-1" , "D2-2" , "D2-3" , "D3-1" , "D3-2" , "D3-3" , "D4"  ,"D5" ,"D6","D7-1" ,"D7-2" , "D10")
dataQC (inflammatory_WT_id , 0.7 , 1 , "inflammatory WT-quantile normalized", inflammatory_WT_days,3) -> inflammatory_WT_id_log_filtered_n # 2002 15 genes


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


### loading receptors and ligands
setwd("~/Documents/Farnush github/Fabio Rossi/cell cell comunication/")
read.delim("receptor-id.txt" , header = TRUE) -> receptors
#receptors <- receptors %>% unlist() %>% as.character()
read.delim("UniqueLigands.txt" , header = TRUE) -> ligands
ligands <- ligands %>% unlist() %>% as.character()

# generating receptors entrez id so we can get information about them from kegg
receptors$Gene.Symbol %>% as.character() -> receptor_symbol
receptors$Gene.ID %>% as.character() -> receptor_id
#sapply(receptors, function(x){which (genes_symbol_entrez_n == toupper(x))}) -> receptors_entrez_idx
sapply(receptor_id, function(x){return (paste0("mmu:",x %>% unname()))} ) %>% unname() -> kegg_entries_receptors

i= 1
(gene_entry = kegg_entries_receptors [i])
(symbol = receptor_symbol[i])
#length(kegg_entries_receptors
receptor_res = c()
#length(kegg_entries_receptors)
for (i in 1:length(kegg_entries_receptors))
{
  print(kegg_entries_receptors[i])
  print(receptor_symbol[i])
  downstreamGenes(kegg_entries_receptors[i] , receptor_symbol[i] ) -> res
  receptor_res = c (receptor_res , res)
}
receptor_res_average <- receptor_res
receptor_res_complete <- receptor_res

setwd ("~/Documents/Farnush GitHub/Fabio Rossi/cell cell comunication/")
save (receptor_res_complete , file = "receptor_keg_hclust/hclust_complete.RData")
save (receptor_res_average , file = "receptor_keg_hclust_average/hclust_average.RData")
load("receptor_keg_hclust/hclust_complete.RData")


receptor_res_average -> cluster_res
# bug : can not partiton the res to different cell types
hierarchicalClusteringAnalysis <- function (cluster_res)
{
  
  (which (cluster_res == "NA") %>% unname() -> idx)
  # for average , the idx 1601  is symbol  so we should remove it
  #idx  = c (idx , 1601)
  # for complete
  idx = c (idx , 1391)
  cluster_res [ -idx] -> cluster_res
  seq (from = 1 , to =length(cluster_res)  , by = 5) -> ec_wt_idx
  seq (from = 2 , to =length(cluster_res)  , by = 5) -> ec_ko_idx
  seq (from = 3 , to =length(cluster_res)  , by = 5) -> fap_wt_idx
  seq (from = 4 , to =length(cluster_res)  , by = 5) -> fap_ko_idx
  seq (from = 5 , to =length(cluster_res)  , by = 5) -> symbols_idx
  cluster_res [ ec_wt_idx] -> ec_wt_res
  cluster_res [ec_ko_idx] -> ec_ko_res
  cluster_res [fap_wt_idx] -> fap_wt_res
  cluster_res [fap_ko_idx] -> fap_ko_res
  cluster_res [symbols_idx] -> symbols
  plotDistK_genes (ec_wt_res , "EC wild type")
  plotDistK_genes (ec_ko_res , "EC kocked out")
  plotDistK_genes (fap_wt_res , "FAP wild type")
  plotDistK_genes (fap_ko_res , "FAP knocked out")
  (findSmallClusters (ec_wt_res , 10) -> ec_wt_small)
  (findSmallClusters (ec_ko_res , 10) -> ec_ko_small )
  (findSmallClusters (fap_wt_res , 10) -> fap_wt_small )
  (findSmallClusters (fap_ko_res , 10) -> fap_ko_small)
  intersect( ec_ko_small$receptor_idx %>% unlist() , intersect( fap_wt_small$receptor_idx %>% unlist(),
  fap_ko_small$receptor_idx %>% unlist() ) )
  #38  67 172 231 286 294
}

plotDistK_genes <- function(cell_res , title)
{
  cell_res %>% unlist() -> cell_res
  (which (cell_res == "NA") -> idx)
  if (length(idx) > 0){
    cell_res [ -idx] -> cell_res
  }
  #cell_res %>% unlist() -> cell_res
  cell_res[seq(1,length(cell_res), by= 2)] %>% unname() %>% as.numeric() -> cluster_size
  cell_res[seq(2,length(cell_res), by= 2)] %>% unname() %>% as.numeric() -> gene_size
  par(mfrow = c(1,2))
  hist (cluster_size , breaks = 25 , col = "turquoise" , main = title)
  hist (gene_size , breaks = 25 , col = "lightseagreen" , main = title)
  return (mean (cluster_size))
  
}

findSmallClusters <- function(cell_res , t)
{
  # searcgh for symbols that have small size clusters in any pathway in any of 4 cell types
  #( using ec_wt_res , ec_ko_res , fap_wt_res , fap_ko_res)
  # visualize their clusters in pdf files
#   cell_res %>% unlist() -> cell_res
#   which (cell_res == "NA") -> idx
#   if (length(idx) > 0){
#     cell_res [ -idx] -> cell_res
#   }
  idx = c() # storing which receptpr has small cluster
  idx_inner = c() # storing which pathway of the receptor has the small cluster
  for (i in 1:length(cell_res))
  {
    cell_res[[i]] [seq (1 , length(cluster_res[[i]]) , by = 2)] %>% unname() -> cluster_sizes
    which (cluster_sizes < t) -> cluster_small_idx
    if (length(cluster_small_idx) > 0)
    {
      idx = c(idx , i)
      idx_inner = c (idx_inner , cluster_small_idx[1])
    }
  }
  ret = list (idx , idx_inner)
  names (ret) <- c ("receptor_idx" , "pathway_idx")
  return(ret)
}
downstreamGenes <- function (gene_entry , symbol)
{
  gene_entry_inf <- try(keggGet(gene_entry), silent=TRUE)
  if ('try-error' %in% class(gene_entry_inf)) {print ("NA")}
  if ('try-error' %in% class(gene_entry_inf)) {return ("NA")}
  gene_entry_inf[[1]]$PATHWAY %>% names() %>% unique()-> gene_entry_pathways
  if (gene_entry_pathways %>% length() == 0) {print ("NA")}
  if (gene_entry_pathways %>% length() == 0) {return ("NA")}
  sapply(gene_entry_pathways, function(x){getPathwayGenesEntrez(x)}) -> gene_entry_pathways_genes
  #if (gene_entry_pathways %>% length() == 1)
   # gene_entry_pathways_genes = list (gene_entry_pathways_genes)
  
  # for each pathway, have a plot
  #file = "~/Documents/Farnush github/Fabio Rossi/cell cell comunication/receptor_keg_plots/"
  file = "~/Documents/Farnush github/Fabio Rossi/cell cell comunication/receptor_keg_hclust/"
  #file = "~/Documents/Farnush github/Fabio Rossi/cell cell comunication/receptor_keg_hclust_average/"
  path = paste (file , paste0 (symbol , ".pdf") , sep = "/")
  pdf(file=path)
  EC_wt_res = c ()
  EC_ko_res = c ()
  FAP_wt_res = c ()
  FAP_ko_res = c ()
  for ( i in 1:length(gene_entry_pathways_genes))
  {
    print (paste (length(gene_entry_pathways_genes) , i , sep = "-"))
    paste (gene_entry_pathways[i],symbol, sep = "/") -> t
    hierarchical_clustering(gene_entry_pathways_genes[[i]] %>% as.character() , 
    EC_WT_log_filtered_n_rep1 , 3 , paste (t , "EC_WT" , sep = "/") )  -> r1
    EC_wt_res = c(EC_wt_res , r1)
    hierarchical_clustering(gene_entry_pathways_genes[[i]] %>% as.character() , 
    EC_KO_log_filtered_n_rep1 , 3 , paste (t , "EC_KO" , sep = "/")) -> r2
    EC_ko_res = c(EC_ko_res , r2)
    hierarchical_clustering(gene_entry_pathways_genes[[i]] %>% as.character() ,
    FAP_WT_log_filtered_n_rep1 , 3 , paste (t , "FAP_WT" , sep = "/")) -> r3
    FAP_wt_res = c( FAP_wt_res , r3)
    hierarchical_clustering(gene_entry_pathways_genes[[i]] %>% as.character() , 
    FAP_KO_log_filtered_n_rep1 , 3 , paste (t , "EC_KO" , sep = "/")) -> r4
    FAP_ko_res = c(FAP_ko_res , r4)
#     plotCluster(EC_WT_log_filtered_n_rep1 , gene_entry_pathways_genes[[i]] %>% toupper(),paste (t , "EC_WT" , sep = "/") , EC_WT_days_rep1,3) -> p1
#     plotCluster(EC_KO_log_filtered_n_rep1 , gene_entry_pathways_genes[[i]]%>% toupper(),paste (t , "EC_KO" , sep = "/"), EC_KO_days_rep1,3) -> p2
#     plotCluster(FAP_WT_log_filtered_n_rep1 , gene_entry_pathways_genes[[i]]%>% toupper(),paste (t , "FAP_WT" , sep = "/"), FAP_WT_days_rep1,3) -> p3
#     plotCluster(FAP_KO_log_filtered_n_rep1 , gene_entry_pathways_genes[[i]]%>% toupper(),paste (t , "EC_KO" , sep = "/"), FAP_KO_days_rep1,3) -> p4
#     multiplot(p1, p2 , p3 , p4 , cols = 1)
    #if (r1 == "NA" || r2 == "NA" || r3 == "NA" || r4 == "NA")
     # contin
#     plot(r1)
#     plot(r2)
#     plot(r3)
#     plot(r4)
  }
  dev.off()
  list (EC_wt_res , EC_ko_res , FAP_wt_res , FAP_ko_res , symbol) -> ret
  names (ret) <- c ("EC_WT" , "EC_KO" , "FAP_WT" , "FAP_KO" , "receptor symbol")
  return (ret)
  #return (gene_entry_pathways_genes)
}

getPathwayGenesEntrez <- function (pathway)
{
  keggGet(pathway) -> pathway_inf
  if (! is.null(pathway_inf[[1]]$GENE))
  {
    ii <- seq (from= 2 , to = length(pathway_inf[[1]]$GENE) , by = 2)
    strs <- pathway_inf[[1]]$GENE[ii]
    sapply(strs, function(x){ return (substr ( x , start = 1 ,stop = unlist(gregexpr (";" , x ))-1 ) %>% toupper())}) %>% unname() -> res
    return (as.data.frame (res))
  } else {
    return ("NA")
  }
}

for (i in 1: length(kegg_entries))
{
  # each kegg entry can have multiple pathways (we do not know which is the one with the gene as the receptor!)
  downstreamGenes (kegg_entries[i]) -> entry_pathway_genes
  
  
}

###### Paolo KEGG 
setwd ("~/Documents/Farnush github/Fabio Rossi/paolo/Kegg Time Clip results/")
# EC-wt
PlotCellSpecificKeggPathways <- function (file_name)
{
  read.delim(paste (file_name,"paths2genes.txt" , sep = "-")) -> table
  path = "~/Documents/Farnush github/Fabio Rossi/paolo/"
  pdf(file=paste (path, paste (file_name , "specific_pathways.pdf" , sep = "-") , sep = "/"))
  for ( i in 1: dim(table)[1])
  {
    strsplit (table [i , 3] %>% as.character() ,";" ) %>% unlist() %>% toupper() -> cluster
    
    plotCluster(EC_WT_log_filtered_n_rep1 , cluster , "EC_WT"  , EC_WT_days_rep1,3) -> p1
    plotCluster(EC_KO_log_filtered_n_rep1 , cluster , "EC_KO" , EC_KO_days_rep1,3) -> p2
    plotCluster(FAP_WT_log_filtered_n_rep1 , cluster,"FAP_WT" , FAP_WT_days_rep1,3) -> p3
    plotCluster(FAP_KO_log_filtered_n_rep1 , cluster,"FAP_KO", FAP_KO_days_rep1,3) -> p4
    
    grid.arrange(p1,p2,p3,p4, top = table[i,1]%>% as.character(),
                 layout_matrix = matrix(c(1,2,3,4), ncol=1, byrow=TRUE))
    #multiplot(p1, p2 , p3 , p4 , cols = 1)
  }
  dev.off()
}

PlotCellSpecificKeggPathways ("EC-wt")
PlotCellSpecificKeggPathways ("EC-ko")
PlotCellSpecificKeggPathways ("FAP-wt")
PlotCellSpecificKeggPathways ("FAP-ko")








