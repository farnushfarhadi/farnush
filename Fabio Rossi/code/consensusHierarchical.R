#### concatenate data 
logTransform (EC_WT , 3) -> EC_WT_log
logTransform(EC_damaged , 3) -> EC_KO_log
logTransform (FAP_WT , 3) -> FAP_WT_log
logTransform (FAP_damaged , 3) -> FAP_KO_log

EC_WT_log_filtered_n$tracking_id %>% as.character() -> EC_WT_genes
EC_damaged_log_filtered_n$tracking_id %>% as.character() -> EC_KO_genes
FAP_WT_log_filtered_n$tracking_id %>% as.character() -> FAP_WT_genes
FAP_damaged_log_filtered_n$tracking_id %>% as.character() -> FAP_KO_genes

intersect(EC_WT_genes , intersect(EC_KO_genes , intersect(FAP_WT_genes , FAP_KO_genes))) -> common_genes
getTableWithCommonGenes <- function(table , common_genes)
{
  table$tracking_id %>% as.character() -> genes
  sapply(common_genes, function(x){which (genes == x)}) -> idx
  return (table [idx %>% unlist(), ])
}

# 7156 genes
getTableWithCommonGenes (EC_WT_log , common_genes) -> EC_WT_log_c
getTableWithCommonGenes (EC_KO_log , common_genes) -> EC_KO_log_c
getTableWithCommonGenes (FAP_WT_log , common_genes) -> FAP_WT_log_c
getTableWithCommonGenes (FAP_KO_log , common_genes) -> FAP_KO_log_c

QuantileNormalize (EC_WT_log_c , 3) -> EC_WT_log_c_n 
QuantileNormalize (EC_KO_log_c , 3) -> EC_KO_log_c_n 
QuantileNormalize (FAP_WT_log_c , 3) -> FAP_WT_log_c_n 
QuantileNormalize (FAP_KO_log_c , 3) -> FAP_KO_log_c_n 

stdRow_table (EC_WT_log_c_n , 3) -> EC_WT_log_n_std
stdRow_table(EC_KO_log_c_n , 3) -> EC_KO_log_n_std
stdRow_table (FAP_WT_log_c_n , 3) -> FAP_WT_log_n_std
stdRow_table (FAP_KO_log_c_n , 3) -> FAP_KO_log_n_std

cbind (EC_WT_log_n_std[ , c(3,4,7,10,11,12,13,14,16)] , EC_KO_log_n_std[ , - c(1,2)],
FAP_WT_log_n_std[ , c(3,4,5,8,11,13,14,15,16,17)] , FAP_KO_log_n_std[ , c(3,4,5,6,8,9,10,11)] ) -> all_EC_FAP_rep1_log_std_filtered_n

all_EC_FAP_rep1_days <- c(paste("EC_WT",EC_WT_days_rep1,sep = "-") ,paste("EC_KO" , EC_KO_days_rep1,sep="-") ,
paste ("FAP_WT" ,FAP_WT_days_rep1, sep = "-") ,paste ("FAP_KO", FAP_KO_days_rep1,sep = "-") )
colnames(all_EC_FAP_rep1_log_std_filtered_n) <- all_EC_FAP_rep1_days
cbind(EC_WT_log_c[, c(1,2)] , all_EC_FAP_rep1_log_std_filtered_n  ) -> all_EC_FAP_rep1_log_std_filtered_n

# plotExpressionDist (allSamples , -0.5 , -0.3)
# threshold is -0.3 when we standardize first

#filterLowExpressedGenes(all_EC_FAP_rep1_log_std , -0.3 , 1) -> all_EC_FAP_rep1_log_std_filtered
hierarchical_clustering(gene_entry_pathways_genes[[i]] %>% as.character() , 
EC_WT_log_filtered_n_rep1 , 3 , paste (t , "EC_WT" , sep = "/") )  -> r1

downstreamGenes_consensus <- function (gene_entry , symbol, file)
{# for a given gene, it finds the kegg pathways that have the gene in 
  # and then cluster them hierarchialy for different cell types
  gene_entry_inf <- try(keggGet(gene_entry), silent=TRUE)
  if ('try-error' %in% class(gene_entry_inf)) {print ("NA")}
  if ('try-error' %in% class(gene_entry_inf)) {return ("NA")}
  gene_entry_inf[[1]]$PATHWAY %>% names() %>% unique()-> gene_entry_pathways
  if (gene_entry_pathways %>% length() == 0) {print ("NA")}
  if (gene_entry_pathways %>% length() == 0) {return ("NA")}
  sapply(gene_entry_pathways, function(x){getPathwayGenesEntrez(x)}) -> gene_entry_pathways_genes
  path = paste (file , paste0 (symbol , ".pdf") , sep = "/")
  pdf(file=path)
  
  res = c()
  for ( i in 1:length(gene_entry_pathways_genes))
  {
    print (paste (length(gene_entry_pathways_genes) , i , sep = "-"))
    (paste (gene_entry_pathways[i],symbol, sep = "/") -> t)
    hierarchical_clustering(gene_entry_pathways_genes[[i]] %>% as.character() , 
    all_EC_FAP_rep1_log_std_filtered_n , 3 , paste (t , "all" , sep = "/") )  -> r1
    res = c(res , r1)
  }
  dev.off()
  return (res)
}

file = "~/Farnush/farnush/Fabio Rossi/cell cell comunication/consensus clustering/average 60/"
#length(kegg_entries_receptors)
for (i in 1:10)
{# do the downstream analysis for alll receptors 
  print (paste ("******", i , sep = ":  "))
  print(kegg_entries_receptors[i])
  print(receptor_symbol[i])
  downstreamGenes_consensus(kegg_entries_receptors[i] , receptor_symbol[i] , file)
}
