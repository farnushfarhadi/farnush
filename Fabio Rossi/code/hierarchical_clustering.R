hierarchical_clustering <- function (genes , table , jj  , title)
{
  #genes <- genes %>% unlist() %>% unname()
  genes_2 <- 0
  ids <- toupper (table$tracking_id)
  (sapply(genes, function(x) {which ( ids== x) })   -> cluster_idx  )
  if  (length(unlist (cluster_idx)) == 0) {
    print ("there is no available gene in the table")
    print(title)
    return("NA")}else if (length(unlist (cluster_idx)) > 2){
    table_genes <- table [ Reduce(union , cluster_idx), -c(1:(jj-1))] %>% as.matrix() 
    genes_2 <- 1
    #constructing dist matrix
    cor (t(as.matrix(table_genes[ , -c(1:(jj-1))]) ) , t(as.matrix(table_genes[ , -c(1:(jj-1))])) ) -> dist
    dist_m = 1 - dist
    #dist_m = 1 - abs (dist)
    colnames(dist_m) = rownames(dist_m) = table$tracking_id[Reduce(union , cluster_idx)] %>% as.character()
    #hclust(as.dist( dist_m ) , method  = "complete" ) -> res_hclust
    hclust(as.dist( dist_m ) , method  = "average" ) -> res_hclust
      par(cex=0.4, mar=c(5, 8, 4, 1))
      plot(res_hclust , xlab = "" , main = title , hang = -1)
      abline(h = 0.25, lty = 1 , col = "red")
      cutree(res_hclust , h = 0.25) %>% unlist() %>% unique() %>% length() -> k
      # 1.2 , 2 - 1.1 , 1.9 - 1 , 1.8
      text(1.4, 1.2 , paste ("clusters:" ,k, sep = " ") , col = "red")
      abline(h = 0.4, lty = 1 , col = "blue")
      cutree(res_hclust , h = 0.4) %>% unlist() %>% unique() %>% length() -> k_
      text(1.4, 1.1 , paste ("clusters:" ,k_, sep = " ") , col = "blue")
      text(1.4, 1 , paste ("number of genes:" ,dim(dist_m)[1], sep = " ") , col = "mediumorchid4")
      num = dim(dist_m)[1]
    }
    if (genes_2 == 0)
    {
      num = 2;
      k_ = 1;
    }
      
    ret = c (k_ , num)
    names (ret) = c ("clusters" , "number of genes")
    # when plotting return cutree O.W. return ret
    return (ret)
    #return (cutree(res_hclust , h = 0.4) )
}


# cutree(tree, k = NULL, h = NULL)

#testing 
i = 4
strsplit (t [i , 3] %>% as.character() ,";" ) %>% unlist() %>% toupper() -> cluster
(sapply(cluster, function(x) {which ( ids== x) })   -> cluster_idx  )
table_genes <- table [ Reduce(union , cluster_idx), -c(1:(jj-1))] %>% as.matrix() 

#constructing dist matrix
cor (t(as.matrix(table_genes[ , -c(1:(jj-1))]) ) , t(as.matrix(table_genes[ , -c(1:(jj-1))])) ) -> dist
dist_m = 1 - dist
dist_m = 1 - abs (dist)
colnames(dist_m) = rownames(dist_m) = table$tracking_id[Reduce(union , cluster_idx)] %>% as.character()

hclust(as.dist( dist_m ) , method  = "complete" ) -> res_hclust
par(cex=0.4, mar=c(5, 8, 4, 1))
plot(res_hclust , xlab = "")
cutree(res_hclust , h = 0.4)
cutree(res_hclust , k = 4)

plotGroupsOfHierarchicalclustering <- function (table , symbol , pathway_idx ,dayNames, jj ,cellType)
{
  which (receptor_symbol == symbol) -> idx
  (gene_entry = kegg_entries_receptors [idx])
  (symbol = receptor_symbol[idx])
  gene_entry_inf <- try(keggGet(gene_entry), silent=TRUE)
  gene_entry_inf[[1]]$PATHWAY %>% names() %>% unique()-> gene_entry_pathways
  
  getPathwayGenesEntrez (gene_entry_pathways [ pathway_idx]) %>% unlist() %>% as.character()-> genes
  hierarchical_clustering(genes , table , jj , "") -> r
  r %>% unlist() %>% unique() %>% length() -> k
  plist = list ()
  for ( i in 1:k)
  {
    which (r == i) %>% names() -> g
    plotCluster(table , g %>% toupper(), paste(i , cellType,sep = "/") , dayNames,jj) -> p
    plist [[i]] <- p
  }
  do.call("grid.arrange", c(plist, ncol=3))
  
}

# common clusters 38  67 172 231 286 294
file = "~/Farnush/farnush/Fabio Rossi/cell cell comunication/receptor_keg_hclust_average_60/common_TP_small/"
file = "~/Farnush/farnush/Fabio Rossi/cell cell comunication/receptor_keg_hclust_complete_60/common_TP_small/"
file = "~/Farnush/farnush/Fabio Rossi/cell cell comunication/receptor_keg_hclust_complete_60/common_TP_small_original/"

for (i in 1:length(common))
{
  print (i)
  (kk = common[i])
  
  (s = symbols[kk] %>% as.character())
  which (receptor_symbol == s) -> i 
  (receptor_symbol[i] -> symbol)
  print(symbol)
  pdf (file = paste0 (paste0 (file , symbol) , ".pdf" ) )
  (idx1 = which (ec_wt_small$receptor_idx == kk))
  (idx2 = which (ec_ko_small$receptor_idx == kk))
  (idx3 = which (fap_wt_small$receptor_idx == kk))
  (idx4 = which (fap_ko_small$receptor_idx == kk))
  plotGroupsOfHierarchicalclustering (table = EC_WT_log_filtered_n_rep1 , symbol , 
  ec_wt_small$pathway_idx[idx1] , EC_WT_days_rep1,3 , "EC_WT")
  plotGroupsOfHierarchicalclustering (table = EC_damaged_log_filtered_n , symbol , 
  ec_ko_small$pathway_idx[idx2] , EC_damaged_days,3 , "EC_KO")
  plotGroupsOfHierarchicalclustering (table = FAP_WT_log_filtered_n_rep1 , symbol , 
  fap_wt_small$pathway_idx[idx3] , FAP_WT_days_rep1,3 , "FAP_WT")
  plotGroupsOfHierarchicalclustering (table = FAP_KO_log_filtered_n_rep1 , symbol , 
  fap_ko_small$pathway_idx[idx4] , FAP_KO_days_rep1,3 , "FAP_KO")
  dev.off()  
  
}

plotTemporalPatterns <- function (table, jj , dayNames , smallClusterInf , file)
{ # get the result of small clusters for a specific cell type
  # plot the temporal patterns of the groups of the hierarchical clustering
  for ( i in 1:length(smallClusterInf$receptor_idx ))
  {
    idx <- smallClusterInf$receptor_idx
    p_idx <- smallClusterInf$pathway_idx
    symbols[idx[i]] %>% as.character() -> s
    path = paste (file , paste0 (s , ".pdf") , sep = "/")
    pdf(file=path)
    plotGroupsOfHierarchicalclustering (table = table , s , p_idx[i] , dayNames,jj)
    dev.off()
  }  
}

file = "~/Documents/Farnush GitHub/Fabio Rossi/cell cell comunication/fap-wt-smallclusers-TP/"
file = "~/Documents/Farnush GitHub/Fabio Rossi/cell cell comunication/fap-ko-smallclusers-TP/"
plotTemporalPatterns (FAP_WT_log_filtered_n_rep1 , 3, FAP_WT_days_rep1 , fap_wt_small, file)
plotTemporalPatterns (FAP_KO_log_filtered_n_rep1 , 3, FAP_KO_days_rep1 , fap_ko_small, file)


