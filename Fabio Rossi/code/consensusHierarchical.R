
# construcing replicates
EC_WT_rep1 = EC_WT[ , c(1,2,3,4,7,10,11,12,13,14,16)]
EC_WT_days_rep1 <- c("D0" , "D2" , "D3"  , "D4" , "D5" , "D6" , "D7" , "D10"  , "D14" )
EC_KO_rep1 = EC_damaged
EC_KO_days_rep1 <- c("D0" , "D1" , "D2" , "D3" , "D5" , "D6" , "D7" , "D10")
FAP_WT_rep1 = FAP_WT[ , c(1,2,3,4,5,8,11,13,14,15,16,17)]
FAP_WT_days_rep1 <- c("D0" ,"D1" , "D2" , "D3"  , "D4" , "D5" , "D6" , "D7" , "D10"  , "D14" )
MP_WT_rep1 = muscleProgenitors_WT[ , c(1,2,3,4,6,7,8)]
MP_WT_days_rep1 <- c("D1" , "D2" , "D3"  , "D5" , "D7" , "D10"  )
MP_KO_rep1 = muscleProgenitors_damaged[ , c(1,2,3,5,6,7,8)]
MP_KO_days_rep1 <- c("D0" , "D3" , "D4" , "D5" , "D6"  , "D10")
IC_WT_rep1 = inflammatory_WT [ , c(1,2,4,7,10,11,12,13,15)]
IC_WT_days_rep1 <- c("D1" , "D2" , "D3"  , "D4" , "D5" , "D6" , "D7" , "D10"  )
# low exp and normalizing
dataQC (EC_WT_rep1 , 0.7 , 0.5 , "" , EC_WT_days_rep1 , 3) -> EC_WT_rep1_l_f_n #11296 11 genes
dataQC (EC_KO_rep1 , 0.7 , 0.5 , "" , EC_KO_days_rep1 , 3) -> EC_KO_rep1_l_f_n #11241 10 genes
dataQC (FAP_WT_rep1 , 0.7 , 0.5 , "" , FAP_WT_days_rep1 , 3) -> FAP_WT_rep1_l_f_n #11348 12 genes
dataQC (FAP_KO_rep1 , 0.7 , 0.5 , "" , FAP_KO_days_rep1 , 3) -> FAP_KO_rep1_l_f_n #11084 10 genes
dataQC (MP_WT_rep1 , 0.7 , 0.5 , "" , MP_WT_days_rep1 , 2) -> MP_WT_rep1_l_f_n #11781 7 genes
dataQC (MP_KO_rep1 , 0.7 , 0.5 , "" , MP_KO_days_rep1 , 2) -> MP_KO_rep1_l_f_n #11558 7 genes
dataQC( IC_WT_rep1 , 0.7 , 0.5 , "" , IC_WT_days_rep1 , 2) -> IC_WT_rep1_l_f_n # 10489 9 genes
# building cells with common genes

EC_WT_rep1_l_f_n$tracking_id %>% as.character() %>% toupper() -> EC_WT_genes
EC_KO_rep1_l_f_n$tracking_id %>% as.character() %>% toupper() -> EC_KO_genes
FAP_WT_rep1_l_f_n$tracking_id %>% as.character() %>% toupper() -> FAP_WT_genes
FAP_KO_rep1_l_f_n$tracking_id %>% as.character() %>% toupper() -> FAP_KO_genes
MP_WT_rep1_l_f_n$tracking_id %>% as.character() %>% toupper() -> MP_WT_genes
MP_KO_rep1_l_f_n$tracking_id %>% as.character() %>% toupper() -> MP_KO_genes
IC_WT_rep1_l_f_n$tracking_id %>% as.character() %>% toupper() -> IC_WT_genes


intersect(EC_WT_genes , intersect(EC_KO_genes , 
         intersect(FAP_WT_genes , 
         intersect (FAP_KO_genes , 
         intersect (MP_WT_genes , 
         intersect (MP_KO_genes , IC_WT_genes)))) )) -> common_genes
# 9961 genes for FAP and EC
# 8988 for all 

getTableWithCommonGenes <- function(table , common_genes)
{
  table$tracking_id %>% as.character() %>% toupper() -> genes
  sapply(common_genes, function(x){which (genes == x)}) -> idx
#   idx = c()
#   for (i in 1:length(common_genes))
#   {
#     which(genes)
#   }
  #return (table [Reduce(union , idx), ])
  return(Reduce(union , idx))
}

# 9998 genes
# 8988 for all 
# getTableWithCommonGenes (EC_WT_rep1_l_f_n , common_genes) -> EC_WT_rep1_l_f_n_c
# getTableWithCommonGenes (EC_KO_rep1_l_f_n , common_genes) -> EC_KO_rep1_l_f_n_c
# getTableWithCommonGenes (FAP_WT_rep1_l_f_n , common_genes) -> FAP_WT_rep1_l_f_n_c
# getTableWithCommonGenes (FAP_KO_rep1_l_f_n , common_genes) -> FAP_KO_rep1_l_f_n_c
# getTableWithCommonGenes (MP_WT_rep1_l_f_n , common_genes) -> MP_WT_rep1_l_f_n_c
# getTableWithCommonGenes (MP_KO_rep1_l_f_n , common_genes) -> MP_KO_rep1_l_f_n_c
getTableWithCommonGenes (IC_WT_rep1_l_f_n , common_genes) -> IC_WT_rep1_l_f_n_c_idx
# since they change in size , I use the IC (some ha s 9018 and some has 9019)
EC_WT_rep1_l_f_n[IC_WT_rep1_l_f_n_c_idx , ] -> EC_WT_rep1_l_f_n_c
EC_KO_rep1_l_f_n[IC_WT_rep1_l_f_n_c_idx , ] -> EC_KO_rep1_l_f_n_c
FAP_WT_rep1_l_f_n [ IC_WT_rep1_l_f_n_c_idx, ] -> FAP_WT_rep1_l_f_n_c
FAP_KO_rep1_l_f_n [ IC_WT_rep1_l_f_n_c_idx , ] -> FAP_KO_rep1_l_f_n_c
MP_WT_rep1_l_f_n [ IC_WT_rep1_l_f_n_c_idx , ] -> MP_WT_rep1_l_f_n_c
MP_KO_rep1_l_f_n [ IC_WT_rep1_l_f_n_c_idx , ] -> MP_KO_rep1_l_f_n_c
IC_WT_rep1_l_f_n [ IC_WT_rep1_l_f_n_c_idx , ] -> IC_WT_rep1_l_f_n_c
# 9018 after all (for genes copy)  (8988)


first_two <- c ("tracking_id" , "locus")
colnames(EC_WT_rep1_l_f_n_c) = c (first_two , EC_WT_days_rep1)
colnames(EC_KO_rep1_l_f_n_c) = c (first_two , EC_KO_days_rep1)
colnames(FAP_WT_rep1_l_f_n_c) = c (first_two , FAP_WT_days_rep1)
colnames (FAP_KO_rep1_l_f_n_c) = c (first_two , FAP_KO_days_rep1)
colnames(MP_WT_rep1_l_f_n_c) = c(first_two[1] , MP_WT_days_rep1)
colnames(MP_KO_rep1_l_f_n_c) = c(first_two[1] , MP_KO_days_rep1)
colnames(IC_WT_rep1_l_f_n_c) = c(first_two[1] , IC_WT_days_rep1)


# QuantileNormalize (EC_WT_log_c , 3) -> EC_WT_log_c_n 
# QuantileNormalize (EC_KO_log_c , 3) -> EC_KO_log_c_n 
# QuantileNormalize (FAP_WT_log_c , 3) -> FAP_WT_log_c_n 
# QuantileNormalize (FAP_KO_log_c , 3) -> FAP_KO_log_c_n 
# 
# stdRow_table (EC_WT_log_c_n , 3) -> EC_WT_log_n_std
# stdRow_table(EC_KO_log_c_n , 3) -> EC_KO_log_n_std
# stdRow_table (FAP_WT_log_c_n , 3) -> FAP_WT_log_n_std
# stdRow_table (FAP_KO_log_c_n , 3) -> FAP_KO_log_n_std
# 
# cbind (EC_WT_log_n_std[ , c(3,4,7,10,11,12,13,14,16)] , EC_KO_log_n_std[ , - c(1,2)],
# FAP_WT_log_n_std[ , c(3,4,5,8,11,13,14,15,16,17)] , FAP_KO_log_n_std[ , c(3,4,5,6,8,9,10,11)] ) -> all_EC_FAP_rep1_log_std_filtered_n
# 
# all_EC_FAP_rep1_days <- c(paste("EC_WT",EC_WT_days_rep1,sep = "-") ,paste("EC_KO" , EC_KO_days_rep1,sep="-") ,
# paste ("FAP_WT" ,FAP_WT_days_rep1, sep = "-") ,paste ("FAP_KO", FAP_KO_days_rep1,sep = "-") )
# colnames(all_EC_FAP_rep1_log_std_filtered_n) <- all_EC_FAP_rep1_days
# cbind(EC_WT_log_c[, c(1,2)] , all_EC_FAP_rep1_log_std_filtered_n  ) -> all_EC_FAP_rep1_log_std_filtered_n
# 

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
  #pdf(file=path)
  
  res = list()
  for ( i in 1:length(gene_entry_pathways_genes))
  {
    print (paste (length(gene_entry_pathways_genes) , i , sep = "-"))
    (paste (gene_entry_pathways[i],symbol, sep = "/") -> t)
    consensus_hclust(gene_entry_pathways_genes[[i]] %>% as.character() , 1)  -> r1
    res[[i]] <- r1
  }
  #dev.off()
  names(res) = gene_entry_pathways
  return (res)
}


consensus_hclust <- function (genes , avg)
{
  genes_2 <- 0
  ids <- EC_WT_rep1_l_f_n_c$tracking_id %>% as.character() %>% toupper()
  (sapply(genes, function(x) {which ( ids== x) })   -> cluster_idx  )
  if  (length(unlist (cluster_idx)) == 0) {
    print ("there is no available gene in the table")
    print(title)
    return("NA")}else if (length(unlist (cluster_idx)) > 2){
      table_genes_ec_wt <- EC_WT_rep1_l_f_n_c [ Reduce(union , cluster_idx), -c(1,2)] %>% as.matrix() 
      table_genes_ec_ko <- EC_KO_rep1_l_f_n_c [ Reduce(union , cluster_idx), -c(1,2)] %>% as.matrix() 
      table_genes_fap_wt <- FAP_WT_rep1_l_f_n_c [ Reduce(union , cluster_idx), -c(1,2)] %>% as.matrix() 
      table_genes_fap_ko <- FAP_KO_rep1_l_f_n_c [ Reduce(union , cluster_idx), -c(1,2)] %>% as.matrix() 
      
      genes_2 <- 1
      #constructing dist matrix
      cor (t(as.matrix(table_genes_ec_wt[ , -c(1,2)]) ) , 
      t(as.matrix(table_genes_ec_wt[ , -c(1,2)])) ) -> dist_ec_wt
      
      cor (t(as.matrix(table_genes_ec_ko[ , -c(1,2)]) ) , 
      t(as.matrix(table_genes_ec_ko[ , -c(1,2)])) ) -> dist_ec_ko
      
      cor (t(as.matrix(table_genes_fap_wt[ , -c(1,2)]) ) , 
      t(as.matrix(table_genes_fap_wt[ , -c(1,2)])) ) -> dist_fap_wt
      
      cor (t(as.matrix(table_genes_fap_ko[ , -c(1,2)]) ) , 
      t(as.matrix(table_genes_fap_ko[ , -c(1,2)])) ) -> dist_fap_ko
      dist <- (dist_ec_wt + dist_ec_ko + dist_fap_wt + dist_fap_ko)/4
      dist_m = 1 - dist
      #dist_m = 1 - abs (dist)
      colnames(dist_m) = rownames(dist_m) = EC_WT_rep1_l_f_n_c$tracking_id[Reduce(union , cluster_idx)] %>% as.character()
      if (avg)
      {
        hclust(as.dist( dist_m ) , method  = "average" ) -> res_hclust  
      }else{
        hclust(as.dist( dist_m ) , method  = "complete" ) -> res_hclust 
      }
      
      # par(cex=0.4, mar=c(5, 8, 4, 1))
      # plot(res_hclust , xlab = "" , main = "EC_FAP" , hang = -1)
      # abline(h = 0.25, lty = 1 , col = "red")
      # cutree(res_hclust , h = 0.25) %>% unlist() %>% unique() %>% length() -> k
      # # 1.2 , 2 - 1.1 , 1.9 - 1 , 1.8
      # text(1.4, 1.2 , paste ("clusters:" ,k, sep = " ") , col = "red")
      # abline(h = 0.4, lty = 1 , col = "blue")
       cutree(res_hclust , h = 0.4) %>% unlist() %>% unique() %>% length() -> k_
      # text(1.4, 1.1 , paste ("clusters:" ,k_, sep = " ") , col = "blue")
      # text(1.4, 1 , paste ("number of genes:" ,dim(dist_m)[1], sep = " ") , col = "mediumorchid4")
       num = dim(dist_m)[1]
    }
  if (genes_2 == 0)
  {
    num = 2;
    k_ = 1;
    #return(list (colnames(dist_m)))
    #return (list (unlist(cluster_idx) ))
  }
  
  ret = c (k_ , num)
  names (ret) = c ("clusters" , "number of genes")
  # when plotting return cutree O.W. return ret
  return (list (ret) )
  #return (list (cutree(res_hclust , h = 0.4) ) )
}

file = "~/Farnush/farnush/Fabio Rossi/cell cell comunication/consensus clustering/average 60/"
file = "~/Documents/Farnush GitHub/Fabio Rossi/cell cell comunication/consensus clustering/average 60/"
#length(kegg_entries_receptors)
res = list()
names = c()
for (i in 1:length(kegg_entries_receptors))
{# do the downstream analysis for alll receptors 
  print (paste ("******", i , sep = ":  "))
  print(kegg_entries_receptors[i])
  print(receptor_symbol[i])
  names = c(names , receptor_symbol[i])
  res[[i]] <-  downstreamGenes_consensus(kegg_entries_receptors[i] , receptor_symbol[i] , file)
}
names (res) <- names
save (res , file = "consensus clustering/average 60/all_inf_res.RData")
save (res , file = "consensus clustering/average 60/number_cluster.RData")

##############################################################################3
hierarchicalClusteringAnalysis <- function (cluster_res)
{
  
  (which (cluster_res == "NA") %>% unname() -> idx)
  cluster_res [ -idx] -> cluster_res
  plotDistK_genes (cluster_res , "")
  (findSmallClusters (cluster_res , 5 , 1) -> c_small)
  
}

makeTableOfClusteringRes <- function (all_inf_res)
{
  toWrite = c()
  (which (all_inf_res == "NA") -> idx )
  all_inf_res [ -idx] -> all_inf_res
  for (i in 1: length(all_inf_res))
  {
    (recep_name <- all_inf_res[i] %>% names())
    for (j in 1: (all_inf_res[[i]] %>% length()) )
    {
      print (paste ( i , j , sep = ":"))
      (all_inf_res[[i]] [[j]] %>% unlist() -> r)
      (r %>% unlist() %>% unique() %>% length() -> K)
      (path <- all_inf_res[[i]][j] %>% names())
      for ( k in 1:K)
      {
        (which (r == k) %>% names() -> g)
        if (length(g) > 2)
          toWrite = rbind (toWrite , c ( recep_name, path , paste(g, collapse = ",")))
        #plotCluster(table , g %>% toupper(), paste(i , cellType,sep = "/") , dayNames,jj) -> p
        #plist [[i]] <- p
      }     
    }
  }
  colnames(toWrite) = c ("receptor" , "pathway_id" , "genes_in_cluster")
  
  write.table(toWrite , file = "consensus clustering/average 60/table_clustering_atleast3.txt" ,row.names = FALSE)
}


makeTableOfClusteringRes (res)

# 
# # prepare cells for finding the ligand 
# dataQC_ligand(EC_WT,3) -> EC_WT_ligand
# dataQC_ligand(EC_damaged,3) -> EC_KO_ligand
# dataQC_ligand(FAP_WT,3) -> FAP_WT_ligand
# dataQC_ligand(FAP_damaged,3) -> FAP_damaged_ligand

read.delim("interactions.txt" , header = TRUE) -> interactions
read.table (file = "consensus clustering/average 60/table_clustering.txt" , header = TRUE) -> rec_path_clusters
read.table (file = "consensus clustering/average 60/table_clustering_atleast2.txt" , header = TRUE) -> rec_path_clusters


cellToCell_ligand_clusters <- function(th)
{
  plist = list()
  list_i = 1
  khar = 1
  cellTypes <- c ("EC_WT" , "EC_KO" , "FAP_WT" , "FAP_KO")
  #dim(interactions)[1]
  for (i in 401: 878)
  {
    (interactions[i, ]$ligand_symbol %>% as.character() -> li)
    (interactions[i, ]$receptor_symbol %>% as.character() -> rec)
    (which (rec_path_clusters$receptor == rec ) -> rec_idx_table)
    if (length(rec_idx_table) > 0)
    {
      rec_path_clusters$genes_in_cluster [ rec_idx_table ] -> clusters
      rec_path_clusters$pathway_id [ rec_idx_table ] -> paths
      #print (i)
      for (j in 1:length(clusters))
      {
        
        (strsplit(clusters[j] %>% as.character() , ",") %>% unlist() -> cluster_symbol)
        if (length(cluster_symbol) > 0 )
        {# not recording the group with only one gene
          for (j1 in 1:length(cellTypes))
          {# cluster cell
            for (j2 in 1: length(cellTypes))
            {# li cell
              #print (paste (j, paste (j1 , j2 , sep = ":") , sep = ":") )
              (s <- corLigand_cluster ( cluster_symbol , li %>% toupper() , cellTypes[j1] , cellTypes[j2] ))
              if (s >= th)
              {
                print (s)
                print (paste ( i , j , sep = ":"))
                khar = khar + 1
                plist [[list_i]] <- c (li ,cellTypes[j2], rec , clusters[j] %>% as.character() , cellTypes[j1],paths[j] %>% as.character(), s)
                #print (paste0(c (li ,cellTypes[j2], rec , clusters[j] %>% as.character() , cellTypes[j1], s)))
                list_i = list_i + 1
              }
            }
          } 
        }
      } 
    }
  }
  return (plist)
  # test if clusters and the ligand would correlate
}

cellToCell_ligand_clusters (0.6) -> cell_cell_06
#cellToCell_ligand_clusters (0.8) -> cell_cell_08_part2
#cellToCell_ligand_clusters (0.8) -> cell_cell_08_ta800
#cell_cell_08_ta400 <- c (cell_cell_08 , cell_cell_08_part2)
cell_cell_08 = c (cell_cell_08_ta400 , cell_cell_08_ta800)
save (cell_cell_08_ta400 , file = "consensus clustering/average 60/ligand_downstream_ta400.RData")


findSymbol <- function (symbol, table)
{
  which (table$tracking_id %>% toupper() == symbol) -> idx
  if (length(idx) > 1)
    idx = idx[1]
  return (idx)
}

corTwoTablesByCommonCols <- function(table , vec)
{ # this function removes the days that are not available in either tables. 
  # then correlate based on the remaining days 
  # examples : D2 D3 D5 D7 and D2 D4 D5 D7 would be D2 D5 D7
  n1 <- colnames(table)
  n2 <- colnames(vec)
  (intersect(n1 , n2 ) -> n_c)
  if (dim(table)[1] == 1){
    (cor (table [ , n_c] , vec[ , n_c] ) -> m) 
  }else{
    (cor (t(table [ , n_c]) , vec[ , n_c] ) -> m)
  }
  
  return (mean(m) )
}

#cluster_cell and ligand_cell could be "EC_WT" , "EC_KO" , "FAP_KO" and "FAP_WT"
corLigand_cluster <- function ( cluster_symbol , ligand_symbol , cluster_cell , ligand_cell )
{
  # I find the ligands through the common genes 
  
  # finding cluster right patterns 
  if (cluster_cell == "EC_WT")
  {
    sapply(cluster_symbol %>% toupper(), findSymbol , EC_WT_rep1_l_f_n_c ) -> cluster_idx
    cluster_exp <- EC_WT_rep1_l_f_n_c [ unlist(cluster_idx) , -c(1,2)] %>% as.matrix()
  }else if (cluster_cell == "EC_KO"){
    sapply(cluster_symbol %>% toupper(), findSymbol , EC_KO_rep1_l_f_n_c ) -> cluster_idx
    cluster_exp <- EC_KO_rep1_l_f_n_c [ unlist(cluster_idx) , -c(1,2) ] %>% as.matrix()
  }else if (cluster_cell == "FAP_WT"){
    sapply(cluster_symbol %>% toupper(), findSymbol , FAP_WT_rep1_l_f_n_c ) -> cluster_idx
    cluster_exp <- FAP_WT_rep1_l_f_n_c [ unlist(cluster_idx) , -c(1,2)] %>% as.matrix()
  }else if (cluster_cell == "FAP_KO"){
    sapply(cluster_symbol %>% toupper(), findSymbol , FAP_KO_rep1_l_f_n_c ) -> cluster_idx
    cluster_exp <- FAP_KO_rep1_l_f_n_c [ unlist(cluster_idx) , - c(1,2)] %>% as.matrix()
  } else {
    print ("Wrong cluster cell")
    return("NA")
  }
  
  
  # finding ligand right pattern
  if (ligand_cell == "EC_WT")
  {
    (which (EC_WT_rep1_l_f_n_c$tracking_id %>% toupper() == ligand_symbol) -> idx)
    if (length(idx) == 0)
      return (-1)
    EC_WT_rep1_l_f_n_c [idx , - c(1,2)] %>% as.matrix() -> ligand_exp 
  } else if (ligand_cell == "EC_KO"){
    (which (EC_KO_rep1_l_f_n_c$tracking_id %>% toupper() ==ligand_symbol) -> idx)
    if (length(idx) == 0)
      return (-1)
    EC_KO_rep1_l_f_n_c [idx , - c(1,2)] %>% as.matrix() -> ligand_exp
  } else if (ligand_cell == "FAP_WT"){
    (which (FAP_WT_rep1_l_f_n_c$tracking_id %>% toupper() ==ligand_symbol) -> idx)
    if (length(idx) == 0)
      return (-1)
    FAP_WT_rep1_l_f_n_c [idx , - c(1,2)] %>% as.matrix() -> ligand_exp
  } else if (ligand_cell == "FAP_KO"){
    (which (FAP_KO_rep1_l_f_n_c$tracking_id %>% toupper() ==ligand_symbol) -> idx)
    if (length(idx) == 0)
      return (-1)
    FAP_KO_rep1_l_f_n_c [ idx , - c (1,2)] %>% as.matrix() -> ligand_exp
  } else {
      print ("Wrong ligand cell")
      return (-1)
  }
  
  (cor_score <- corTwoTablesByCommonCols (cluster_exp , ligand_exp ))
  return(abs (cor_score) )
    
}

l = c()
for (i in 1: length(cell_cell_08))
{
  strsplit(cell_cell_08 [[i]] [4] , ",") %>% unlist %>% length() -> cluster_l
  l = c ( l , cluster_l)
}
hist(l[which(l>5)] , breaks = 50 , col = "green" )
hist(l[which(3<l)] , breaks = 50 , col = "green" )
which(l==2) %>% length()
which(l>10) %>% length()
which(l>15) %>% length()

which (l >10) -> l10
pdf (file = "consensus clustering/average 60/li_clusters_80_original.pdf")
for (i in 1: length(l10) )
{
  # c (li ,cellTypes[j2], rec , clusters[j] %>% as.character() , cellTypes[j1],paths[j] %>% as.character(), s)
  (cell_cell_08_ta400[[l10[i]]] -> c)
  cluster <- strsplit(c[4] , ",") %>% unlist() %>% toupper()
li <- c[1] %>% toupper()
  if (c[2] == "EC_WT")
  {
    table_l = EC_WT_rep1_l_f_n_c
    days_l = EC_WT_days_rep1
  }else if (c[2] == "EC_KO"){
    table_l = EC_KO_rep1_l_f_n_c
    days_l = EC_KO_days_rep1
  }else if (c[2] == "FAP_WT"){
    table_l = FAP_WT_rep1_l_f_n_c
    days_l = FAP_WT_days_rep1
  }else {
    table_l = FAP_KO_rep1_l_f_n_c
    days_l = FAP_KO_days_rep1
  } 
  
  if (c[5] == "EC_WT")
  {
    table_c = EC_WT_rep1_l_f_n_c
    days_c = EC_WT_days_rep1
  }else if (c[5] == "EC_KO"){
    table_c = EC_KO_rep1_l_f_n_c
    days_c = EC_KO_days_rep1
  }else if (c[5] == "FAP_WT"){
    table_c = FAP_WT_rep1_l_f_n_c
    days_c = FAP_WT_days_rep1
  }else {
    table_c = FAP_KO_rep1_l_f_n_c
    days_c = FAP_KO_days_rep1
  }
  plotCluster(table_c ,  cluster, paste (c[3] , c[5], sep = "/") , days_c , 3) -> p1
  plotCluster(table_l ,  li, paste (c[1] , c[2] , sep = "/") , days_l , 3) -> p2
  multiplot(p1 , p2 , cols =1 )
}
dev.off()
# jj = 3
# 
# cor (t(as.matrix(all_EC_FAP_rep1_log_std_filtered_n[ , -c(1:(jj-1))]) ) , 
# t(as.matrix(all_EC_FAP_rep1_log_std_filtered_n[ , -c(1:(jj-1))])) ) -> dist_all
# 
# cor (t(as.matrix(EC_WT_log_c_n[ , -c(1:(jj-1))]) ) , 
# t(as.matrix(EC_WT_log_c_n[ , -c(1:(jj-1))])) ) -> dist_ec_wt
# 
# cor (t(as.matrix(EC_KO_log_c_n[ , -c(1:(jj-1))]) ) , 
# t(as.matrix(EC_KO_log_c_n[ , -c(1:(jj-1))])) ) -> dist_ec_ko
# 
# cor (t(as.matrix(FAP_WT_log_c_n[ , -c(1:(jj-1))]) ) , 
#      t(as.matrix(FAP_WT_log_c_n[ , -c(1:(jj-1))])) ) -> dist_fap_wt
# 
# cor (t(as.matrix(FAP_KO_log_c_n[ , -c(1:(jj-1))]) ) , 
#      t(as.matrix(FAP_KO_log_c_n[ , -c(1:(jj-1))])) ) -> dist_fap_ko


# plotExpressionDist (allSamples , -0.5 , -0.3)
# threshold is -0.3 when we standardize first

#filterLowExpressedGenes(all_EC_FAP_rep1_log_std , -0.3 , 1) -> all_EC_FAP_rep1_log_std_filtered