hierarchical_clustering <- function (genes , table , jj  , title)
{
  #genes <- genes %>% unlist() %>% unname()
  ids <- toupper (table$tracking_id)
  (sapply(genes, function(x) {which ( ids== x) })   -> cluster_idx  )
  if  (length(unlist (cluster_idx)) == 0) {
    print ("there is no available gene in the table")
    print(title)
    return("NA")}else if (length(unlist (cluster_idx)) > 2){
    table_genes <- table [ Reduce(union , cluster_idx), -c(1:(jj-1))] %>% as.matrix() 
      
    #constructing dist matrix
    cor (t(as.matrix(table_genes[ , -c(1:(jj-1))]) ) , t(as.matrix(table_genes[ , -c(1:(jj-1))])) ) -> dist
    dist_m = 1 - dist
    #dist_m = 1 - abs (dist)
    colnames(dist_m) = rownames(dist_m) = table$tracking_id[Reduce(union , cluster_idx)] %>% as.character()
    hclust(as.dist( dist_m ) , method  = "average" ) -> res_hclust
      par(cex=0.4, mar=c(5, 8, 4, 1))
      plot(res_hclust , xlab = "" , main = title , hang = -1)
      abline(h = 0.25, lty = 1 , col = "red")
      cutree(res_hclust , h = 0.25) %>% unlist() %>% unique() %>% length() -> k
      # 1.2 , 2 - 1.1 , 1.9 - 1 , 1.8
      text(1.4, 1.2 , paste ("clusters:" ,k, sep = " ") , col = "red")
      abline(h = 0.5, lty = 1 , col = "blue")
      cutree(res_hclust , h = 0.5) %>% unlist() %>% unique() %>% length() -> k_
      text(1.4, 1.1 , paste ("clusters:" ,k_, sep = " ") , col = "blue")
      text(1.4, 1 , paste ("number of genes:" ,dim(dist_m)[1], sep = " ") , col = "mediumorchid4")
    }
    ret = c (k , dim(dist_m)[1])
    names (ret) = c ("clusters" , "number of genes")
    return (ret)
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


