hierarchical_clustering <- function (genes , table , jj  , title)
{
  genes <- genes %>% unlist() %>% unname()
  ids <- toupper (table$tracking_id)
  (sapply(genes, function(x) {which ( ids== x) })   -> cluster_idx  )
  if  (length(unlist (cluster_idx)) == 0) {
    print ("there is no available gene in the table")
    print(title)
    return("NA")}else{
    table_genes <- table [ Reduce(union , cluster_idx), -c(1:(jj-1))] %>% as.matrix() 
      
    #constructing dist matrix
    cor (t(as.matrix(table_genes[ , -c(1:(jj-1))]) ) , t(as.matrix(table_genes[ , -c(1:(jj-1))])) ) -> dist
    dist_m = 1 - dist
    #dist_m = 1 - abs (dist)
    colnames(dist_m) = rownames(dist_m) = table$tracking_id[Reduce(union , cluster_idx)] %>% as.character()
    }
  
    hclust(as.dist( dist_m ) , method  = "complete" ) -> res_hclust
    par(cex=0.4, mar=c(5, 8, 4, 1))
    plot(res_hclust , xlab = "" , main = title , hang = -1)
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


