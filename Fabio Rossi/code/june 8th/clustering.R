library (Rmisc)
method <- function (table , threshold )
{
  #print (paste ("Analyzing group" , group_number))
  
  #compute the correlation for the genes in the cell
  cor (t(as.matrix(table[ , -c(1,2)]) ) , t(as.matrix(table[ , -c(1,2)])) ) -> cor_matrix 
  genes_l <- dim(table)[1]
  #genes_l <- l
  gene_ID <- table$tracking_id %>% as.character()
  colnames(cor_matrix) <- gene_ID
  #colnames(cor_matrix) <- gene_ID[1:l]
  rownames(cor_matrix) <- gene_ID
  #rownames(cor_matrix) <- gene_ID[1:l]
  
  #initialization
  clusters = list();
  (clusters$set1 = c(gene_ID[1]))
  if (dim(cor_matrix)[1] > 1)
  {
    for (i in 2: genes_l)
    {
      print(i)
      (newFeature <- gene_ID[i])
      inserted = 0 ;
      for (j in 1: length (clusters))
      {
        (set <- paste0("set",j)) # clusters[set]
        featuresInSet <- unname (unlist (clusters[set]) )
        if (all (cor_matrix [ newFeature , featuresInSet] > threshold ) ) 
        {
          # if this feature showed high correlation with all features in the set, add it to the set
          print ("inserted")
          inserted = 1;
          (clusters[[set]] <- c(clusters[[set]] , newFeature))
          # letting the feature be insertd in all the groups that could have it
          break;
        }
      }
      # otherwise, define a new set
      if (!inserted)
      {
        print("new cluster")
        (set <- paste0("set",length(clusters)+1) )
        (clusters[set] <- c(newFeature) )
      }
    }
    subgroups_num <- sapply(clusters , function(x) {length (x)})   
    setsToPlot <- clusters[ which (subgroups_num > 1)]
    return (clusters)
    # setToPlots contain only subgroups with more than one feature
    # numberlist <- 1:length(setsToPlot)
    # lapply(numberlist, plotter , setsToPlot = setsToPlot , group = group , group_number= group_number)
    # 
    #subgroups contain ALL subgroups with one or more than one features
    
  }
  #return (subgroups)
 }
 
 
 plotCluster <- function (table , cluster , title)
 {
   sapply(cluster, function(x) {which (table$tracking_id == x)}) -> cluster_idx  
   t = table [ cluster_idx %>% unname(), -c(1,2)] %>% as.matrix() 
   cluster.means <- apply(t, 1, mean)
   cluster.stdevs <- apply(t, 1, sd)
   
   for (i in 1:dim (t)[1])
   {
     t[i,] <-  (t[i,] - cluster.means[i]) / cluster.stdevs[i]
   }
   n  = table[cluster_idx %>% unname(),1] %>% as.character();
   normal.table = data.frame( gene = n , t);
   normalized.intreshape <- melt(normal.table, id.var = "gene", variable_name = "Time_point")
   
   feature_number <- dim (normal.table ) [1]
   p <- ggplot(normalized.intreshape, aes(x = Time_point, y = value, colour = as.character(gene)) ) +
     geom_point() +
     geom_line(aes(group = as.character(gene))) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1) , legend.position="none") +
     ggtitle (paste (feature_number , " genes - " ,  title))
   return (p)
 }
 
 
 method (EC_damaged_log_filtered , 0.95 ) -> res95
 res95_num <- sapply(res95 , function(x) {length (x)}) 
 
 method (EC_damaged_log_filtered , 0.90 ) -> res90
 res90_num <- sapply(res90 , function(x) {length (x)}) 
 
 res95_num  %>% order() %>% tail (30) -> idx
 res95_num [idx]
 
 plotCluster (EC_damaged_log , res95[[357]] , "95%") -> p1
 plotCluster (EC_damaged_log , res95[[112]] , "95%") -> p2
 plotCluster (EC_damaged_log , res95[[596]] , "95%") -> p3 # P1 P3 have genes in common 12
 plotCluster (EC_damaged_log , res95[[472]] , "95%") -> p4
 multiplot (p1 , p2 , p3 , p4 , cols = 2)
 
 
 method (EC_damaged_log_filtered , 0.95 ) -> res95_notRep
 res95_notRep_num <- sapply(res95_notRep , function(x) {length (x)}) 
 
 method (EC_damaged_log_filtered , 0.90 ) -> res90_notRep
 res90_notRep_num <- sapply(res90_notRep , function(x) {length (x)}) 
 
 method (FAP_damaged_log_filtered , 0.95 ) -> FAPres95_notRep
 res95_notRep_num  %>% order() %>% tail (30) -> idx
 res95_notRep_num [idx]
 plotCluster (EC_damaged_log , res95_notRep[[112]] , "95%") -> p1
 plotCluster (EC_damaged_log , res95_notRep[[47]] , "95%") -> p2
 plotCluster (EC_damaged_log , res95_notRep[[91]] , "95%") -> p3
 plotCluster (EC_damaged_log , res95_notRep[[5]] , "95%") -> p4
 plotCluster (EC_damaged_log , res95_notRep[[541]] , "95%") -> p5
 plotCluster (EC_damaged_log , res95_notRep[[635]] , "95%") -> p6
 plotCluster (EC_damaged_log , res95_notRep[[180]] , "95%") -> p7
 plotCluster (EC_damaged_log , res95_notRep[[122]] , "95%") -> p8
 multiplot (p1 , p2 , p3 , p4 , cols = 2  )
 multiplot (p5 , p6 , p7 , p8 , cols = 2  )
 
 res90_notRep_num  %>% order() %>% tail (30) -> idx
 res90_notRep_num [idx]
 
 kmeanpp <- function(table )
 {
   
 }