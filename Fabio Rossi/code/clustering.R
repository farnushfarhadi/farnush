
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

 plotCluster <- function (table , cluster , title , dayNames)
 {
   (sapply(cluster, function(x) {which (table$tracking_id == x)}) -> cluster_idx  )
   t = table [ Reduce(union , cluster_idx), -c(1,2)] %>% as.matrix() 
   # if (z)
   # {
   #   cluster.means <- apply(t, 1, mean)
   #   cluster.stdevs <- apply(t, 1, sd)
   #   
   #   for (i in 1:dim (t)[1])
   #   {
   #     t[i,] <-  (t[i,] - cluster.means[i]) / cluster.stdevs[i]
   #   }
   # }
   colnames(t) = dayNames
   
   n  = table[Reduce(union , cluster_idx),1] %>% as.character();
   normal.table = data.frame( gene = n , t);
   normalized.intreshape <- melt(normal.table, id.var = "gene", variable_name = "Time_point")
   
   feature_number <- dim (normal.table ) [1]
   p1 <- ggplot(normalized.intreshape, aes(x = Time_point, y = value, colour = as.character(gene)) ) +
     geom_point() +
     geom_line(aes(group = as.character(gene))) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1) , legend.position="none") +
     ggtitle (paste (feature_number , " genes - " ,  paste (title , "original" , sep = "-")) )
   
   cluster.means <- apply(t, 1, mean)
   cluster.stdevs <- apply(t, 1, sd)
   
   for (i in 1:dim (t)[1])
   {
     t[i,] <-  (t[i,] - cluster.means[i]) / cluster.stdevs[i]
   }
   
   n  = table[Reduce(union , cluster_idx),1] %>% as.character();
   normal.table = data.frame( gene = n , t);
   normalized.intreshape <- melt(normal.table, id.var = "gene", variable_name = "Time_point")
   
   feature_number <- dim (normal.table ) [1]
   p2 <- ggplot(normalized.intreshape, aes(x = Time_point, y = value, colour = as.character(gene)) ) +
     geom_point() +
     geom_line(aes(group = as.character(gene))) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1) , legend.position="none") +
     ggtitle (paste (feature_number , " genes - " ,  paste (title , "standardized" , sep = "-")) )
   
   multiplot(p1 , p2 , rows = 2)
   #return (p)
 }
 
 
 method (EC_damaged_log_filtered , 0.95 ) -> res95
 res95_num <- sapply(res95 , function(x) {length (x)}) 
 
 method (EC_damaged_log_filtered , 0.90 ) -> res90
 res90_num <- sapply(res90 , function(x) {length (x)}) 
 
 res95_num  %>% order() %>% tail (30) -> idx
 res95_num [idx]
 
 

 
 method (EC_damaged_log_filtered , 0.90 ) -> res90_notRep
 res90_notRep_num <- sapply(res90_notRep , function(x) {length (x)}) 
 
res90_notRep_num  %>% order() %>% tail (30) -> idx
res90_notRep_num [idx]

plotTopdf <- function (path , table , cluster_res , title , num , dayNames)
{
  cluster_res_num <- sapply(cluster_res , function(x) {length (x)}) 
  cluster_res_num  %>% order() %>% tail (num) -> idx
  
  pdf(file=path) 
  for (i in num:1)
  {
    print (plotCluster( table , cluster_res[[  idx[i] ]]  , title , dayNames) )
  }
  dev.off()
}
somePDFPath = "~/Desktop/Fabio-summer 2016/code/EC_clustering95notRep.pdf"


#################

#### EC_damaged
#method (EC_damaged_log_filtered , 0.95 ) -> res95_notRep
EC_damaged_days <- c("D0" , "D1" , "D2" , "D3" , "D5" , "D6" , "D7" , "D10")
load ("../../code/clustering res/method1/90EC_damaged_noRep.RData")
load("../../code/clustering res/method1/95EC_damaged_noRep.RData")
res95_notRep_num <- sapply(res95_notRep , function(x) {length (x)}) 
res95_notRep_num  %>% order() %>% tail (50) -> idx
res95_notRep_num [idx]
logTransform(EC_damaged) -> EC_damaged_log
plotCluster(EC_damaged_log , res95_notRep[[tail(idx ,1)]] , "standardized" , 1) -> pn
plotCluster(EC_damaged_log , res95_notRep[[tail(idx ,1)]] , "original" , 0) -> p
multiplot(pn , p , cols = 2)
plotTopdf (somePDFPath ,EC_damaged_log , res95_notRep , "95%" , 100 )

#### FAP_damaged
method (FAP_damaged_log_filtered , 0.90 ) -> res90_FAP_noRep
FAP_res90_notRep_num <- sapply(res90_FAP_noRep , function(x) {length (x)}) 
FAP_res90_notRep_num  %>% order() %>% tail (30) -> idx
FAP_res90_notRep_num [idx]
somePDFPath = "~/Desktop/Fabio-summer 2016/code/FAP_clustering90notRep.pdf"
plotTopdf (somePDFPath ,FAP_damaged_log , res90_FAP_noRep , "90%" , 50 )
setwd("~/Desktop/Fabio-summer 2016/code/clustering res/")

## combning clusters?

#### EC WT
method (EC_WT_log_filtered_n , 0.90 ) -> EC_WT_m1_90
EC_WT_m1_90_size <- sapply(EC_WT_m1_90 , function(x) {length (x)}) 
EC_WT_m1_90_size  %>% order() %>% tail (100) -> idx
EC_WT_m1_90_size [idx]
EC_WT_Path = "../../code/clustering res/method1/EC/EC_WT_90.pdf"
EC_WT_days <- c("D0" , "D2-1" , "D2-2" , "D2-3" , "D3-1" , "D3-2" , "D3-3" , "D4" , "D5" , "D6" , "D7" , "D10-1" , "D10-2" , "D14" )
plotTopdf (EC_WT_Path ,EC_WT_log_filtered_n , EC_WT_m1_90 , "90%" , 100 ,  EC_WT_days)
EC_WT_log_filtered_n_avg = data.frame()

