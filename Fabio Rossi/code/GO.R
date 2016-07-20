source("http://bioconductor.org/biocLite.R")
biocLite("topGO")
biocLite("preprocessCore")
install.packages("GSA")
library (GSA)
library(topGO)
library(preprocessCore)
library(ALL)
GSA.read.gmt("../../GO/msigdb.v5.1.symbols.gmt") -> genes_annot
load ("../../code/clustering res/95EC_damaged_noRep.RData")
res95_notRep_num <- sapply(res95_notRep , function(x) {length (x)}) 
res95_notRep_num %>% order() %>% tail (30) -> idx
res95_notRep_num[idx]
res95_notRep[[112]]
write.table( res95_notRep[[112]] , file ="../../genes.txt", quote = FALSE , row.names = FALSE )
plotCluster(EC_damaged_log_filtered ,res95_notRep[[112]] , "hello" )
GOdata <- new("topGOdata", ontology = "MF", allGenes = res95_notRep[[112]],annot = annFUN.gene2GO, gene2GO = geneID2GO)

Reduce (union  , genes_annot$genesets[c(1:length(genes_annot$genesets))]) -> all_genes_annot
table$tracking_id %>% as.character() -> all_genes_data
#intersect (all_genes_annot , all_genes_data) -> all_genes
intersect (toupper (all_genes_data) , all_genes_annot ) -> all_genes   #8261



GO_analysis_clusters <- function (clusters , all_genes , genes_annot , table , go_t_l , go_t_h , c_size_l)
{
  # filtering very small or very large GO terms 
  sapply (genes_annot$genesets , length ) -> genes_annot_l
  (which (genes_annot_l < go_t_l) -> l_go_t_l ) %>% length()
  (which (genes_annot_l > go_t_h) -> l_go_t_h ) %>% length()
  genes_annot$genesets[- c(l_go_t_l , l_go_t_h)] -> genes_annot$genesets
  genes_annot$geneset.names[- c(l_go_t_l , l_go_t_h)] -> genes_annot$geneset.names
  
  #filtering small clusters clusters
  sapply(unique(clusters) , function(x){which (clusters == x) %>% length()}) -> clusters_size
  which(clusters_size < c_size_l) -> small_idx
  if (length(small_idx) >1 )
  {
    unique(clusters)[-small_idx] -> clusters_toTest
  }else{
    unique(clusters) -> clusters_toTest
  }
  
  res = matrix (nrow = clusters_toTest %>% length() , ncol = genes_annot$geneset.names %>% length())
  rownames(res) = clusters_toTest %>% as.character()
  colnames(res) = genes_annot$geneset.names
  # length (unique(clusters))
  multiple_tests = 0 ;
  for ( i in clusters_toTest %>% head(5))
  {
    intersect ( all_genes , table$tracking_id [which (clusters==i)] %>% as.character() %>% toupper() ) -> genes1 # genes of the cluster
    print (i)
    for (j in 1: (genes_annot$geneset.names %>% length()) )
    {
      (intersect ( all_genes , genes_annot$genesets[j] %>%  unlist() ) -> genes2 )# genes of the GO term
      (cYgN <- setdiff(genes1 , genes2) %>% length() )
      (cYgY <- intersect (genes1 , genes2) %>% length() ) 
      (cNgY <- setdiff(genes2 , genes1) %>% length() )
      (cNgN <- setdiff(all_genes , union(genes1 , genes2)) %>% length() )
      Ctable <- matrix(c(cYgY , cYgN , cNgY , cNgN) , nrow = 2 , byrow = TRUE)
      colnames(Ctable) = c("GO_y" , "GO _n")
      rownames(Ctable) = c("c_y" , "c_n")
      chisq.test(ftable(Ctable))-> test
      #print(paste(i , j , sep = " "))
      #print(test$p.value)
      if (is.na(test$p.value) )
      {
        res[i,j] = 0
      } else{
        if ( test$p.value< 0.05/ (length(clusters_toTest) * length(genes_annot$genesets)) )
        {
          res[i,j] = 1
        }else{
          res[i,j] =  0
        }
      }
      
    }
  }

  return (res)
}

GO_analysis_clusters (EC_WT_train1_kmean_Euc_norm[[1]] , all_genes , genes_annot , EC_WT_train1 , 10 , 1000 , 10 ) -> kmean74_Euc_norm_GO
clusters = EC_WT_train1_kmean_Euc_norm[[1]]
table = EC_WT_train1
go_t_l = 10 
go_t_h = 1000 
c_size_l = 10


apply(res_GO, 1, function(x){which(x==1) %>% length()}) -> res_GO_l
res_GO_l %>% unname() %>% hist()

colnames (res_GO)[which (res_GO[1, ] == 1) ] -> pekh
rownames(res_GO)[1] 
table$tracking_id [which (clusters==54)] %>% as.character() %>% toupper() -> g 
write.table(g , file= "../../pekh.txt" , row.names = FALSE, col.names = FALSE , quote = FALSE)


intersect(all_genes_annot , toupper(all_genes_annot)) %>% length()
all_genes_annot %>% length()
setdiff(all_genes_annot , intersect(all_genes_annot , toupper(all_genes_annot))) -> r
