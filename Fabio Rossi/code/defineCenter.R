isPeak <- function (idx , expr , th)
{
  if ( (expr [idx] - expr[idx-1] > th) & (expr [idx] - expr[idx+1] > th ))
    return (TRUE)
  return(FALSE)
}

isBottom <- function (idx , expr , th)
{
  if ( (expr [idx-1] - expr[idx] > th) & (expr [idx+1] - expr[idx] > th) )
    return (TRUE)
  return(FALSE)
}

cluster_EC <- function(gene_symbol, th)
{
  (which (EC_WT_log_filtered_n_rep1$tracking_id %>% toupper() == gene_symbol) -> idx)
  if (length(idx)==0)
    return ("NA")
  if (length(idx) >1)
    idx = idx[1]
  EC_WT_log_filtered_n_rep1[idx , -c(1,2)] %>% as.numeric() -> expr 
  if (isPeak(3,expr , th) & isBottom(7 , expr , th))
  {
    #print("1")
    return("1")
  }
  if (isBottom(3,expr , th) & isBottom(7 , expr, th))
  {
    return("2")
  }
  if (isBottom(3,expr , th) & isPeak(7 , expr , th))
  {
    return("3")
  }
  if (isPeak(3,expr , th) & isPeak(7 , expr, th))
  {
    return("4")
  }
  return("NA")
}

downstreamGenes_v2 <- function (gene_entry , symbol, file)
{# for a given gene, it finds the kegg pathways that have the gene in 
  # and then cluster them hierarchialy for different cell types
  gene_entry_inf <- try(keggGet(gene_entry), silent=TRUE)
  if ('try-error' %in% class(gene_entry_inf)) {print ("NA")}
  if ('try-error' %in% class(gene_entry_inf)) {return ("NA")}
  gene_entry_inf[[1]]$PATHWAY %>% names() %>% unique()-> gene_entry_pathways
  if (gene_entry_pathways %>% length() == 0) {print ("NA")}
  if (gene_entry_pathways %>% length() == 0) {return ("NA")}
  sapply(gene_entry_pathways, function(x){getPathwayGenesEntrez(x)}) -> gene_entry_pathways_genes
  
  # for each pathway, have a plot
  #file = "~/Documents/Farnush github/Fabio Rossi/cell cell comunication/receptor_keg_plots/"
  #file = "~/Documents/Farnush github/Fabio Rossi/cell cell comunication/receptor_keg_hclust/"
  #file = "~/Farnush/farnush/Fabio Rossi/cell cell comunication/receptor_keg_hclust_average_60/"
  #file = "~/Documents/Farnush github/Fabio Rossi/cell cell comunication/receptor_keg_hclust_average/"
  path = paste (file , paste0 (symbol , ".pdf") , sep = "/")
  pdf(file=path)
  
  for ( i in 1:length(gene_entry_pathways_genes))
  {
    print (paste (length(gene_entry_pathways_genes) , i , sep = "-"))
    (paste (gene_entry_pathways[i],symbol, sep = "/") -> t)
    sapply(gene_entry_pathways_genes[[i]] %>% as.character(), cluster_EC,th = 0.3) -> res
    plist = list ()
    toPlot = 1
    for ( j in 1:4)
    {
      which (res == j) %>% names() -> g
      print (length(g) )
      if (length(g) <2)
      {
        toPlot = 0
        break;
      }
      plotCluster(EC_WT_log_filtered_n_rep1 , g %>% toupper(), 
      paste(gene_entry_pathways[i],j,sep="/") , EC_WT_days_rep1,jj) -> p
      plist [[j]] <- p
    }
    if (toPlot)
      do.call("grid.arrange", c(plist, ncol=2))
  }
  dev.off()
}

file = "~/Farnush/farnush/Fabio Rossi/cell cell comunication/define centers/"
#length(kegg_entries_receptors)
for (i in 102:length(kegg_entries_receptors))
{# do the downstream analysis for alll receptors 
  print (paste ("******", i , sep = ":  "))
  print(kegg_entries_receptors[i])
  print(receptor_symbol[i])
  downstreamGenes_v2(kegg_entries_receptors[i] , receptor_symbol[i] , file)
}
