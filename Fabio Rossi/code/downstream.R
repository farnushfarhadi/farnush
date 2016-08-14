setwd("~/Documents/Farnush github/Fabio Rossi/cell cell comunication/")
#read.table("FLAT_FILES_072010/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt" , header = FALSE ) -> ppi_human
read.delim("HPRD/FLAT_FILES_072010/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt" , header = FALSE ) -> ppi_human
ppi_human[ , 1:6]   -> ppi_human
colnames(ppi_human) <- c ("interactor_1_geneSymbol","interactor_1_hprd_id","interactor_1_refseq_id",
                          "interactor_2_geneSymbol","interactor_2_hprd_id","interactor_2_refseq_id")


read.delim("ligand-receptor.txt" , header = TRUE) -> ligand_receptor 

ligand_receptor$Receptor..Symbols. %>% as.character() %>% unique()-> receptor # 378 receptors

plotLigand_ReceptorDownstream <- function (receptor_gene , table1 , table2 , t1 , t2 , dayNames1 , dayNames2, jj1 , jj2)
{
  # extracting genes from ppi 
  which (ppi_human$interactor_1_geneSymbol == receptor_gene) -> idx1
  ppi_human$interactor_2_geneSymbol[idx1] %>% as.character() -> g1
  which (ppi_human$interactor_2_geneSymbol == receptor_gene) -> idx2
  ppi_human$interactor_1_geneSymbol[idx2] %>% as.character() -> g2
  
  g <- union(g1 , g2)
  
  # removing ligands from g to get downstream genes
  which(ligand_receptor$Receptor..Symbols. == receptor_gene) -> idx
  ligand_receptor$Ligand..Symbol. [idx] -> ligands
  ds_g <- setdiff(g , ligands)
  
  # plotting downstream genes and ligands
  plotCluster(table1 , ds_g , paste ("receptor downstream" , t1 , sep = "-") , dayNames1 , jj1) -> p1
  plotCluster(table2 , ligands , paste ("ligand" , t2 , sep = "-") , dayNames1 , jj2) -> p2
  
}