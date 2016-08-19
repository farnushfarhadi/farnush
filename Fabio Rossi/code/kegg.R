source("https://bioconductor.org/biocLite.R")
biocLite("KEGGREST")
library(KEGGREST)


# procedure:
keggGet("mmu:22339") -> pekh 
pekh[[1]]$PATHWAY %>% names() -> p
keggGet(p[1]) -> path
path[[1]]$GENE # it will give me both kegg entry id (entrez) and all gene symbol, I will use entrez


# converting gene ID to entrez ID 
library (org.Hs.eg.db)
genes_symbol_entrez <- as.list(org.Hs.egALIAS2EG)
genes_symbol_entrez_n <- names (genes_symbol_entrez)

