source("http://bioconductor.org/biocLite.R") 
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
install.packages("WGCNA")
library (WGCNA)

workingDir = ".";
setwd(workingDir);
options(stringsAsFactors = FALSE )

setwd("../data/saved in R /")
load (file = "EC_WT.RData")
load (file = "EC_damaged.RData")
load (file = "FAP_WT.RData")
load (file = "FAP_damaged.RData")
load (file = "inflammatory_WT.RData")
load (file = "muscleProgenitors_WT.RData")
load (file = "muscleProgenitors_damaged.RData")

#EC_all <- cbind (EC_WT , EC_damaged [ , - c(1,2)])
prepareDataWGCNA <- function (wt , ko , idx)
{
  all <- cbind (wt , ko [ , - c(1,2)])
  filterLowExpressedGenes(logTransform(all) , 0.7 , 1) -> all_log_filtered
  
  #remove isoforms
  sapply(all_log_filtered$tracking_id %>% unique() %>% as.character(), function(x) {which (all_log_filtered$tracking_id == x)} ) -> hh 
  k = c()
  for (i in 1:length(hh))
  {
    k = c(k , hh[[i]][1])
  }
  # remove dupliated rows
  wt_expr <- all_log_filtered [ k, idx:dim(wt)[2]] # 8711 14
  damaged_expr <- all_log_filtered [ k, - c(1:dim(wt)[2])] # 8711 8
  
  rownames(wt_expr ) <- all_log_filtered$tracking_id %>% as.character() %>% unique()
  colnames(wt_expr) <- colnames(all_log_filtered) [idx:dim(wt)[2]]
  rownames(damaged_expr ) <- all_log_filtered$tracking_id %>% as.character() %>% unique()
  colnames(damaged_expr) <- colnames(all_log_filtered) [( - c(1:dim(wt)[2]))]
  
  
  
  wt_expr <- as.data.frame(t(wt_expr) )
  damaged_expr <- as.data.frame(t(damaged_expr))
  
  return (list (wt_expr , damaged_expr))
  
}

### soft thresholding 
#datExprAll <- prepareDataWGCNA(EC_WT , EC_damaged , 3) # 8711
#datExprAll <- prepareDataWGCNA(FAP_WT , FAP_damaged , 3) # 7985
datExprAll <- prepareDataWGCNA(muscleProgenitors_WT , muscleProgenitors_damaged , 2) # 8563

wt_expr <- datExprAll[[1]]
damaged_expr <- datExprAll[[2]]

powers = c(c(1:10), seq(from = 12, to=20, by=2))

powers = c(1:25)
powers = seq (from = 115900 , to = 116020 , by = 5)
# Call the network topology analysis function

sft_wt = pickSoftThreshold(wt_expr, powerVector = powers, verbose = 5)
sft_damaged = pickSoftThreshold(damaged_expr, powerVector = powers, verbose = 5)

save (sft_wt ,file= "../../WGCNA/moulePreservation/EC/sft_EC_wt.RData")
save (sft_damaged , file =  "../../WGCNA/moulePreservation/EC/sft_EC_ko.RData")

save (sft_wt ,file= "../../WGCNA/moulePreservation/FAP/sft_FAP_wt.RData")
save (sft_damaged , file =  "../../WGCNA/moulePreservation/FAP/sft_FAP_ko.RData")

save (sft_wt ,file= "../../WGCNA/moulePreservation/muscle/sft_muscle_wt.RData")
save (sft_damaged , file =  "../../WGCNA/moulePreservation/muscle/sft_muscle_ko.RData")

# load(file = "../../WGCNA/moulePreservation/FAP/sft_FAP_wt.RData")
# load(file = "../../WGCNA/moulePreservation/FAP/sft_FAP_ko.RData")
load(file = "../../WGCNA/moulePreservation/muscle/sft_muscle_wt.RData")
load(file = "../../WGCNA/moulePreservation/muscle/sft_muscle_ko.RData")

#FAP 25 - EC 13
plot_networkTopology (sft_wt)
#FAP 22 - EC 16
plot_networkTopology (sft_damaged) 

net_wt = blockwiseModules(wt_expr, power = 25, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "femaleMouseTOM", verbose = 3)

net_damaged = blockwiseModules(damaged_expr, power = 22 , TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
                             numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "femaleMouseTOM", verbose = 3)

save (net_wt ,file= "../../WGCNA/moulePreservation/EC/net_EC_wt.RData")
save (net_wt ,file= "../../WGCNA/moulePreservation/EC/net_EC_ko.RData")
load(file = "../../WGCNA/moulePreservation/EC/net_EC_wt.RData")
load(file = "../../WGCNA/moulePreservation/EC/net_EC_ko.RData")


save (net_wt ,file= "../../WGCNA/moulePreservation/FAP/net_FAP_wt.RData")
save (net_damaged , file =  "../../WGCNA/moulePreservation/FAP/net_FAP_ko.RData")
load(file = "../../WGCNA/moulePreservation/FAP/net_FAP_wt.RData")
load(file = "../../WGCNA/moulePreservation/FAP/net_FAP_ko.RData")

clusteringDendogramGenes(net_wt)
clusteringDendogramGenes(net_damaged)

saveIt (net_wt , "../../WGCNA/moulePreservation/EC/network_wt.RData")
saveIt (net_damaged , "../../WGCNA/moulePreservation/EC/network_ko.RData")

saveIt (net_wt , "../../WGCNA/moulePreservation/FAP/network_wt.RData")
saveIt (net_damaged , "../../WGCNA/moulePreservation/FAP/network_ko.RData")


##### loading the networks and preparing data for module preservation
nSets = 2;
multiExpr = list();
multiExpr[[1]] = list(data = wt_expr);
multiExpr[[2]] = list(data = damaged_expr);
setLabels = c("WT", "CCR2_KO");
names(multiExpr) = setLabels

x = load("../../WGCNA/moulePreservation/EC/network_wt.RData")
x = load("../../WGCNA/moulePreservation/FAP/network_wt.RData")
color_wt <- moduleColors
x = load("../../WGCNA/moulePreservation/FAP/network_ko.RData")
x = load("../../WGCNA/moulePreservation/EC/network_ko.RData")
color_damaged <- moduleColors
colorList = list(color_wt, color_damaged);
names(colorList) = setLabels;

## calculate module preservation statistics
# 82 min EC
# 63 min FAP
system.time( {
  mp = modulePreservation(multiExpr, colorList,referenceNetworks = c(1:2),loadPermutedStatistics = FALSE,nPermutations = 200,verbose = 3)
});

save(mp , file = "../../WGCNA/moulePreservation/EC/mp.RData")
load (file = "../../WGCNA/moulePreservation/EC/mp.RData")

save(mp , file = "../../WGCNA/moulePreservation/FAP/mp.RData")
load (file = "../../WGCNA/moulePreservation/FAP/mp.RData")

magenta_genes <- read.table("../../WGCNA/magenta_genes_EC.txt" , header = FALSE) %>% unlist() %>% unname()
days = c("D0" , "D1" , "D2"  , "D3" , "D5" , "D6" , "D7" , "D10" )
plotCluster(EC_damaged_log_filtered , magenta_genes, "preserved in EC CCR2 KO",days )
plotCluster(EC_WT_log_filtered_n_1 ,magenta_genes ,"preserved in EC WT",EC_WT_days1  )
# 
corAverage ( filterLowExpressedGenes( logTransform(EC_WT), 0.7 , 1), magenta_genes , 0) -> s
s[lower.tri(s, diag = FALSE) ]-> s
par (mfrow = c(1,1))
hist (s , breaks = 100 , col = "blue" , main = "Pearson cor - most preserved module in EC WT")
#

corAverage ( filterLowExpressedGenes( logTransform(EC_damaged), 0.7 , 1), magenta_genes , 0) -> s
s[lower.tri(s, diag = FALSE) ]-> s
par (mfrow = c(1,1))
hist (s , breaks = 100 , col = "blue" , main = "Pearson cor - most preserved module in EC CCR2 KO")

gy_genes <- read.table("../../WGCNA/gy_genes_EC.txt" , header = FALSE) %>% unlist() %>% unname()
days = c("D0" , "D1" , "D2"  , "D3" , "D5" , "D6" , "D7" , "D10" )
plotCluster(EC_damaged_log_filtered , gy_genes, "preserved in EC CCR2 KO",days )
plotCluster(EC_WT_log_filtered_n_1 ,gy_genes ,"preserved in EC WT",EC_WT_days1  )

setwd("../../WGCNA/tutorial")
sampleTree = hclust(dist (t(EC_WT[-c(1,2)]) ) , method = "average")
# sampleTree = hclust(dist (QuantileNormalize(logTransform(t(EC_WT[-c(1,2)]))) ) , method = "average")
# sampleTree = hclust(dist (logTransform(t(EC_WT[-c(1,2)]))) , method = "average")
sizeGrWindow(12, 9)
par (cex = 0.6)
par (mar = c(0,4,2,0))
plot (sampleTree , main = "sample clustering to detect outlier" , sub="" , xlab = "" , cex.lab = 1.5 , cex.axis=1.5 , cex.main = 2)
plot (sampleTree , main = "sample clustering to detect outlier" )

clusteringDendogramGenes <- function(net)
{
  # open a graphics window
  sizeGrWindow(12, 9)
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
}

# changed h in abline
plot_networkTopology <- function (sft)
{
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.8,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}

saveIt <- function(net , file)
{
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];
  save(MEs, moduleLabels, moduleColors, geneTree, 
       file = file)
}