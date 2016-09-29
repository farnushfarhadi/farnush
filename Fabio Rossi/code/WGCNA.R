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
prepareDataWGCNA <- function (wt , ko , idx , wt_d , ko_d)
{
  all <- cbind (wt , ko [ , - c(1:(idx-1))])
  filterLowExpressedGenes(logTransform(all , idx) , 0.7 , 0.5 , 1) -> all_log_filtered
  QuantileNormalize(all_log_filtered , idx) -> all_log_filtered
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
  colnames(wt_expr) <- paste ("wt" , wt_d , sep = "," )
  rownames(damaged_expr ) <- all_log_filtered$tracking_id %>% as.character() %>% unique()
  colnames(damaged_expr) <-  paste ("ko" , ko_d , sep = "," )
  
  
  
  wt_expr <- as.data.frame(t(wt_expr) )
  damaged_expr <- as.data.frame(t(damaged_expr))
  
  return (list (wt_expr , damaged_expr))
  
}

### soft thresholding 
#datExprAll <- prepareDataWGCNA(EC_WT , EC_damaged , 3) # 8711
#datExprAll <- prepareDataWGCNA(FAP_WT , FAP_damaged , 3) # 7985
datExprAll <- prepareDataWGCNA(muscleProgenitors_WT , muscleProgenitors_damaged , 2) # 8563
datExprAll <- prepareDataWGCNA(EC_WT_rep1 , EC_KO_rep1 , 3 , EC_WT_days_rep1 , EC_KO_days_rep1 ) 
datExprAll <- prepareDataWGCNA(FAP_WT_rep1 , FAP_KO_rep1 , 3 , FAP_WT_days_rep1 , FAP_KO_days_rep1 ) 


wt_expr <- datExprAll[[1]]
damaged_expr <- datExprAll[[2]]
save(wt_expr , file = "EC/wt_expr.RData")
save(damaged_expr , file = "EC/damaged_expr.RData")

save(wt_expr , file = "FAP/wt_expr.RData")
save(damaged_expr , file = "FAP/damaged_expr.RData")
powers = c(c(1:10), seq(from = 12, to=30, by=2))

powers = c(1:25)
#powers = seq (from = 115900 , to = 116020 , by = 5)
# Call the network topology analysis function
load ("FAP/wt_expr.RData")
load("FAP/damaged_expr.RData")
sft_wt = pickSoftThreshold(wt_expr, powerVector = powers, verbose = 5)
sft_damaged = pickSoftThreshold(damaged_expr, powerVector = powers, verbose = 5)

setwd("~/Documents/Farnush GitHub/Fabio Rossi/WGCNA/moulePreservation/")
save (sft_wt ,file= "EC/sft_EC_wt.RData")
save (sft_damaged , file =  "EC/sft_EC_ko.RData")

save (sft_wt ,file= "FAP/sft_FAP_wt.RData")
save (sft_damaged , file =  "FAP/sft_FAP_ko.RData")

save (sft_wt ,file= "muscle/sft_muscle_wt.RData")
save (sft_damaged , file =  "muscle/sft_muscle_ko.RData")

# load(file = "../../WGCNA/moulePreservation/FAP/sft_FAP_wt.RData")
# load(file = "../../WGCNA/moulePreservation/FAP/sft_FAP_ko.RData")
load(file = "../../WGCNA/moulePreservation/muscle/sft_muscle_wt.RData")
load(file = "../../WGCNA/moulePreservation/muscle/sft_muscle_ko.RData")

load(file = "EC/sft_EC_wt.RData")
load(file = "EC/sft_EC_ko.RData")
load(file = "FAP/sft_FAP_wt.RData")
load(file = "FAP/sft_FAP_ko.RData")
#FAP 28  0.85 - EC 26 0.80
plot_networkTopology (sft_wt , 0.80)
#FAP 16 0.85 - EC 22 0.80
plot_networkTopology (sft_damaged , 0.80) 

load(file = "EC/wt_expr.RData")
load(file = "EC/damaged_expr.RData")
net_wt = blockwiseModules(wt_expr, power = 26, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "femaleMouseTOM", verbose = 3)

net_damaged = blockwiseModules(damaged_expr, power = 22 , TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
                             numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "femaleMouseTOM", verbose = 3)

save (net_wt ,file= "EC/net_EC_wt.RData")
save (net_damaged ,file= "EC/net_EC_ko.RData")
load(file = "EC/net_EC_wt.RData")
load(file = "EC/net_EC_ko.RData")


save (net_wt ,file= "FAP/net_FAP_wt.RData")
save (net_damaged , file =  "FAP/net_FAP_ko.RData")
load(file = "FAP/net_FAP_wt.RData")
load(file = "FAP/net_FAP_ko.RData")

clusteringDendogramGenes(net_wt)
clusteringDendogramGenes(net_damaged)

saveIt (net_wt , "EC/network_wt.RData")
saveIt (net_damaged , "EC/network_ko.RData")

saveIt (net_wt , "FAP/network_wt.RData")
saveIt (net_damaged , "FAP/network_ko.RData")


##### loading the networks and preparing data for module preservation
nSets = 2;
multiExpr = list();
multiExpr[[1]] = list(data = wt_expr);
multiExpr[[2]] = list(data = damaged_expr);
setLabels = c("WT", "CCR2_KO");
names(multiExpr) = setLabels

load("EC/network_wt.RData")
load("FAP/network_wt.RData")
color_wt <- moduleColors
load("EC/network_ko.RData")
load("FAP/network_ko.RData")
color_damaged <- moduleColors
colorList = list(color_wt, color_damaged);
names(colorList) = setLabels;

## calculate module preservation statistics
# EC 128 min
#    FAP 93 min
system.time( {
  mp = modulePreservation(multiExpr, colorList,referenceNetworks = c(1:2),loadPermutedStatistics = FALSE,nPermutations = 200,verbose = 3)
});

save(mp , file = "EC/mp.RData")
load (file = "EC/mp.RData")

save(mp , file = "FAP/mp.RData")
load (file = "FAP/mp.RData")

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
plot_networkTopology <- function (sft , th)
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
  abline(h=th,col="red")
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

load ("EC/mp.RData")
load("EC/net_EC_ko.RData")
load ("EC/net_EC_wt.RData")
load ("EC/damaged_expr.RData")
load("EC/wt_expr.RData")
print_z_wgcna (mp , ref = 2 , test =1 , 6 ) -> res # test modules found in ko to be preserved in wt - 21
net_damaged$colors %>% unique() %>% length() # 31 modules detected in ec ko
length(color_damaged %>% unique()) # 31
saveModules(net = net_damaged , table = damaged_expr ,modulePresRes = res ) -> ec_notPresInWT

print_z_wgcna (mp , ref =1 , test =2 , 10 )
net_wt$colors %>% unique() %>% length()
load ("FAP/mp.RData")
load("FAP/net_FAP_ko.RData")
load ("FAP/net_FAP_wt.RData")
load ("FAP/damaged_expr.RData")
load("FAP/wt_expr.RData")

print_z_wgcna (mp , ref = 2 , test =1 , 6 ) -> res # modules found in ko
net_damaged$colors %>% unique() %>% length()

print_z_wgcna (mp , ref = 1 , test =2 , 6 ) -> res # modules found in ko
net_wt$colors %>% unique() %>% length()

saveModules(net = net_damaged , table = damaged_expr ,modulePresRes = res ) -> ec_notPresInWT

saveModules <- function (net , table , modulePresRes)
{
  (colors_not_pres <- rownames(modulePresRes[order (modulePresRes$Zsummary.pres),]) )# sorting by z summary  
  notPresModules <- c()
  l <- net$colors %>% unique()
  c  <- labels2colors (l)
  (sapply(colors_not_pres, function(x){which (c == x)}) -> idx )
  for (i in 1:length(idx))
  {
    which (net$colors == l[idx[[i]]]) -> g_idx
    colnames ( table ) [g_idx] -> genes
    notPresModules = rbind(notPresModules , c (paste (genes , collaps = ",") , length(genes) , modulePresRes$Zsummary.pres[i]) )
  }
  return (notPresModules)
}

### show module preservation plots
 

print_z_wgcna <- function (mp , ref , test , notPres)
{
  statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
  statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
  # Compare preservation to quality:
  res <- cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
               signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))
  #print( res )
  # Module labels and module sizes are also contained in the results
  
  (modColors = rownames(mp$preservation$observed[[ref]][[test]]) )
  (moduleSizes = mp$preservation$Z[[ref]][[test]][, 1] )
  # leave grey and gold modules out
  plotMods = !(modColors %in% c("grey", "gold"));
  # Text labels for points
  text = modColors[plotMods];
  # Auxiliary convenience variable
  plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
  # Main titles for the plot
  mains = c("Preservation Median rank", "Preservation Zsummary");
  # Start the plot
  sizeGrWindow(10, 5);
  #pdf(fi="Plots/BxHLiverFemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
  par(mfrow = c(1,2))
  par(mar = c(4.5,4.5,2.5,1))
  for (p in 1:2)
  {
    min = min(plotData[, p], na.rm = TRUE);
    max = max(plotData[, p], na.rm = TRUE);
    # Adjust ploting ranges appropriately
    if (p==2)
    {
      if (min > -max/10) min = -max/10
      ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
    } else
      ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
    plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
         main = mains[p],
         cex = 2.4,
         ylab = mains[p], xlab = "Module size", log = "x",
         ylim = ylim,
         xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
    labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
    # For Zsummary, add threshold lines
    if (p==2)
    {
      abline(h=0)
      abline(h=2, col = "blue", lty = 2)
      abline(h=10, col = "darkgreen", lty = 2)
    }
  }
  which (res$Zsummary.pres < notPres) -> idx
  print (modColors)
  #print (length(modColors))
  print (length(text))
  return (res [ idx,])
}


