# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
# Load the package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#Read in the female liver data set
# The following 3421 probe set were arrived at using the following steps
#1) reduce to the 8000 most varying, 2) 3600 most connected, 3) focus on unique genes
file = bzfile("data/cnew_liver_bxh_f2female_8000mvgenes_p3600_UNIQUE_tommodules.xls.bz2");
dat0=read.table(file, header=TRUE)
names(dat0)
# this contains information on the genes
datSummary=dat0[,c(1:8,144:150)]
# the following data frame contains
# the gene expression data: columns are genes, rows are arrays (samples)
datExprFemale <- t(dat0[,9:143])
no.samples <- dim(datExprFemale)[[1]]
dim(datExprFemale)
# Set the columns names to probe names
colnames(datExprFemale) = datSummary$substanceBXH
# This module assignment was obtained by Ghazalpour et al
colorsFemale = dat0$module

setwd("~/Documents/Farnush github/Fabio Rossi/WGCNA/data/")
file = bzfile("LiverMaleFromLiverFemale3600.csv.bz2");
data = read.csv(file, header = TRUE);
datExprMale = t(data[, substring(colnames(data), 1, 3)=="F2_"]);
colnames(datExprMale) = data$substanceBXH


#### calculation of module preservation

setLabels = c("Female", "Male");
multiExpr = list(Female = list(data = datExprFemale), Male = list(data = datExprMale));
multiColor = list(Female = colorsFemale);

# 45 min
system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} );

# Save the results
save(mp, file = "../code/modulePreservation.RData");
load (file = "../code/modulePreservation.RData")

### show module preservation plots
ref = 2
test = 1
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
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
# If plotting into a file, close it
dev.off();

