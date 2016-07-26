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


file = bzfile("data/Dataset 1 (network construction).csv.bz2");
dat1 = read.csv(file, header=T)
dim(dat1)
names(dat1);
datExpr=data.frame(t(dat1[dat1$Brain_variant_H>0,2:39]))

indexHuman=c(19:36)
indexChimp=c(1:18)
# Number of data sets that we work with
nSets = 2;
# Object that will contain the expression data
multiExpr = list();
multiExpr[[1]] = list(data = datExpr[indexHuman, ]);
multiExpr[[2]] = list(data = datExpr[indexChimp, ]);
# Names for the two sets
setLabels = c("Human", "Chimp");
# Important: components of multiExpr must carry identificating names
names(multiExpr) = setLabels
# Display the dimensions of the expression data (if you are confused by this construct, ignore it):
lapply(multiExpr, lapply, dim)


x = load("data/HumanChimp-OldhamAnalysis-colorHuman-colorChimp-inNetwork.RData")
# Create an object (list) holding the module labels for each set:
colorList = list(colorHuman, colorChimp);
# Components of the list must be named so that the names can be matched to the names of multiExpr
names(colorList) = setLabels;
