datExprAll <- prepareDataWGCNA(EC_WT_rep1[1:2000,] , EC_KO_rep1[1:2000,] , 3 , EC_WT_days_rep1 , EC_KO_days_rep1 ) 
wt_expr <- datExprAll[[1]]
damaged_expr <- datExprAll[[2]]
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft_wt = pickSoftThreshold(wt_expr, powerVector = powers, verbose = 5) # 26
sft_damaged = pickSoftThreshold(damaged_expr, powerVector = powers, verbose = 5) # 9
plot_networkTopology (sft_wt , 0.75) # 28
plot_networkTopology (sft_damaged , 0.75) # 9

net_wt = blockwiseModules(wt_expr, power = 18, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "femaleMouseTOM", verbose = 3)

net_damaged = blockwiseModules(damaged_expr, power = 9 , TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
                               numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "femaleMouseTOM", verbose = 3)

