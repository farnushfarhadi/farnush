### heatmap for all samples 
cbind (EC_WT_rep1_l_f_n_c[ , -c(1,2)] , EC_KO_rep1_l_f_n_c[ , -c(1,2)] , FAP_WT_rep1_l_f_n_c[,-c(1,2)] , FAP_KO_rep1_l_f_n_c[,-c(1,2)] , 
       IC_WT_rep1_l_f_n_c[,-1] , MP_WT_rep1_l_f_n_c[,-1] , MP_KO_rep1_l_f_n_c[,-1]) -> allSamples


cbind (EC_WT_rep1[ , -c(1,2)] , EC_KO_rep1[ , -c(1,2)] , FAP_WT_rep1[,-c(1,2)] , FAP_KO_rep1[,-c(1,2)] , 
       IC_WT_rep1[,-1] , MP_WT_rep1[,-1] , MP_KO_rep1[,-1]) -> allSamples

cbind (EC_WT[ , -c(1,2)] , EC_damaged[ , -c(1,2)] , FAP_WT[,-c(1,2)] , FAP_damaged[,-c(1,2)] , 
       inflammatory_WT[,-1] , muscleProgenitors_WT[,-1] , muscleProgenitors_damaged[,-1]) -> allSamples
#heatmap of all samples 
cbind (EC_WT_rep1_l_f_n_c[, c(1,2)] , allSamples) -> allSamplestoHeat

cbind (EC_WT_rep1[, c(1,2)] , allSamples) -> allSamplestoHeat
cbind (EC_WT[, c(1,2)] , allSamples) -> allSamplestoHeat
logTransform(allSamplestoHeat , 3) -> allSamplestoHeat_log
filterLowExpressedGenes(allSamplestoHeat_log , 0.7 , 0.5 , 3) -> allSamplestoHeat_filtered # 2263 15 genes
QuantileNormalize (allSamplestoHeat_filtered , 3) -> allSamplestoHeat_filtered_n 
allSamples_names = c(paste("EC_WT",EC_WT_days_rep1),paste("EC_KO",EC_KO_days_rep1),
                     paste("FAP_WT",FAP_WT_days_rep1),paste("FAP_KO",FAP_KO_days_rep1),
                     paste("IC_WT",IC_WT_days_rep1),
                     paste("MP_WT",MP_WT_days_rep1),paste("MP_KO",MP_KO_days_rep1))

allSamples_names = c(paste("EC_WT",EC_WT_days),paste("EC_KO",EC_damaged_days),
                     paste("FAP_WT",FAP_WT_days),paste("FAP_WT",FAP_damaged_days),
                     paste("IC_WT",inflammatory_WT_days),
                     paste("MP_WT",muscleProgenitors_WT_days),paste("mus_KO",muscleProgenitors_damaged_days))

ss_corHeatmap (allSamplestoHeat_filtered_n , "all samples" ,allSamples_names,3)

png(filename='ic wt.png', width=1000, height=1000)
heatmap.2(cor_matrix,  scale = NULL , trace="none",  dendrogram = "column" ,
col=cols, cexCol=1.2, cexRow=1 , margins= c (8,8), srtCol=90 , 
main = toWrite )
graphics.off()
