library (splines)


#load("../../code/clustering res/90FAP_damaged_noRep.RData")
load("../../code/clustering res/90FAP_damaged_noRep.RData")
res90_FAP_noRep[[116]] -> genes
res90_FAP_noRep[[187]] -> genes # good 
res90_FAP_noRep[[282]] -> genes # 
plotCluster (FAP_damaged_log_filtered , genes , "FAP group 116")
sapply (genes , function(x) {which (FAP_damaged_log_filtered$tracking_id == x)} ) -> genes_i
genes_i %>% unname() -> i 
FAP_damaged_log_filtered [ i , -c(1,2)] -> fap
(fap[ , 4] + fap [,5]) /2 -> fap[,4]
fap[ , -5] -> fap

fap [1:5,] -> fap 
fap_t = c(0, 1, 2, 3, 4, 5, 6 , 10)
matplot( fap_t , t(fap) , type="l", lwd=2 )



data = data.frame(time = fap_t , geneExpression = as.numeric (fap[3,]) )
ispl <- interpSpline( data$time, data$geneExpression )
plot( predict( ispl, seq( 0, 14, length.out = 50 ) ), type = "l" )
points (data$time , data$geneExpression)

ispl <- polySpline(interpSpline( geneExpression ~ time,  data, bSpline = TRUE))
#plot(ispl )
plot( predict( ispl, seq( 0, 14, length.out = 50 ) ), type = "l" )
points (data$time , data$geneExpression)


## smooth spline 
attach (data)
plot(time, geneExpression, main = "data(cars)  &  smoothing splines")
smooth.spline(time, geneExpression, df = 5) -> res
predict (res , x = c(7 , 14))
lines(res, lty = 2, col = "red")
lines(smooth.spline(time, geneExpression, df = 7), lty = 2, col = "blue")
lines(smooth.spline(time, geneExpression, df = 4), lty = 2, col = "black")


(lm ( geneExpression ~ bs (time , 4)) -> m )
plot (geneExpression ~ time , pch = 17)
lines (geneExpression ~ time )
lines (predict (m) ~ time , lty = 2)


# table = table of available time points - has one replicates for each time point
# t_new = vector of time inclusing all available and missing time points , this is used for imputing 
# t_available = vector of time including time points available in our data 
# FAP : we have 0,1,2,3,4,5,6,10
# hide 1 and 5 day -> impute -> evaluate by comparing to real values 

fitSmoothingSpline <- function (table , idx_miss , df )
{
   
  #par(mfrow = c(3,3))
  t = c(0,1,2,3,4,5,6,10)
  # standardize row 
  # for (i in 1:dim(table)[1])
  # {
  #   table[i,] %>% as.matrix() %>% as.numeric() -> d
  #   table[i,] =  (d - mean (d) ) / sd (d) ;
  # }
  table_in = table [, -idx_miss]
  y = table [, idx_miss] %>% as.matrix() %>% as.numeric()
  (t_available = t[-idx_miss])
  (t_new = t[idx_miss])
  y_ = c()
  dim(table_in)[1]
  for (i in 1:dim(table_in)[1] )
  {

    (data = data.frame(time = t_available , geneExpression = as.numeric (as.matrix (table_in[i,]) ) ) )
    (lm ( geneExpression ~ bs (time , df) , data = data) -> m )
    (new <- data.frame(time = setdiff (t_new , t_available) ))
    (predict(m , newdata = new) %>% unname() -> imputed_exp )
    y_ = c(y_ , imputed_exp)
    y[i]
    #
    if (i %% 500 == 0)
    {
      print (i)
      time = seq(0 , 14 , by = 0.05 )
      new <- data.frame( time )
      predict(m , newdata = new) -> spline_model
      plot (spline_model ~ time , pch = '.' , col = "blue" , main = paste (i , paste0 ("day" , idx_miss-1)  , sep = "-") )
      (data = data.frame(time = t , geneExpression = as.numeric (as.matrix (table[i,]) ) ) )
      points (geneExpression[-idx_miss] ~ time[-idx_miss] ,data = data, pch = 16  , col = "green")
      points (geneExpression[idx_miss] ~ time[idx_miss] ,data = data, pch = 16  , col = "red")
    }
  }
  
  (r2 <- 1 - sum ((y_ - y)^2) / sum ((y-mean(y))^2) )
  
  return (r2)
}
pdf(file="../../code/splineDays.pdf")
sapply(c(1:8) , fitSmoothingSpline , table = FAP_damaged_noRep_log_filtered , df = 4) -> day_r2
dev.off()

filterLowExpressedGenes(FAP_damaged_noRep_log , 0.7 , 1) -> FAP_damaged_noRep_log_filtered # 8394 genes 


table = FAP_damaged_noRep_log_filtered 
df = 4


