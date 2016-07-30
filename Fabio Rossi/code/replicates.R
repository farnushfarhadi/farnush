drawReplicates <- function (table , cellTypeCode , gene)
{
  which(table$tracking_id == gene) -> idx
  
  if(cellType == "EC_WT")
  {
    table[idx , -c(1,2)] -> g  
    plot ( y= g[,c(1,2,5,8,9,10,11,12,14)],x = c(0,2,3,4,5,6,7,10,14) , type = "l")
    points( y= g[,c(2:4)],x = c(2,2,2) )
    points( y= g[,c(5:7)],x = c(3,3,3) )
    points( y= g[,c(12,13)],x = c(10,10) )
  }
}