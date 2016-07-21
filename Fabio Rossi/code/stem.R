#STEM
read.table("../../stem/EC.txt" , header = TRUE , sep= "\t") -> stem_dat
stem_dat[-lowData_idx , ] -> stem_dat
for (i in 1:dim(stem_dat)[2])
{
  as.matrix(stem_dat) -> stem_dat
  (which (stem_dat [ , i] == "FAIL")  -> idx)
  if (length(idx) >0)
    stem_dat [idx , i] <- ""
}
write.table (stem_dat , file = "~/Desktop/internship-ubc/stem/EC_damaged_processed.txt"  , sep = "\t" , quote =FALSE , row.names = FALSE)
