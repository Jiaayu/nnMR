setwd("C:/Users/Zz/Desktop/免疫细胞")
setwd('C:/Users/Zz/Desktop/nnMR/nnMR')
library(data.table)
library(devtools)
library(usethis)
a <- fread('731immune_cells_info.csv',header = T)
transform_immune <- a
dict_immune <- a[,c(1,2)]
dict_immune <- setNames( dict_immune$TRAIT,dict_immune$ID)

usethis::use_data(transform_immune, overwrite = TRUE)
usethis::use_data(dict_immune, overwrite = TRUE)
setwd('C:/Users/Zz/Desktop/nnMR/nnMR')
document()
