library(nnMR)
# 保存代码
library(devtools)
setwd('C:\\Users\\Zz\\Desktop\\nnMR\\nnMR')
document()



library(profvis)
setwd('C:\\Users\\Zz\\Desktop\\gastroduodenal')
profvis({
  nn_LDSC_auto(exposure_path='D:\\1400\\原文件\\3',exposure_source=1400,result_path='LDSC3.csv',
             outcome_path='finngen_R9_K11_GASTRODUOULC.gz',outcome_source='finn',outcome_case=329603)
  })



Sys.setenv(CLI_PROGRESS_DISABLE = TRUE)
nn_LDSC_auto(exposure_path='D:\\1400\\原文件\\3',exposure_source=1400,result_path='LDSC3.csv',
             outcome_path='finngen_R9_K11_GASTRODUOULC.gz',outcome_source='finn',outcome_case=329603)



   library(readxl)
library(writexl)
dt1 <- read.csv('LDSC_data.csv')[,-1]
dt2 <- read_xlsx('1400(5e-06)-finnUlcer-translated.xlsx')
dt <- merge(dt1,dt2,by.x = 'h2.trait',by.y = 'exposure',all.x = T)
write_xlsx(dt,'LDSC_data_filtered.xlsx')



dd <- readxl::read_xlsx('D:\\1400\\1400血清代谢物欧洲人群总表.xlsx')[,c(1,4)]
dd$study_ID <- paste0(dd$study_ID,'.tsv.gz')
dd$initial_sample_size <- gsub(" European", "", dd$initial_sample_size)
size <- as.numeric(dd[dd$study_ID=='GCST90200663.tsv.gz','initial_sample_size'])
class(size)


