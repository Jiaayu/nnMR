#' @export
nn_LDSC_auto <- function(exposure_path,exposure_source,result_path,
                         outcome_path,outcome_source,outcome_case){
  # setwd('C:\\Users\\Zz\\Desktop\\饮食-代谢-ibs\\finnIBS')
  outcome_name <- 'Finn-IBS'


  # 加载存档
  if(file.exists(result_path)){
    checkpoint <- as.vector(read.csv(result_path)[['h2.trait']])
  }else{
    checkpoint <- c()
  }


  # 读取结局数据
  outcome_data<-data.table::fread(file = outcome_path,header = T,showProgress=F)
  if (outcome_source=='finn'){
    outcome_data <- outcome_data[,-c(6,8,12:13)]
    outcome_data$N <- outcome_case
    colnames(outcome_data) <- c("chr", "pos","ref","alt","rsid","pval","beta","se",
                                "eaf","N")
  }

  if (outcome_source=='ieu'){
    outcome_data$N <- outcome_data$N_CASE+outcome_data$N_CONTROL
    colnames(outcome_data) <- c("rsid", "pval","chr","pos","alt","ref","eaf","beta",
                                "se","N1","N2","N")
  }


  # 数据处理代码...
  outcome_data <- nnMR::nn_read_outcome_data(filename = outcome_data,
                                    sep = ',',
                                    snp_col = 'rsid',
                                    effect_allele_col = 'alt',
                                    other_allele_col = 'ref',
                                    eaf_col = 'eaf',
                                    beta_col = 'beta',
                                    se_col = 'se',
                                    pval_col = 'pval',
                                    samplesize_col = 'N',
                                    chr_col = 'chr',
                                    pos_col = 'pos')

  file_list <- list.files(exposure_path)
  total_num <- length(file_list)

  for (i in 1:length(file_list)){
    file_id <- file_list[i]

    cat('(',i,'/',total_num,') 当前分析：',file_id,sep='')
    if (file_id %in% checkpoint){
      cat('...已分析\n')
      next
    }
    # exposure_path <- paste0('D:\\486\\原始文件\\',file_id,'.txt.gz')
    # 读取暴露数据
    cat('...读取')
    exposure_data <- data.table::fread(paste0(exposure_path,'/',file_id),fill = T,showProgress=F)
    if (exposure_source==1400){
      dd <- readxl::read_xlsx('D:\\1400\\1400血清代谢物欧洲人群总表.xlsx')[,c(1,4)]
      dd$study_ID <- paste0(dd$study_ID,'.tsv.gz')
      dd$initial_sample_size <- gsub(" European", "", dd$initial_sample_size)
      size <- as.numeric(dd[dd$study_ID==file_id,'initial_sample_size'])
      exposure_data$N <- size
      colnames(exposure_data) <- c("chr", "pos","alt","ref","eaf","beta","se","pval","rsid","N")
    }else if (exposure_source==486){
      exposure_data <- exposure_data[,-c(5,6,7,11:15,19:34)]
      colnames(exposure_data) <- c("rsid", "alt","ref","eaf","beta","se","pval","N",
                                   "chr","pos")
    }else if(exposure_source=='rewrite')
      colnames(exposure_data) <- c("rsid", "alt","ref","eaf","beta","se","pval","N",
                                 "chr1","pos1","chr","pos")

    exposure_data <- nn_read_exposure_data(filename = exposure_data,
                                           sep = ',',
                                           snp_col = 'rsid',
                                           effect_allele_col = 'alt',
                                           other_allele_col = 'ref',
                                           eaf_col = 'eaf',
                                           beta_col = 'beta',
                                           se_col = 'se',
                                           pval_col = 'pval',
                                           samplesize_col = 'N',
                                           chr_col = 'chr',
                                           pos_col = 'pos')

    # LDSC分析
    cat('...LDSC')
    res <- nn_LDSC(exposure_data,outcome_data,an='EUR',
                   ld='C:\\Users\\Zz\\Desktop\\lddsc\\eur',
                   wld='C:\\Users\\Zz\\Desktop\\lddsc\\eur')
    if (!is.list(res)){
      new_row <- as.data.frame(list('h2.trait'=file_id,'h2.mean_chisq'=NA,
                               'h2.lambda_gc'=NA,'h2.intercept'=NA,'h2.intercept_se'=NA,
                               'h2.ratio'=NA,'h2.ratio_se'=NA,'h2.h2_observed'=NA,'h2.h2_observed_se'=NA,
                               'h2.h2_Z'=NA,'h2.h2_p'=NA,'rg.rg'=NA,'rg.rg_se'=NA,'rg.rg_p'=NA))
    cat('...P-value:NA\n')
    }else{
      h2 <- as.data.frame(res[1])[1,]
      rg <- as.data.frame(res[2])[,-2]
      new_row <- merge(h2,rg,by.x = 'h2.trait',by.y = 'rg.trait1')
      new_row$h2.trait <- file_id
      cat('...P-value:',round(rg$rg.rg_p[1],3),'\n',sep = '')
    }
    if (file.exists(result_path)){
      existing_data <- read.csv(result_path,header = T)[,-1]
      new_data <- rbind(existing_data,new_row)
      write.csv(new_data,result_path)
    }else{
      write.csv(new_row,result_path)
    }

  }
}


