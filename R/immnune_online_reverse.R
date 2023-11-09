#' @param exposure_data 必须为数据框格式的暴露数据
#' @param output_path 输出结果文件的名字，默认为result_reverse.csv，如果使用多线程分析则每个线程输出文件名字必须不同（否则会覆盖）
#' @param multiprocess 是否启用多线程，默认关闭
#' @param start_num 多线程的分析起点和终点，如果要使用多线程，则分别设置每个线程的起点和终点，比如线程1设置为1和350，线程2设置为351和731（总共731）
#' @param end_num 多线程的分析起点和终点，如果要使用多线程，则分别设置每个线程的起点和终点，比如线程1设置为1和350，线程2设置为351和731（总共731）
#' @param exposure_pval 筛选暴露的P值，默认为1e-5
#' @param clump_r2 clump使用的r2和kb，默认为0.1和500
#' @param clump_kb clump使用的r2和kb，默认为0.1和500
#' @title Auto run 731 immune cells MR
#' @examples
#' suppressMessages(immnune_online_reverse(exposure_data,output_path='result_reverse.csv',
#'                                         multiprocess=FALSE,start_num=NA,end_num=NA,
#'                                         exposure_pval=5e-8,clump_r2=0.001,clump_kb=10000'))
#'
#'
#' @export
#'
#'

# setwd('C:/Users/Zz/Desktop/test_immune')
# data <- data.table::fread('finngen_R9_K11_GASTRODUOULC.gz')
#
# out_dat <- data
# out_dat$samplesize <- 350062
# out_dat$phenotype <- 'gastro'
# out_dat <- TwoSampleMR::format_data(dat = out_dat,type = 'exposure',
#                                     snp_col = 'rsids',chr_col = '#chrom',pos_col = 'pos',
#                                     effect_allele_col = 'alt',other_allele_col = 'ref',
#                                     pval_col = 'pval',beta_col = 'beta',se_col = 'sebeta',
#                                     eaf_col = 'af_alt',samplesize_col = 'samplesize',phenotype_col = 'phenotype')
# exposure_data <- out_dat
# library(nnMR)
# exposure_pval=5e-8
# clump_r2=0.001
# clump_kb=10000
# output_path <- 'result.csv'
immnune_online_reverse <- function(exposure_data,output_path='result_reverse.csv',
                         multiprocess=FALSE,start_num=NA,end_num=NA,
                         exposure_pval=5e-8,clump_r2=0.001,clump_kb=10000){

  # 待传入
  start_time <- proc.time()[3]
  exposure_name <- exposure_data$exposure[1]

  # 设置前置数据
  data("dict_immune")
  data("transform_immune")
  file_list <- transform_immune[[1]]
  if (multiprocess){
    file_list <- file_list[start_num:end_num]
  }
  # 记录进度
  z <- 1
  zz <- length(file_list)

  # 读取存档
  if (file.exists(output_path)){
    save_id <<- fread(output_path)[[1]]
  }else{
    save_id <<- NULL
  }



  #==========================
  # 处理暴露
  # 读取暴露并筛选
  exp_dat <- subset(exposure_data,pval.exposure<exposure_pval)
  exp_dat <- nnMR::filter_F(exp_dat)
  # 在线CLUMP
  while (TRUE) {
    tryCatch({
      exp_dat_clump <- TwoSampleMR::clump_data(exp_dat,clump_kb = clump_kb,clump_r2 = clump_r2)
      break
    },error=function(error){
      cat('error 502, retrying...\n')
    })
  }
  cat('CLUMP后还剩',nrow(exp_dat_clump),'个SNP\n')
  # 添加exposure名字
  exp_dat_clump$exposure <- exposure_name

  # 开始循环
  for (file_id in file_list){
    skip <- FALSE
    # file_id <- file_list[1]
    # file_id <- 'GCST90001500'
    outcome_name <- dict_immune[[file_id]]
    cat('(',z,'/',zz,') Analysing',' \'',exposure_name,'\'',' to \'',outcome_name,'\'',sep = '')
    # 判断是否已经做过
    if (file_id %in% save_id){
      z <- z+1
      cat('...已跳过\n')
      next
    }



    #==========================
    # 处理结局
    file_id_gwas <- paste0('ebi-a-',file_id)
    # 在线提取结局
    while (TRUE) {
      tryCatch({
        out_dat <- TwoSampleMR::extract_outcome_data(snps = exp_dat_clump$SNP,outcomes = file_id_gwas)
        break
      },error=function(error){
        cat('error 502, retrying...\n')
      })
    }

    if (nrow(out_dat)<3){
      cat('结局snp不足，已跳过\n')
      z <- z+1
      next
    }



    # harmonise
    dat_harmonised <- TwoSampleMR::harmonise_data(exposure_dat = exp_dat_clump,
                                                  outcome_dat = out_dat)
    # 如果harmonise后SNP少于3个，则跳过
    if (nrow(dat_harmonised)<3){
      cat('harmonise后SNP不足，已跳过\n')
      z <- z+1
      next
    }
    # MR分析
    res <- TwoSampleMR::mr(dat_harmonised,method_list = c('mr_egger_regression','mr_weighted_median','mr_ivw'))
    reg <- TwoSampleMR::generate_odds_ratios(res)
    cat('...P-value:',reg$pval[3])
    # 如果P不显著则跳过
    if (reg$pval[3]>0.05){
      cat('...P值不显著，已跳过\n')
      z <- z+1
      next
    }
    cat('***','进行后续分析...')
    # 创建文件夹，写入harmonise文件
    if (!dir.exists('harmonise_reverse')){
      dir.create('harmonise_reverse')
    }
    data.table::fwrite(dat_harmonised,paste0('harmonise_reverse/',file_id,'_harmonised_reverse.csv'))

    # PRESSO分析
    tryCatch({
      presso <<- TwoSampleMR::run_mr_presso(dat_harmonised,NbDistribution = 5000)
      outlier_corrected_p <<- presso[[1]]$`Main MR results`$`P-value`[2]
      global_p <<- presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
      outliers <<- presso[[1]]$`MR-PRESSO results`$`Distortion Test`$`outliers Indices`
      if (is.null(outliers)){
        presso_dat <<- data.frame(list('exposure'=file_id,'PRESSO_outlier_corrected_p'=outlier_corrected_p,
                                       'PRESSO_global_p'=global_p,'PRESSO_outliers'=NA))
      }else{
        presso_dat <<- data.frame(list('exposure'=file_id,'PRESSO_outlier_corrected_p'=outlier_corrected_p,
                                       'PRESSO_global_p'=global_p,'PRESSO_outliers'=outliers))
      }
    },error=function(error){
      presso_dat <<- data.frame(list('exposure'=file_id,'PRESSO_outlier_corrected_p'=NA,
                                     'PRESSO_global_p'=NA,'PRESSO_outliers'=NA))
    })
    # 整合数据
    reg$exposure <- file_id
    data <- merge(reg,presso_dat,by ='exposure')[c(1:3),]
    ple <- TwoSampleMR::mr_pleiotropy_test(dat_harmonised)[,-c(1:3,6)]
    ple$exposure <- file_id
    colnames(ple)[3] <- 'pleiotropy_pval'


    ht <- TwoSampleMR::mr_heterogeneity(dat_harmonised)[2,-c(1:3,5)]
    ii <- (ht$Q-ht$Q_df)/ht$Q
    ht$exposure <- file_id
    if (ii<0){
      ht$I2 <- 0
    }else{
      ht$I2 <- (ht$Q-ht$Q_df)/ht$Q
    }


    mgx <- merge(ple,ht,by='exposure')
    new <- merge(data,mgx,by='exposure')
    data <- new[,-c(2,3)]

    if (!file.exists(output_path)){
      data.table::fwrite(data,output_path)
    }else{
      existing_data <- data.table::fread(output_path)
      new_data <- rbind(existing_data,data)
      data.table::fwrite(new_data,output_path)
    }

    creating_mr_plots(dat_harmonised=dat_harmonised,res=res)

    z <- z+1
    cat('ok\n')
  }
  # 统计时间
  end_time <- proc.time()[3][[1]]-start_time[[1]]
  t <- round(end_time,1)
  if (t<=60){
    cat('分析结束，用时',t,'秒',sep = '')
  }else{
    tm <- floor(t/60)
    ts <- t-60*tm
    cat('分析结束，用时',tm,'分钟',ts,'秒',sep = '')
    rm(list = ls())
  }
}


