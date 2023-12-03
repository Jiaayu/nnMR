#' @param exposure_path 731种免疫细胞的csv文件所在目录
#' @param outcome_data 必须为数据框格式的结局数据
#' @param output_path 输出结果文件的名字，默认为result.csv，如果使用多线程分析则每个线程输出文件名字必须不同（否则会覆盖）
#' @param multiprocess 是否启用多线程，默认关闭
#' @param start_num 多线程的分析起点和终点，如果要使用多线程，则分别设置每个线程的起点和终点，比如线程1设置为1和350，线程2设置为351和731（总共731）
#' @param end_num 多线程的分析起点和终点，如果要使用多线程，则分别设置每个线程的起点和终点，比如线程1设置为1和350，线程2设置为351和731（总共731）
#' @param exposure_pval 筛选暴露的P值，默认为1e-5
#' @param clump_r2 clump使用的r2和kb，默认为0.1和500
#' @param clump_kb clump使用的r2和kb，默认为0.1和500
#' @param bfile_path clump所需要的文件路径
#' @title Auto run 731 immune cells MR
#' @examples
#' suppressMessages(immnune_auto(exposure_path='D:/immune/csv',
#'                  outcome_data=outcome_data,output_path='result.csv',
#'                  multiprocess=FALSE,start_num=NA,end_num=NA,
#'                  exposure_pval=1e-5,clump_r2=0.1,clump_kb=500,
#'                  bfile_path='D:/clump_pop/EUR'))
#'
#'
#' @export
immnune_auto_reverse <- function(outcome_path='D:/immune/csv',exposure_data,output_path='result_reverse.csv',
                         multiprocess=FALSE,start_num=NA,end_num=NA,
                         exposure_pval=1e-5,clump_r2=0.1,clump_kb=500,bfile_path='D:/clump_pop/EUR'){
  # test
  # setwd('C:/Users/Zz/Desktop/nntest')
  # exposure_path='D:/immune/exposure'
  #
  # out_dat <- data.table::fread('finngen_R9_K11_GASTRODUOULC.gz')
  # out_dat$samplesize <- 350062
  # out_dat$phenotype <- 'GASTRODUOULC'
  # out_dat <- TwoSampleMR::format_data(dat = out_dat,type = 'outcome',
  #                                     snp_col = 'rsids',chr_col = '#chrom',pos_col = 'pos',
  #                                     effect_allele_col = 'alt',other_allele_col = 'ref',
  #                                     pval_col = 'pval',beta_col = 'beta',se_col = 'sebeta',
  #                                     eaf_col = 'af_alt',samplesize_col = 'samplesize',phenotype_col = 'phenotype')
  #
  # outcome_data <- out_dat
  # output_path='result.csv'
  # exposure_pval=1e-5
  # clump_r2=0.1
  # clump_kb=500
  # bfile_path='D:/clump_pop/EUR'
  # library(nnMR)


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
    save_id <<- fread(output_path)[['file_id']]
  }else{
    save_id <<- NULL
  }


  # 处理暴露
  # exposure_data
  exp_dat <- exposure_data
  exp_dat <- subset(exp_dat,pval.exposure<exposure_pval)
  # clump,如果报错则跳过
  tryCatch({
    exp_dat_clump <- exp_dat %>%
      rename(rsid = SNP,
             pval = pval.exposure) %>%
      nnMR::nn_ld_clump(dat = .,
                        clump_kb = clump_kb,
                        clump_r2 = clump_r2,
                        clump_p = 1,
                        plink_bin =plinkbinr::get_plink_exe(),
                        bfile = bfile_path)%>%
      rename(SNP = rsid,
             pval.exposure = pval)
  },error=function(error){
    stop('Clump出错，重试或调整参数\n')
  })
  # 计算F
  exp_dat_clump <- nnMR::filter_F(exp_dat_clump)

  # 开始循环
  for (file_id in file_list){
    skip <- FALSE
    # file_id <- file_list[13]
    # file_id <- 'GCST90001500'
    outcome_name <- dict_immune[[file_id]]
    cat('(',z,'/',zz,') Analysing',' \'',exposure_name,'\'',' to \'',outcome_name,'\'',sep = '')
    # 判断是否已经做过
    if (file_id %in% save_id){
      z <- z+1
      cat('...已跳过\n')
      next
    }
    # 读取结局并筛选
    outcome_data <- data.table::fread(file.path(outcome_path,paste0(file_id,'_outcome.tsv.gz')))

    # 读取结局并提取相关数据
    out_dat <- subset(outcome_data,SNP %in% exp_dat_clump[['SNP']])

    # harmonise
    dat_harmonised <- TwoSampleMR::harmonise_data(exposure_dat = exp_dat_clump,
                                                  outcome_dat = out_dat)
    # 如果harmonise后SNP少于3个，则跳过
    if (nrow(dat_harmonised)<3){
      cat('...harmonise后SNP不足\n')
      skip_data <- data.frame(list('exposure'=exposure_name,'outcome'=outcome_name,'method'=NA,'nsnp'=nrow(dat_harmonised),'b'=NA,
                                   'se'=NA,'pval'=NA,'lo_ci'=NA,'up_ci'=NA,'or'=NA,'or_lci95'=NA,
                                   'or_uci95'=NA,'file_id'=file_id,'egger_intercept'=NA,'pleiotropy_pval'=NA,
                                   'Q'=NA,'Q_df'=NA,'Q_pval'=NA,'I2'=NA))
      if (!file.exists(output_path)){
        data.table::fwrite(skip_data,output_path)
      }else{
        existing_data <- data.table::fread(output_path)
        new_data <- rbind(existing_data,skip_data)
        data.table::fwrite(new_data,output_path)
      }
      z <- z+1
      next
    }
    # MR分析
    res <- TwoSampleMR::mr(dat_harmonised,method_list = c('mr_egger_regression','mr_weighted_median','mr_ivw'))
    reg <- TwoSampleMR::generate_odds_ratios(res)
    cat('...P-value:',reg$pval[3])
    # 如果P不显著则跳过
    if (reg$pval[3]>0.05){
      cat('...P值不显著...')
    }
    # cat('***','进行后续分析...')
    # 创建文件夹，写入harmonise文件
    if (!dir.exists('harmonise')){
      dir.create('harmonise')
    }
    data.table::fwrite(dat_harmonised,paste0('harmonise/',file_id,'_harmonised.csv'))

    # PRESSO分析
    # tryCatch({
    #   presso <<- TwoSampleMR::run_mr_presso(dat_harmonised,NbDistribution = 5000)
    #   outlier_corrected_p <<- presso[[1]]$`Main MR results`$`P-value`[2]
    #   global_p <<- presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
    #   outliers <<- presso[[1]]$`MR-PRESSO results`$`Distortion Test`$`outliers Indices`
    #   if (is.null(outliers)){
    #     presso_dat <<- data.frame(list('exposure'=file_id,'PRESSO_outlier_corrected_p'=outlier_corrected_p,
    #                                    'PRESSO_global_p'=global_p,'PRESSO_outliers'=NA))
    #   }else{
    #     presso_dat <<- data.frame(list('exposure'=file_id,'PRESSO_outlier_corrected_p'=outlier_corrected_p,
    #                                    'PRESSO_global_p'=global_p,'PRESSO_outliers'=outliers))
    #   }
    # },error=function(error){
    #   presso_dat <<- data.frame(list('exposure'=file_id,'PRESSO_outlier_corrected_p'=NA,
    #                                  'PRESSO_global_p'=NA,'PRESSO_outliers'=NA))
    # })
    # 整合数据
    # data <- merge(reg,presso_dat,by ='exposure')[c(1:3),]
    data <- reg
    data$file_id <- file_id
    ple <- TwoSampleMR::mr_pleiotropy_test(dat_harmonised)[,-c(1:3,6)]
    colnames(ple)[3] <- 'pleiotropy_pval'


    ht <- TwoSampleMR::mr_heterogeneity(dat_harmonised)[2,-c(1:3,5)]
    ii <- (ht$Q-ht$Q_df)/ht$Q
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

    nnMR::creating_mr_plots(dat_harmonised=dat_harmonised,res=res,file_id = file_id)

    z <- z+1
    cat('分析完成\n')
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
    cat('进行FDR校正...')

    # output_path <- 'C:/Users/Zz/Desktop/MR数据/待完成/nntest/result.csv'
    nnMR::nnFDR(output_path)

    cat('完成')
    rm(list = ls())
  }
}
