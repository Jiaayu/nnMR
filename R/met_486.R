#' @export
#' @author Jieyu
#' @param exposure_path 暴露文件所在的目录，运行后会遍历该目录的所有文件，逐个与结局跑MR分析
#' @param outcome_path 结局文件的路径（目前只支持单文件）
#' @param p_threshold 筛选的P值，默认5e-8
#' @param kb Clump时的kb值，默认10000
#' @param r2 Clump时的r2，默认0.001
#' @param outcome_source 结局数据的来源，可选'ieu','finn'，会自动格式化结局
met_486 <- function(exposure_path,exposure_name,outcome_path,outcome_name,
                    p_threshold=5e-8,kb=10000,r2=0.001,outcome_source,result_path){

  if (outcome_source=='ieu'){
    outcome_snp <<- 'variant_id'
    outcome_beta<<-'beta'
    outcome_se<<-'standard_error'
    outcome_effect_allele<<-'effect_allele'
    outcome_other_allele='other_allele'
    outcome_eaf<<-'effect_allele_frequency'
    outcome_pval<<-'p_value'
  }else if (outcome_source=='finn'){
    outcome_snp<<-'rsids'
    outcome_beta<<-'beta'
    outcome_se<<-'sebeta'
    outcome_effect_allele<<-'alt'
    outcome_other_allele<<-'ref'
    outcome_eaf<<-'af_alt'
    outcome_pval<<-'pval'
  }else{
    print(outcome_source)
  }
  # 遍历暴露文件
  start_time <<- proc.time()[3]
  if (!dir.exists('harmonise_file')){
    dir.create('harmonise_file')
  }
  if (!dir.exists('plot')){
    dir.create('plot')
  }

  checkpoint_path<<-paste0('486(',p_threshold,')-',outcome_name,'.RData')
  if (file.exists(checkpoint_path)){
    load(checkpoint_path)
  }else{
    checkpoint <<- c()
  }

  final_path=paste0('486(',p_threshold,')-',outcome_name,'.csv')
  all_files <- list.files(path = exposure_path)
  csv_files <<- all_files[grepl("\\.csv$", all_files)]
  if (outcome_source=='other'){
    outcome_row_data <- outcome_path
  }else{
    outcome_row_data <- data.table::fread(outcome_path,header = T,fill = T)
  }

  total_count <- length(csv_files)
  now_count <<- 1
  skip_next <<- FALSE
  # 循环读取csv
  for (file_name in csv_files){


    file_id <- gsub("\\.csv$", "", file_name)

    # 主程序
    if (!(file_id %in% checkpoint)){
      cat('(',now_count,'/',total_count,')正在分析 ',file_id,'...',sep = '')
        exposure_raw_data <- nnMR::nn_read_exposure_data(
        file=paste0(exposure_path,'\\',file_name),
        sep = ",",
        snp_col = "MarkerName",
        beta_col = "Effect",
        se_col = "StdErr",
        effect_allele_col = "Allele1",
        other_allele_col = "Allele2",
        eaf_col = "Freq1",
        pval_col = "P-value",
        samplesize_col = 'TotalSampleSize',
        chr_col = 'Chr',
        pos_col = 'Pos',
        phenotype_col=file_id,
      )
      exposure_raw_data$exposure <- file_id

      # 筛选P
      exp_dat <- exposure_raw_data[exposure_raw_data$pval.exposure<p_threshold,]

      # SNP过少跳过
      if (nrow(exp_dat)<3){
        cat('SNP不足，已跳过\n')
        checkpoint <<- c(checkpoint,file_id)
        save(checkpoint,file = checkpoint_path)
        now_count <<- now_count+1
        next
      }


      # Clump
      tryCatch({
        exp_dat_clumped <- exp_dat %>%
          rename(rsid = SNP,
                 pval = pval.exposure) %>%
          nnMR::nn_ld_clump(dat = .,
                            clump_kb = kb,
                            clump_r2 = r2,
                            clump_p = 1,
                            plink_bin =plinkbinr::get_plink_exe(),
                            bfile = "D:\\clump_pop\\EUR")%>%
          rename(SNP = rsid,
                 pval.exposure = pval)
      },error=function(error){
        skip_next <<- TRUE
      })
      if (skip_next){
        checkpoint <<- c(checkpoint,file_id)
        save(checkpoint,file = checkpoint_path)
        now_count <<- now_count+1
        cat('SNP不足，已跳过\n')
        next
      }

      # SNP过少保存
      if (nrow(exp_dat_clumped)<3){
        checkpoint <<- c(checkpoint,file_id)
        save(checkpoint,file = checkpoint_path)
        now_count <<- now_count+1
        cat('SNP不足，已跳过\n')
        next
      }

      # 筛选F和R2
      exp_dat_clumped <- filter_F(exp_dat_clumped)
      if (nrow(exp_dat_clumped)<3){
        checkpoint <<- c(checkpoint,file_id)
        save(checkpoint,file = checkpoint_path)
        now_count <<- now_count+1
        cat('SNP不足，已跳过\n')
        next
      }

      # 读取结局
      if (outcome_source=='other'){
        snps <- exp_dat_clumped[['SNP']]
        out_dat <- subset(outcome_row_data,SNP  %in% snps)
      }else{
        out_dat <- nnMR::nn_read_outcome_data(
          snps = exp_dat_clumped$SNP,
          filename = outcome_row_data,
          #sep = "\t",
          snp_col = outcome_snp,
          beta_col = outcome_beta,
          se_col = outcome_se,
          effect_allele_col = outcome_effect_allele,
          other_allele_col = outcome_other_allele,
          eaf_col = outcome_eaf,
          pval_col = outcome_pval
        )
      }
      dat_harmonised <- nnMR::nn_harmonise_data(
        exposure_dat =  exp_dat_clumped,
        outcome_dat = out_dat
      )

      if (nrow(dat_harmonised)<3){
        checkpoint <<- c(checkpoint,file_id)
        save(checkpoint,file = checkpoint_path)
        now_count <<- now_count+1
        cat('SNP不足，已跳过\n')
        next
      }


      res <- TwoSampleMR::mr(dat_harmonised,method_list = c('mr_egger_regression','mr_weighted_median','mr_ivw'))
      reg <- generate_odds_ratios(res)
      if (reg[reg$method=='Inverse variance weighted','pval']>0.05){
        checkpoint <<- c(checkpoint,file_id)
        save(checkpoint,file = checkpoint_path)
        now_count <<- now_count+1
        cat('P值不显著\n')
        cat('P:',reg[reg$method=='Inverse variance weighted','pval'],sep = '')
        next
      }

      if (reg[reg$method=='MR Egger','b']*reg[reg$method=='Weighted median','b']>0&
          reg[reg$method=='MR Egger','b']*reg[reg$method=='Inverse variance weighted','b']>0&
          reg[reg$method=='Weighted median','b']*reg[reg$method=='Inverse variance weighted','b']>0){

        write.csv(dat_harmonised,paste0('harmonise_file/',file_id,'_harmonised.csv'))
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

        data <- merge(reg,presso_dat,by='exposure')[c(1:3),-c(2,3)]
        data$outcome <- outcome_name

        # dat_harmonised <- data.table::fread('C:/Users/Zz/Desktop/gastroduodenal/remake/harmonise_file/M01110_harmonised.csv')
        ple <- mr_pleiotropy_test(dat_harmonised)[,-c(1:3,6)]
        ple$exposure <- file_id
        colnames(ple)[3] <- 'plei_pval'


        ht <- mr_heterogeneity(dat_harmonised)[2,-c(1:3,5)]
        ht$i <- (ht$Q-ht$Q_df)/ht$Q
        ht$exposure <- file_id


        mgx <- merge(ple,ht,by='exposure')
        new <- merge(data,mgx,by='exposure')
        data <- new



        if (!file.exists(final_path)){
          data.table::fwrite(data,final_path)
        }else{
          existing_data <- data.table::fread(final_path)
          new_data <- rbind(existing_data,data)
          data.table::fwrite(new_data,final_path)
        }

        res_single <- TwoSampleMR::mr_singlesnp(dat_harmonised)
        res_loo <- TwoSampleMR::mr_leaveoneout(dat_harmonised)
        p1 <- TwoSampleMR::mr_scatter_plot(res, dat_harmonised)

        ggplot2::ggsave(p1[[1]], file=paste0('./plot/',file_id,"_res.pdf"), width=7, height=7)


        #森林图
        p2 <- TwoSampleMR::mr_forest_plot(res_single)
        ##保存图片
        ggplot2::ggsave(p2[[1]], file=paste0('./plot/',file_id,"_forest.pdf"), width=7, height=7)


        ##留一法图
        p3 <- TwoSampleMR::mr_leaveoneout_plot(res_loo)
        ##保存图片
        ggplot2::ggsave(p3[[1]], file=paste0('./plot/',file_id,"_loo.pdf"), width=7, height=7)


        ##漏斗图
        res_single <- TwoSampleMR::mr_singlesnp(dat_harmonised)
        p4 <- TwoSampleMR::mr_funnel_plot(res_single)
        ##保存图片
        ggplot2::ggsave(p4[[1]], file=paste0('./plot/',file_id,"_single.pdf"), width=7, height=7)

      }
      else{
        cat('(',now_count,'/',total_count,') ',file_id,'...已跳过\n',sep = '')
      }
      now_count <<- now_count+1

    }else{
      cat('不同号\n')
      now_count<<-now_count+1
      next
      }

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

