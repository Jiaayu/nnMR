#' @export
met_ice <- function(exposure_path,outcome_path,
                    p_threshold=5e-8,kb=10000,r2=0.001,outcome_source){
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
    stop('error')
  }
  # 遍历暴露文件
  start_time <<- proc.time()[3]
  checkpoint_path<<-paste0('486(',p_threshold,')-',outcome_source,'.RData')
  if (file.exists(checkpoint_path)){
    load(checkpoint_path)
  }else{
    checkpoint <<- c()
  }

  final_path=paste0('486(',p_threshold,')-',outcome_source,'.csv')
  all_files <- list.files(path = exposure_path)
  csv_files <<- all_files[grepl("\\.csv$", all_files)]

  outcome_row_data <- data.table::fread(outcome_path,header = T,fill = T)
  total_count <- length(csv_files)
  now_count <<- 1
  skip_next <<- FALSE
  # 循环读取csv
  for (file_name in csv_files){
    file_id <- gsub("\\.csv$", "", file_name)
    checkpoint <- read.csv(final_path)[[2]]
    # 主程序
    if (!(file_id %in% checkpoint)){
      cat('(',now_count,'/',total_count,')正在分析 ',file_id,'...',sep = '')
      exposure_raw_data <- nn_read_exposure_data(
        file=paste0(exposure_path,'\\',file_name),
        sep = ",",
        snp_col = "rsids",
        beta_col = "Beta",
        se_col = "SE",
        effect_allele_col = "effectAllele",
        other_allele_col = "otherAllele",
        pval_col = "Pval",
        phenotype_col='exposure'
      )

      # 筛选P
      exp_dat <- exposure_raw_data[exposure_raw_data$pval.exposure<p_threshold,]

      # SNP过少跳过
      if (nrow(exp_dat)<3){
        data_row <- data.frame('exposure'=file_id,'outcome'=outcome_source,'p_threshold'=p_threshold,'kb'=kb,'r2'=r2,
                               'NO.SNP'=nrow(exp_dat_clumped),'beta'=NA,
                               'p.value'=NA,'or.95CI'=NA)
        if (!file.exists(final_path)){
          write.csv(data_row,final_path)
        }else{
          existing_data <- read.csv(final_path)[,-1]
          new_data <- rbind(existing_data,data_row)
          write.csv(new_data,final_path)
        }
        checkpoint <<- c(checkpoint,file_id)
        save(checkpoint,file = checkpoint_path)
        now_count <<- now_count+1
        cat('SNP不足，已跳过\n')
        next
      }


      # Clump
      tryCatch({exp_dat_clumped <- exp_dat %>%
        rename(rsid = SNP,
               pval = pval.exposure) %>%
        nn_ld_clump(dat = .,
                          clump_kb = kb,
                          clump_r2 = r2,
                          clump_p = 1,
                          plink_bin =plinkbinr::get_plink_exe(),
                          bfile = "D:\\clump_pop\\EUR")%>%
        rename(SNP = rsid,
               pval.exposure = pval)},
        error=function(e){
          print('error,skipping')
          data_row <- data.frame('exposure'=file_id,'outcome'=outcome_source,'p_threshold'=p_threshold,'kb'=kb,'r2'=r2,
                                 'NO.SNP'=NA,'beta'=NA,
                                 'p.value'=NA,'or.95CI'=NA)
          if (!file.exists(final_path)){
            write.csv(data_row,final_path)
          }else{
            existing_data <- read.csv(final_path)[,-1]
            new_data <- rbind(existing_data,data_row)
            write.csv(new_data,final_path)
          }
          checkpoint <<- c(checkpoint,file_id)
          save(checkpoint,file = checkpoint_path)
          now_count <<- now_count+1
          cat('SNP不足，已跳过\n')
          skip_next <<- TRUE
        }



      )
      if (skip_next) {
        skip_next <<- FALSE
        next
      }

      # SNP过少保存
      if (nrow(exp_dat_clumped)<3){
        data_row <- data.frame('exposure'=file_id,'outcome'=outcome_source,'p_threshold'=p_threshold,'kb'=kb,'r2'=r2,
                               'NO.SNP'=nrow(exp_dat_clumped),'beta'=NA,
                               'p.value'=NA,'or.95CI'=NA)
        if (!file.exists(final_path)){
          write.csv(data_row,final_path)
        }else{
          existing_data <- read.csv(final_path)[,-1]
          new_data <- rbind(existing_data,data_row)
          write.csv(new_data,final_path)
        }
        checkpoint <<- c(checkpoint,file_id)
        save(checkpoint,file = checkpoint_path)
        now_count <<- now_count+1
        cat('SNP不足，已跳过\n')
        next
      }

      # 筛选F
      exp_dat_clumped['F'] <- (exp_dat_clumped['beta.exposure']^2) / (exp_dat_clumped['se.exposure']^2)
      exp_dat_clumped <- exp_dat_clumped[exp_dat_clumped$F >= 10,]

      # 读取结局
      out_dat <- nn_read_outcome_data(
        snps = exp_dat_clumped$SNP,
        filename = outcome_row_data,
        #sep = "\t",
        snp_col = outcome_snp,
        beta_col = outcome_beta,
        se_col = outcome_se,
        effect_allele_col = outcome_effect_allele,
        other_allele_col = outcome_other_allele,
        eaf_col = outcome_eaf,
        pval_col = outcome_pval,
        phenotype_col = 'outcome'
      )
      dat_harmonised <- nn_harmonise_data(
        exposure_dat =  exp_dat_clumped,
        outcome_dat = out_dat
      )
      if (nrow(dat_harmonised)<3){
        data_row <- data.frame('exposure'=file_id,'outcome'=outcome_source,'p_threshold'=p_threshold,'kb'=kb,'r2'=r2,
                               'NO.SNP'=nrow(dat_harmonised),'beta'=NA,
                               'p.value'=NA,'or.95CI'=NA)
        if (!file.exists(final_path)){
          write.csv(data_row,final_path)
        }else{
          existing_data <- read.csv(final_path)[,-1]
          new_data <- rbind(existing_data,data_row)
          write.csv(new_data,final_path)
        }
        checkpoint <<- c(checkpoint,file_id)
        save(checkpoint,file = checkpoint_path)
        now_count <<- now_count+1
        cat('SNP不足，已跳过\n')
        next
      }

      #write.csv(dat_harmonised,paste0(file_id,'_harmonised.csv'))
      res <- nn_mr(dat_harmonised)
      reg <- generate_odds_ratios(res)
      if (res$nsnp[1]<3){
        data_row <- data.frame('exposure'=file_id,'outcome'=outcome_source,'p_threshold'=p_threshold,'kb'=kb,'r2'=r2,
                               'NO.SNP'=nrow(dat_harmonised),'beta'=NA,
                               'p.value'=NA,'or.95CI'=NA)
        if (!file.exists(final_path)){
          write.csv(data_row,final_path)
        }else{
          existing_data <- read.csv(final_path)[,-1]
          new_data <- rbind(existing_data,data_row)
          write.csv(new_data,final_path)
        }
        checkpoint <<- c(checkpoint,file_id)
        save(checkpoint,file = checkpoint_path)
        now_count <<- now_count+1
        cat('SNP不足，已跳过\n')
        next
      }
      s <- reg[reg$method=='Inverse variance weighted','nsnp']
      b <- reg[reg$method=='Inverse variance weighted','b']
      p <- reg[reg$method=='Inverse variance weighted','pval']
      o <- reg[reg$method=='Inverse variance weighted','or']
      lo <- reg[reg$method=='Inverse variance weighted','or_lci95']
      uo <- reg[reg$method=='Inverse variance weighted','or_uci95']
      orr <- paste0(round(o,3),'(',round(lo,3),'-',round(uo,3),')')
      if (p<=0.05){
        cat('P-value：',round(p,3),'*\n',sep = '')
      }else{
        cat('P-value：',round(p,3),'\n',sep = '')
      }
      data_row <- data.frame('exposure'=file_id,'outcome'=outcome_source,'p_threshold'=p_threshold,'kb'=kb,'r2'=r2,'NO.SNP'=s,'beta'=b,
                             'p.value'=p,'or.95CI'=orr)
      if (!file.exists(final_path)){
        write.csv(data_row,final_path)
      }else{
        existing_data <- read.csv(final_path)[,-1]
        new_data <- rbind(existing_data,data_row)
        write.csv(new_data,final_path)
      }
      checkpoint <<- c(checkpoint,file_id)
      save(checkpoint,file = checkpoint_path)



    }
    else{
      cat('(',now_count,'/',total_count,') ',file_id,'...已跳过\n',sep = '')
    }
    now_count <<- now_count+1



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

