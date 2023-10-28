#' @export
myy_ldsc <- function(file_list,out_dat,outfilename,savepath){
  data <- NULL
  for (file_id in file_list){
    dat <- data.table::fread(paste0('D:/486/原始文件/',file_id,'.txt.gz'),fill = T)
    exp_dat <- TwoSampleMR::format_data(dat = dat,
                                        type = 'exposure',phenotype_col = file_id,
                                        snp_col = 'MarkerName',effect_allele_col = 'Allele1',
                                        other_allele_col = 'Allele2',eaf_col = 'FreqSE',
                                        beta_col = 'Effect',se_col = 'StdErr',pval_col = 'P-value',samplesize_col = 'TotalSampleSize',
                                        chr_col = 'Chr',pos_col = 'Pos')
    exp_dat$samplesize.exposure=7284
    ld=MRmyy::ldsc(
      exp_name=file_id,
      out_name='outcome',
      expfilename = NULL,
      outfilename = NULL,
      exp_dat = exp_dat,
      out_dat = out_dat,
      exp_type = "binary",
      out_type = "binary",
      hm3 = "C:/Users/Zz/Desktop/ldsc/eur_w_ld_chr/w_hm3.snplist",
      ld = "C:/Users/Zz/Desktop/ldsc/eur_w_ld_chr",
      wld = "C:/Users/Zz/Desktop/ldsc/eur_w_ld_chr",
      pop_rate = NULL,
      remove_HLA = FALSE,
      run_MRlap = FALSE,
      MAF = 0.01,
      pval = NA,
      kb = NA,
      r2 = NA,
      MR_reverse = NA,
      save_logfiles = FALSE,
      save_path=savepath)
    if (is.null(data)){
      data <- ld
    }else{
      data <- rbind(data,ld)
    }
    fwrite(data,'LDSC.csv')

  }
}
