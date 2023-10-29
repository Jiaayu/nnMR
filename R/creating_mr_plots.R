#' @param dat_harmonised harmonise后的数据
#'
#' @param res 初次MR分析后的数据
#' @param file_id file_id
#'
#' @export
#' @title Auto creating MR plots

creating_mr_plots <- function(dat_harmonised,res,file_id){
  if (!dir.exists('plot')){
    dir.create('plot')
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
