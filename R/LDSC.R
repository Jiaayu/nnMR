#' @export
nn_LDSC<-function(expo,outcome,an,sample_prev=NA,
                  population_prev=NA,ld,wld,chr_filter=c(1:22),n_blocks=200){
  id.o<-outcome$id.outcome[1]
  id.e<-expo$id.exposure[1]
  expo<-expo%>%mutate(Z=beta.exposure/se.exposure)
  expo<-expo%>%select(SNP=SNP,N=samplesize.exposure,Z=Z
                      ,A1=effect_allele.exposure
                      ,A2=other_allele.exposure)
  expo<-as_tibble(expo)
  outcome<-outcome%>%mutate(Z=beta.outcome/se.outcome)
  outcome<-outcome%>%select(SNP=SNP,N=samplesize.outcome,Z=Z
                            ,A1=effect_allele.outcome
                            ,A2=other_allele.outcome)
  outcome<-as_tibble(outcome)

  dat<-list(expo,outcome)
  names(dat)<-c(id.e,id.o)

  rm(expo,outcome)

  res<-try(ldscr::ldsc_rg(dat,ancestry = an,sample_prev=sample_prev,
                          population_prev=population_prev,ld=ld,wld=wld,
                          n_blocks=n_blocks,chr_filter=chr_filter))
  return(res)
}

