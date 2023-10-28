#' @export

filter_F <- function(dat){
    dat$R2<-(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))/
      ((2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))+
         (2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$se.exposure^2)*dat$samplesize.exposure))
    dat$F<- (dat$samplesize.exposure-2)*dat$R2/(1-dat$R2)
    dat<-subset(dat,F>10)
    return(dat)
}
