library(metadynminer)

#weightboltzmann(cvfile, npoints=60, maxfe=100, temp=300, eunits="kJ/mol",
#                  imin=1, imax=NULL, xlim=NULL, ylim=NULL)
#reweightbonomi(cvfile, npoints=60, maxfe=100, temp=300, eunits="kJ/mol",
#               imin=1, imax=NULL, xlim=NULL, ylim=NULL)
#reweighttiwary(cvfile, hillsfile, npoints=60, maxfe=100, nfes=100, temp=300,
#               gamma=10, eunits="kJ/mol", imin=1, imax=NULL, xlim=NULL, ylim=NULL, usefes2=F)

hillsf<-read.hills("../data/HILLS", per=c(T,T))
colvar<-read.table("../data/COLVAR")
nfes <- 100
nbins <- 60
temp <- 300
gamma <- 10
ebtac <- c()
framestosum <- hillsf$size[1]/nfes
suppressWarnings(
  tfes <- fes(hillsf, imax=1, npoints=nbins)-fes(hillsf, imax=1, npoints=nbins)
)
for(i in 1:nfes) {
  tfes <- tfes + fes(hillsf, imin=(i-1)*framestosum+1, imax=i*framestosum, npoints=nbins)
  s1 <- sum(exp(-1000*tfes$fes/8.314/temp))
  s2 <- sum(exp(-1000*tfes$fes/8.314/temp/gamma))
  ebtac<-c(ebtac,s1/s2)
}
icv1 <- ceiling((colvar[,2]+pi)*nbins/2/pi)
icv2 <- ceiling((colvar[,3]+pi)*nbins/2/pi)
bp <- colvar[,4]
ebtacc <- rep(ebtac, each=nrow(colvar)/nfes)
if(length(ebtacc)>nrow(colvar)) ebtacc<-ebtacc[1:nrow(colvar)]
if(length(ebtacc)<nrow(colvar)) ebtacc[length(ebtacc):nrow(colvar)]<-ebtacc[length(ebtacc)]

ofes <- matrix(rep(0, nbins*nbins), nrow=nbins)
for(i in 1:nbins) {
  for(j in 1:nbins) {
    ofes[i,j]<-sum((icv1==i)*(icv2==j)*exp(1000*bp/8.314/temp)/ebtacc)
  }
}
ofes <- -8.314*temp*log(ofes)/1000
ofes <- ofes - min(ofes)
ofes[ofes==Inf]<-100
png("fes.png")
tfes <- fes(hillsf, npoints=nbins)
tfes$fes <- ofes
plot(tfes)
dev.off()

