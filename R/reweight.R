library(metadynminer)
hillsf<-read.hills("HILLS", per=c(T,T))
colvar<-read.table("COLVAR")
nfes <- 100
temp <- 2.49*1000/8.314
gamma <- 10
ebtac <- c()
framestosum <- hillsf$size[1]/nfes
tfes <- fes(hillsf, imax=1)-fes(hillsf, imax=1)
for(i in 1:nfes) {
  tfes <- tfes + fes(hillsf, imin=(i-1)*framestosum+1, imax=i*framestosum)
  s1 <- sum(exp(-1000*tfes$fes/8.314/temp))
  s2 <- sum(exp(-1000*tfes$fes/8.314/temp/gamma))
  ebtac<-c(ebtac,s1/s2)
}
cvsandbp <- data.frame(colvar[,2:4])

#feshist<-function(x) {
#  omat <- matrix(rep(0, 60*60), nrow=60)
#  icv1 <- (x[1,1]+pi)*30/pi
#  icv2 <- (x[1,2]+pi)*30/pi
#  omat[icv1,icv2] <- exp(x[1,3]/8.314/0.3)
#  return(omat)
#}
#fes <- matrix(rep(0, 60*60), nrow=60)
#mlist <- list()
#for(i in 1:nrow(cvsandbp)) {
#  mlist <- list()
#  mlist[[1]] <- fes
#  mlist[[2]] <- feshist(cvsandbp[i,])
#  fes <- Reduce('+', mlist)
#}

for(i in 1:60) {
  for(j in 1:60) {
    

fes <- -8.314*0.3*log(fes)
fes <- fes - min(fes)
image(fes)



