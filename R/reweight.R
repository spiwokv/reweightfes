library(metadynminer)

#weightvanthoff(cvfile, npoints=60, maxfe=100, temp=300, eunits="kJ/mol",
#               imin=1, imax=NULL, xlim=NULL, ylim=NULL)
#reweightbonomi(cvfile, npoints=60, maxfe=100, temp=300, eunits="kJ/mol",
#               imin=1, imax=NULL, xlim=NULL, ylim=NULL)


#' Calculate free energy surface from metadynamics by reweighting algorithm
#' by Tiwary and Parrinello, J. Phys. Chem. B (2015) <doi:10.1021/jp504920s>
#'
#' `reweighttiwary` calculates free energy surface from hills and colvars by
#' Tiwary&Parrinello algorithm.
#'
#' @param cvfile colvarfile object.
#' @param hillsfile hillsfile object.
#' @param nfes number of bias potential sums used to estimate c(t) (default 100).
#' @param temp temperature in Kelvins
#' @param eunit energy units (kJ/mol or kcal/mol, kJ/mol is default)
#' @param gamma bias factor of well-tempered metadynamics (default NULL).
#' @param imax index of a collective variable record from which summation stops (default the rest of hills).
#' @param xlim numeric vector of length 2, giving the CV1 coordinates range.
#' @param ylim numeric vector of length 2, giving the CV2 coordinates range.
#' @param npoints resolution of the free energy surface in number of points (default 60).
#' @param maxfe free energy of point on the output free energy surface with no sampling (default 100).
#' @param usefes2 logical parameter, if TRUE 'fes2' function is used instead of 'fes' (default FLASE).
#' @return fes object.
#'
#' @export
#' @examples
#' tfes<-reweighttiwary(cvs=acealanmeCVs, hills=acealanme, imax=5000)
reweighttiwary(cvs, hills, npoints=60, maxfe=100, nfes=100, temp=300,
               gamma=10, eunits="kJ/mol", imax=NULL, xlim=NULL, ylim=NULL, usefes2=F) {
  if(class(cvs)!="colvarfile") {
    stop("Error: Wrong colvarfile formate")
  }
  if(class(hills)!="hillsfile") {
    stop("Error: Wrong hillsfile formate")
  }
  if(nrow(cvs$cvs)==((hills$size[2]-3)/2)) {
    stop("Error: Different number of collective variables in colvarfile and hillsfile")
  }
  if(eunits=="kJ/mol") {
    beta <- 1000/8.314/temp
  } else {
    if(eunits=="kcal/mol") {
      beta <- 1000/8.314/temp/4.184
    } else {
      stop("Error: invalid energy unit, use kJ/mol or kcal/mol")
    }
  }
  if(nrow(cvs$cvs)==2) {
    minCV1 <- min(cvs$cvs[,1])
    maxCV1 <- max(cvs$cvs[,1])
    minCV2 <- min(cvs$cvs[,2])
    maxCV2 <- max(cvs$cvs[,2])
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    ylims<-c(minCV2-0.05*(maxCV2-minCV2), maxCV2+0.05*(maxCV2-minCV2))
    if(!is.null(xlim)) {xlims<-xlim}
    if((hills$per[1]==T)&is.null(xlim)) {xlims<-hills$pcv1}
    if(!is.null(ylim)) {ylims<-ylim}
    if((hills$per[2]==T)&is.null(ylim)) {ylims<-hills$pcv2}
  }
  if(nrow(cvs$cvs)==1) {
    minCV1 <- min(cvs$cvs[,1])
    maxCV1 <- max(cvs$cvs[,1])
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    if(!is.null(xlim)) {xlims<-xlim}
    if((hills$per[1]==T)&is.null(xlim)) {xlims<-hills$pcv1}
  }
  ebtac <- c()
  if(imax==NULL) {
    ihills <- hills$size[1]
    framestosum <- hills$size[1]/nfes
  } else {
    ihills <- imax*hills$size[1]/length(cvs$times)
    framestosum <- 1:ihills/nfes
  }
  if(usefes2) {
    tfes <- fes2(hills, imax=1, npoints=npoints, xlim=xlims, ylim=ylims) -
            fes2(hills, imax=1, npoints=npoints, xlim=xlims, ylim=ylims)
  } else {
    suppressWarnings(tfes <- fes(hills, imax=1, npoints=npoints, xlim=xlims, ylim=ylims) -
                             fes(hills, imax=1, npoints=npoints, xlim=xlims, ylim=ylims))
  }
  for(i in 1:nfes) {
    if(usefes2) {
      tfes <- tfes + fes2(hills, imin=(i-1)*framestosum+1, imax=i*framestosum,
                          npoints=npoints, xlim=xlims, ylim=ylims)
    } else {
      tfes <- tfes + fes(hills, imin=(i-1)*framestosum+1, imax=i*framestosum,
                         npoints=npoints, xlim=xlims, ylim=ylims)
    s1 <- sum(exp(-beta*tfes$fes))
    s2 <- sum(exp(-beta*tfes$fes/gamma))
    ebtac<-c(ebtac,s1/s2)
  }
  if(cvs$bias!=NULL) {
    if(is.null(imax)) {
      bp <- cvs$bias
    } else {
      bp <- cvs$bias[1:imax]
    }
  } else {
    stop("Error: Input colvarfile does not contain bias potential")
  }
  ebtacc <- rep(ebtac, each=nrow(cvs$cvs)/nfes)
  if(length(ebtacc)>nrow(cvs$cvs)) ebtacc<-ebtacc[1:nrow(cvs$cvs)]
  if(length(ebtacc)<nrow(cvs$cvs)) ebtacc[length(ebtacc):nrow(cvs$cvs)]<-ebtacc[length(ebtacc)]
  if(nrow(cvs$cvs)==2) {
    x<-0:(npoints-1)*(xlims[2]-xlims[1])/(npoints-1)+xlims[1]
    y<-0:(npoints-1)*(ylims[2]-ylims[1])/(npoints-1)+ylims[1]
    if(is.null(imax)) {
      icv1 <- ceiling((cvs$cvs[,1]-xlims[1])*npoints/(xlims[2]-xlims[1]))
      icv2 <- ceiling((cvs$cvs[,2]-ylims[1])*npoints/(ylims[2]-ylims[1]))
    } else {
      icv1 <- ceiling((cvs$cvs[1:imax,1]-xlims[1])*npoints/(xlims[2]-xlims[1]))
      icv2 <- ceiling((cvs$cvs[1:imax,2]-ylims[1])*npoints/(ylims[2]-ylims[1]))
    }
    ofes <- matrix(rep(0, npoints*npoints), nrow=npoints)
    for(i in 1:npoints) {
      for(j in 1:npoints) {
        ofes[i,j]<-sum((icv1==i)*(icv2==j)*exp(beta*bp)/ebtacc)
      }
    }
    ofes <- -log(ofes)/beta
    ofes <- ofes - min(ofes)
    ofes[ofes==Inf]<-maxfe
    cfes<-list(fes=ofes, hills=hills$hillsfile, rows=npoints, dimension=2, per=hills$per, x=x, y=y, pcv1=hills$pcv1, pcv2=hills$pcv2)
    class(cfes) <- "fes"
  }
  if(nrow(cvs$cvs)==1) {
    minCV1 <- min(cvs$cvs[,1])
    maxCV1 <- max(cvs$cvs[,1])
    xlims<-c(minCV1-0.05*(maxCV1-minCV1), maxCV1+0.05*(maxCV1-minCV1))
    if(!is.null(xlim)) {xlims<-xlim}
    if((hills$per[1]==T)&is.null(xlim)) {xlims<-hills$pcv1}
    x<-0:(npoints-1)*(xlims[2]-xlims[1])/(npoints-1)+xlims[1]
    if(is.null(imax)) {
      icv1 <- ceiling((cvs$cvs[,1]-xlims[1])*npoints/(xlims[2]-xlims[1]))
    } else {
      icv1 <- ceiling((cvs$cvs[1:imax,1]-xlims[1])*npoints/(xlims[2]-xlims[1]))
    }
    ofes <- rep(0, npoints)
    for(i in 1:npoints) {
      ofes[i]<-sum((icv1==i)*exp(beta*bp)/ebtacc)
    }
    ofes <- -log(ofes)/beta
    ofes <- ofes - min(ofes)
    ofes[ofes==Inf]<-maxfe
    cfes<-list(fes=ofes, hills=hills$hillsfile, rows=npoints, dimension=1, per=hills$per, x=x, pcv1=hills$pcv1)
    class(cfes) <- "fes"
  }
  return(cfes)
}

