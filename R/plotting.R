#' Plot colvarfile object
#'
#' `plot.colvarfile` plots colvarfile object. For a colvarfile with one collective variable it plots its evolution.
#' For a colvarfile with two collective variables it plots CV1 vs. CV2. For a colvarfile with more collective
#' variables it does not work, but you can still use $times and $cvs with default 'plot' function.
#'
#' @param x colvarfile object.
#' @param ignoretime time in the first column of the colvarfile will be ignored.
#' @param main an overall title for the plot: see 'title'.
#' @param sub a sub title for the plot: see 'title'.
#' @param xlab a title for the x axis: see 'title'.
#' @param ylab a title for the y axis: see 'title'.
#' @param asp the y/x aspect ratio, see 'plot.window'.
#' @param pch plotting 'character', i.e., symbol to use. See 'points'.
#' @param col color code or name, see 'par'.
#' @param bg background (fill) color for the open plot symbols given by
#'        'pch = 21:25'.
#' @param cex character (or symbol) expansion: a numerical vector. This
#'        works as a multiple of 'par("cex")'.
#' @param lwd line width for drawing symbols see 'par'.
#' @param xlim numeric vector of length 2, giving the x coordinates range.
#' @param ylim numeric vector of length 2, giving the y coordinates range.
#' @param axes a logical value indicating whether both axes should be drawn
#'        on the plot.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @examples
#' plot(acealanmeCVs)
plot.colvarfile<-function(x, ignoretime=FALSE,
                          xlab=NULL, ylab=NULL,
                          xlim=NULL, ylim=NULL,
                          main=NULL, sub=NULL,
                          pch=1, col="black", bg="red", cex=1,
                          asp=NULL, lwd=1, axes=TRUE,...) {
  cvfile <-x
  xlims<-NULL
  ylims<-NULL
  if(!is.null(xlim)) {xlims<-xlim}
  if(!is.null(ylim)) {ylims<-ylim}
  if(cvfile$ncvs==1) {
    if(is.null(xlab)) xlab="time"
    if(is.null(ylab)) ylab="CV"
    if(is.null(cvfile$times)) {
      if(cvfile$ncvs==1) {
        times <- 1:length(cvfile$cvs)
      } else {
        times <- 1:nrow(cvfile$cvs)
      }
    } else {
      times <- cvfile$times
    }
    if(ignoretime) {
      plot(seq(from=times[1],by=times[2],length.out=length(times)),
           cvfile$cvs, type="l",
           xlab=xlab, ylab=ylab,
           main=main, sub=sub,
           xlim=xlims, ylim=ylims,
           col=col, cex=cex, lwd=lwd,
           asp=asp, axes=axes)
    } else {
      plot(times, cvfile$cvs, type="l",
           xlab=xlab, ylab=ylab,
           main=main, sub=sub,
           xlim=xlims, ylim=ylims,
           col=col, cex=cex, lwd=lwd,
           asp=asp, axes=axes)
    }
  }
  if(cvfile$ncvs==2) {
    if(is.null(xlab)) xlab="CV1"
    if(is.null(ylab)) ylab="CV2"
    plot(cvfile$cvs[,1], cvfile$cvs[,2], type="p",
         xlab=xlab, ylab=ylab,
         main=main, sub=sub,
         xlim=xlims, ylim=ylims,
         pch=pch, col=col, bg=bg, cex=cex, lwd=lwd,
         asp=asp, axes=axes)
  }
  if(cvfile$ncvs>2) {
    cat("plot.colvarfile does not work for a colvarfile with more collective\n")
    cat("variables than 2, but you can still use $times and $cvs with default\n")
    cat("plot function.\n")
  }
}

#' Plot points for colvarfile object
#'
#' `points.colvarfile` plots points for colvarfile object. For a colvarfile with one
#' collective variable it plots its evolution. For a colvarfile with two collective
#' variables it plots CV1 vs. CV2. For a colvarfile with more collective
#' variables it does not work, but you can still use $times and $cvs with default
#' 'plot' and 'points' function.
#'
#' @param x colvarfile object.
#' @param ignoretime time in the first column of the colvarfile will be ignored.
#' @param pch plotting 'character', i.e., symbol to use. See 'points'.
#' @param col color code or name, see 'par'.
#' @param bg background (fill) color for the open plot symbols given by
#'        'pch = 21:25'.
#' @param cex character (or symbol) expansion: a numerical vector. This
#'        works as a multiple of 'par("cex")'.
#' @param lwd line width for drawing symbols see 'par'.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @examples
#' plot(acealanmeCVs)
#' points(acealanmeCVs, col="red")
points.colvarfile<-function(x, ignoretime=FALSE,
                           pch=1, col="black", bg="red", cex=1,
                           lwd=1, ...) {
  cvfile <-x
  if(cvfile$ncvs==1) {
    if(is.null(cvfile$times)) {
      if(cvfile$ncvs==1) {
        times <- 1:length(cvfile$cvs)
      } else {
        times <- 1:nrow(cvfile$cvs)
      }
    } else {
      times <- cvfile$times
    }
    if(ignoretime) {
      points(seq(from=times[1],by=times[2],length.out=length(times)),
             cvfile$cvs, col=col, cex=cex, lwd=lwd)
    } else {
      points(times, cvfile$cvs, col=col, cex=cex, lwd=lwd)
    }
  }
  if(cvfile$ncvs==2) {
    points(cvfile$cvs[,1], cvfile$cvs[,2],
           pch=pch, col=col, bg=bg, cex=cex, lwd=lwd)
  }
  if(cvfile$ncvs>2) {
    cat("points.colvarfile does not work for a colvarfile with more collective\n")
    cat("variables than 2, but you can still use $times and $cvs with default\n")
    cat("plot/points function.\n")
  }
}

#' Plot lines for colvarfile object
#'
#' `lines.colvarfile` plots lines for colvarfile object. For a colvarfile with one
#' collective variable it plots its evolution. For a colvarfile with two collective
#' variables it plots CV1 vs. CV2. For a colvarfile with more collective
#' variables it does not work, but you can still use $times and $cvs with default
#' 'plot' and 'lines' function.
#'
#' @param x colvarfile object.
#' @param ignoretime time in the first column of the colvarfile will be ignored.
#' @param col color code or name, see 'par'.
#' @param lwd line width for drawing symbols see 'par'.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @examples
#' plot(acealanmeCVs)
#' lines(acealanmeCVs, col="red")
lines.colvarfile<-function(x, ignoretime=FALSE,
                           lwd=1, col="black",...) {
  cvfile <-x
  if(cvfile$ncvs==1) {
    if(is.null(cvfile$times)) {
      if(cvfile$ncvs==1) {
        times <- 1:length(cvfile$cvs)
      } else {
        times <- 1:nrow(cvfile$cvs)
      }
    } else {
      times <- cvfile$times
    }
    if(ignoretime) {
      lines(seq(from=times[1],by=times[2],length.out=length(times)),
            cvfile$cvs, col=col, lwd=lwd)
    } else {
      lines(times, cvfile$cvs, col=col, lwd=lwd)
    }
  }
  if(cvfile$ncvs==2) {
    lines(cvfile$cvs[,1], cvfile$cvs[,2],
          col=col, lwd=lwd)
  }
  if(cvfile$ncvs>2) {
    cat("lines.colvarfile does not work for a colvarfile with more collective\n")
    cat("variables than 2, but you can still use $times and $cvs with default\n")
    cat("plot/lines function.\n")
  }
}

#' Plot evolution of bias potential in colvarfile object
#'
#' `plotbias` plots evolution of bias potential in colvarfile.
#'
#' @param colvar colvarfile object.
#' @param ignoretime time in the first column of the HILLS file will be ignored.
#' @param main an overall title for the plot: see 'title'.
#' @param sub a sub title for the plot: see 'title'.
#' @param xlab a title for the x axis: see 'title'.
#' @param ylab a title for the y axis: see 'title'.
#' @param asp the y/x aspect ratio, see 'plot.window'.
#' @param col color code or name, see 'par'.
#' @param lwd line width for drawing symbols see 'par'.
#' @param xlim numeric vector of length 2, giving the x coordinates range.
#' @param ylim numeric vector of length 2, giving the y coordinates range.
#' @param axes a logical value indicating whether both axes should be drawn
#'        on the plot.
#'
#' @export
#' @examples
#' plotbias(acealanmeCVs)
plotbias<-function(cvfile, ignoretime=FALSE,
                   xlab=NULL, ylab=NULL,
                   xlim=NULL, ylim=NULL,
                   main=NULL, sub=NULL,
                   col="black", asp=NULL, lwd=1, axes=TRUE) {
  if(class(cvfile)=="colvarfile") {
    if(cvfile$bias==NULL) {
      stop("Error: Input colvarfile does not contain any bias potential")
    }
    if(is.null(xlab)) xlab="time"
    if(is.null(ylab)) ylab="bias potential"
    if(is.null(cvfile$times)) {
      if(cvfile$ncvs==1) {
        times <- 1:length(cvfile$cvs)
      } else {
        times <- 1:nrow(cvfile$cvs)
      }
    } else {
      times <- cvfile$times
    }
    if(ignoretime) {
      plot(seq(from=times[1],by=times[2],length.out=length(times)),
           cvfile$bias, type="l",
           xlab=xlab, ylab=ylab,
           main=main, sub=sub,
           col=col, lwd=lwd,
           asp=asp, axes=axes)
    } else {
      plot(times, cvfile$bias,
           xlab=xlab, ylab=ylab,
           main=main, sub=sub,
           col=col, lwd=lwd,
           asp=asp, axes=axes)
    }
  } else {
    stop("Error: Function plotbias requires object colvarfile as an input")
  }
}

