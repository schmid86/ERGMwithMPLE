
# S3 method
plot.MPLEergm <- function(object, vars.per.page=3, cex.axis=1, cex.lab=1, cex.main=1, lwd=2){

  
  stat.diff.mat <- t(object$obs.sim.statistics) - object$observed.statistics

  par(mfrow=c(vars.per.page, 2))
  for(i in 1:ncol(object$obs.sim.statistics)){

    plot(stat.diff.mat[i,], type="l", main= paste("Trace of ", colnames(object$obs.sim.statistics)[i]),
         ylab=" ", xlab="Iterations", cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab )
    abline(h=0, col="red", lwd=lwd)
    hist(stat.diff.mat[i,], main= paste("Density of ", colnames(object$obs.sim.statistics)[i]),
         ylab=" ", xlab=" ", cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab)
    abline(v=0, col="red", lwd=lwd)

  }

}


