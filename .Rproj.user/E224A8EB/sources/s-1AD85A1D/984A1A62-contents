

# dependencies: ergm, doParallel, foreach
#library(ergm)
#library(doParallel)
#library(foreach)


MPLEergm <- function(formula, R= 500, control.mple = control.ergm(drop=TRUE),
                     method="bootstrap",
                     control.sim = control.simulate.ergm(MCMC.interval=1024),
                     checkforseparateddata=TRUE,
                     control.boot=control.ergm(drop=TRUE), nowonit=100, cpu.cores=1, constraints =~.){


  if(checkforseparateddata==TRUE){

    X <- ergmMPLE(formula)
    data <- as.data.frame(cbind(X$response, X$predictor))
    colnames(data)[1]<- "y"
    w <- X$weights

    # update formula, so the intercept in glm is not calculated twice
    # form.glm <- update(formula, .~.-1)

    model <- glm(y ~.-1, data=data, family=binomial, weights=w)
  } # end if checkforseparateddata


  # fit model using MPLE
  mple <- ergm(formula, estimate="MPLE" , control=control.mple, constraints = constraints)

  # get the MPLE
  theta.mple<- coef(mple)

  observed.statistics <- summary(formula)

  # get number of variables in model
  num.variables <- length(coef(mple))

  #simulate network from mple
  cat("Simulating ", R , " networks\n")
  sim.mple <- simulate(mple, nsim=R, control=control.sim ) # simulate network



  if(method=="bootstrap"){
  # create empty matrix to store mple of bootstrap samples
  boot.mple.mat <- matrix(0, R, num.variables)
  colnames(boot.mple.mat) <- names(coef(mple))

  # create empty matrix to store network statistics of bootstrap samples
  boot.stat.mat <- matrix(0, R, num.variables)
  colnames(boot.stat.mat) <- names(coef(mple))

  # replace response network in formula
  tt <- terms(formula)
  new.formula <- reformulate(attr(tt, "term.labels"), "sim.mple[[i]]")

  cat("Estimating MPLE of simulated networks\n")

  # if only one CPU core is used

  if(cpu.cores==1){
    for(i in 1:R){

      if(i %% nowonit == 0)
      {
        cat("Now on iteration ", i, "\n")
      } # end if

      boot.mple<- suppressMessages(ergm(new.formula, estimate="MPLE" , control=control.boot, constraints=constraints))

      boot.mple.mat[i,] <- coef(boot.mple)

      boot.stat.mat[i,] <- summary(new.formula)
    } # end for


  } else{ # parallel

    registerDoParallel(cores=cpu.cores)

    par.results <-foreach(i=1:R, .export='ergm', .combine=rbind)%dopar%{



      tt <- terms(formula)
      sim.mple.net <- sim.mple[[i]]
      par.formula <- reformulate(attr(tt, "term.labels"), "sim.mple.net")

      boot.mple<- suppressMessages(ergm(par.formula, estimate="MPLE" , control=control.boot))

      g1 <- coef(boot.mple)

      g2 <- summary(par.formula)

      g<- c(g1, g2)

      g

    } # end foreach

    # split par.results into boot.mple.mat and boot.stat.mat

    boot.mple.mat <- par.results[,1:num.variables ]
    boot.stat.mat <- par.results[,(num.variables+1):(2*num.variables) ]

  } # end else # parallel

  # create empty matrix
  Result.matrix <- matrix(0, num.variables, 3)
  rownames(Result.matrix) <- names(coef(mple))
  colnames(Result.matrix) <- c("MPLE", "2.5%", "97.5%")

  Result.matrix[,1] <- coef(mple)

  for(i in 1:num.variables){
    Result.matrix[i,2] <- quantile(boot.mple.mat[,i], probs = 0.025)
    Result.matrix[i,3] <- quantile(boot.mple.mat[,i], probs = 0.975)

  }

  return.list <- list(estimate="bootstrap", boot.interval = Result.matrix , boot.mple.matrix=boot.mple.mat , MPLE= mple, formula=formula,
                      observed.statistics=observed.statistics, obs.sim.statistics=boot.stat.mat)
  class(return.list) <- "MPLEergm"

  return(return.list)

  }else{  # end if method=bootstrap
    if(method=="godambe"){

      # get change statistics
      net.sim.dat <- ergmMPLE(formula)

      ##################################################
      # caluclation of J(theta) = -\delta u(theta,y)

      cat("Calculating Sensitivity Matrix\n")

      resp<- net.sim.dat[[1]]
      X<- net.sim.dat[[2]]
      weight <- net.sim.dat[[3]]


      J <- matrix(0,num.variables,num.variables)


      for(h in 1:length(weight)){
        for(k in 1:num.variables){
          for(l in 1:num.variables){

            J[k,l] <- J[k,l] + weight[h]*(exp(sum(theta.mple*X[h,]))/(1+exp(sum(theta.mple*X[h,])))^2 * X[h,k]* X[h,l])

          } # for l
        }  # for k

      } # for h

      J.inv <- solve(J)

      ##########################################################################
      # calculation of V(theta) = Var(u(theta,y)) using the sim.num networks

      cat("Estimating Variability Matrix\n")
      net.stat <- matrix(0, nrow=length(sim.mple), ncol=num.variables)
      colnames(net.stat) <- names(coef(mple))
      u.data <- matrix(0,nrow=length(sim.mple), ncol=num.variables)


      for(i in 1:length(sim.mple)){

        # replace response network in formula
        tt <- terms(formula)
        new.formula <- reformulate(attr(tt, "term.labels"), "sim.mple[[i]]")

        dat<- ergmMPLE(new.formula )

        net.stat[i,]<- summary(new.formula)


        resp<- dat[[1]]
        X<- dat[[2]]
        weight <- dat[[3]]


        u.single <- rep(0, num.variables)
        for(h in 1:length(weight)){
          for(k in 1:num.variables){

            u.single[k]<- u.single[k] + weight[h]*(X[h,k]*(resp[h] - exp(sum(theta.mple*X[h,]))/(1+exp(sum(theta.mple*X[h,]))) ))

          } # for k
        } # for h

        u.data[i,] <- u.single


      } # end for i


      # calculate V.hat by estimating sd

      u.mean <- colMeans(u.data)
      #u.mean


      u.sum <- matrix(0,num.variables,num.variables)
      for(i in 1: nrow(u.data)){

        u.diff <- u.data[i,]- u.mean
        u.sum <- u.sum + u.diff%*%t(u.diff)

      }

      u.sum.n <- u.sum/(nrow(u.data)-1)
      #u.sum.n


      G <- J.inv%*%u.sum.n%*%J.inv

      sd.godambe <- rep(0, num.variables)
      for(i in 1:num.variables){
        sd.godambe[i] <- sqrt(G[i,i])
      } # end for i


      # create empty matrix
      Result.matrix <- matrix(0, num.variables, 5)
      rownames(Result.matrix) <- names(coef(mple))
      colnames(Result.matrix) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "")

      Result.matrix[,1] <- coef(mple)
      Result.matrix[,2] <- sd.godambe
      z <- coef(mple)/sd.godambe
      Result.matrix[,3] <- z
      p <- 2*pnorm(-abs(z))
      Result.matrix[,4] <- p
      Result.matrix <- as.data.frame(Result.matrix)

      # stars
      stars<- rep(0, num.variables)
      for(i in 1:length(p)){
        if(p[i]>0.1){
          stars[i]<- ''
        }else{
          if(p[i]>0.05){
            stars[i] <- '.'
          }else{
            if(p[i]>0.01){
              stars[i] <- '*'
            }else{
              if(p[i]>0.001){
                stars[i] <- '**'
              }else{
                stars[i] <- '***'
              }
            }
          }
        }

      }

      Result.matrix[,5] <- stars
      names(Result.matrix)[5]<- ""


      return.list <- list(estimate="godambe" , coefficient = coef(mple) , MPLE= mple, formula=formula,
                          observed.statistics=observed.statistics, G=G, invHess = J.inv,
                          sd.godambe=sd.godambe, obs.sim.statistics = net.stat, summary.results = Result.matrix)
      class(return.list) <- "MPLEergm"

      return(return.list)



    }# end if method==godambe
  } # end else

}


#m2 <- ergm.mple(faux.mesa.high ~ edges + nodematch("Sex")+gwesp(0.25,fixed=TRUE), R=100,
#               control.mple=control.ergm(MCMC.samplesize=2000), method="godambe",
#               control.sim = control.simulate.ergm(MCMC.burnin=5000, MCMC.interval=10000),
#               cpu.cores=1)


# starting matrix
#set.seed(123123)
#A <- matrix(rbinom(81, 1,0.1), 9,9)
#diag(A)<- 0
#A <- as.network.matrix(A, directed=FALSE)

#summary(A~edges+triangles)
#plot(A)

#m3 <- ergm.mple(A~edges+triangles, R=100)
