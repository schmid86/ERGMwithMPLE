# goodness of fit for bootERGM object

gof.MPLEergm<- function(object, coef=NULL, GOF=NULL,
                        constraints=NULL, control=control.gof.ergm(), verbose=FALSE){

  g1 <- gof(object$MPLE, coef=coef, GOF=GOF, constraints=constraints,
            control=control, verbose=verbose)
}
