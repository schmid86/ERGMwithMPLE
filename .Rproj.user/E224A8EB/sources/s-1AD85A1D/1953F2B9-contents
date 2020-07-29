summary.MPLEergm<- function(object){
  
  
  if(object$estimate =="bootstrap"){
  
  cat("================================================\n")
  cat("          Parametric Bootstrap Results          \n")
  cat("================================================\n")
  cat(" \n")
  cat("Formula: ")
  print(object$formula)
  cat(" \n")
  cat("Bootstrap Samples: ", nrow(object$boot.mple.matrix), "\n")
  cat(" \n")
  cat("95% Bootstrap Confidence Intervals: \n")
  cat(" \n")
  print(object$boot.interval)
  
  }else{
    if(object$estimate == "godambe"){
      
      cat("================================================\n")
      cat("              Summary of Model Fit              \n")
      cat("================================================\n")
      cat(" \n")
      cat("Formula: ")
      print(object$formula)
      cat(" \n")
      cat("Simulated Networks: ", nrow(object$obs.sim.statistics), "\n")
      cat(" \n")
      cat("Maximum Pseudolikelihood Results: \n")
      cat(" \n")
      print(object$summary.results)
      cat("---\n")
      cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 '' 1")
    }
  }
  
}


