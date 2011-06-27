my.IC <- function(penden.env) {
  values <- -eigen(get("Derv2.cal",penden.env))$values
  rho <- min(values[values>1e-09])
  d <- dim(get("Derv2.cal",penden.env))[1]
  #mytrace2 <- sum(diag(solve(get("Derv2.pen",penden.env)+rho*diag(d))%*%(get("Derv2.cal",penden.env)+rho*diag(d))))
  #mytrace <- sum(diag(solve(get("Derv2.pen",penden.env))%*%(get("Derv2.cal",penden.env)))) #AIC2=(2*get("log.like",penden.env)+2*mytrace2)
  mytrace <- sum(diag(my.positive.definite.solve(get("Derv2.pen",penden.env))%*%(get("Derv2.cal",penden.env))))
  #return(list(myAIC=(-2*get("log.like",penden.env)+2*mytrace),mytrace=mytrace,myBIC=(2*get("log.like",penden.env)+mytrace*log(get("n",penden.env)))))
  assign("AIC",-2*get("log.like",penden.env)+2*mytrace,penden.env)
  assign("BIC",2*get("log.like",penden.env)+mytrace*log(get("n",penden.env)),penden.env)
  assign("mytrace",mytrace,penden.env)
}
