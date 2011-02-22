pencopula <- function(data,d=3,D=d,q=1,base="B-spline",max.iter=20,plot.bsp=FALSE,
                    lambda=NULL,pen.order=2,adapt.grid=FALSE,add=TRUE,
                    alpha=0,symmetric=TRUE,data.frame=parent.frame()) {

  library(quadprog)
  library(fda)
  library(spam)
  library(lattice)

  penden.env <- new.env()
  assign("frame",data.frame,penden.env)
  assign("Y",data,penden.env)
  
  if(is.matrix(data)|is.data.frame(data)) den.form <-  pendenForm(penden.env)
  else stop(" 'data' is not a matrix or a data frame. Create 'data' column by column.")

  if(is.null(d)) d <- 3 # Anzahl Halbierungen des Intervals [0,1]
  if(is.null(D)) D <- d # max. Hierachiestufe

  if(D<d) d <- D
  if(base=="Bernstein") {
    D <- 2*d
    pen.order <- 1
  }
  assign("D",D,penden.env)
  assign("d",d,penden.env)
  
  assign("add",add,penden.env)
  assign("pen.order",pen.order,penden.env)
  assign("alpha",alpha,penden.env)
  assign("symmetric",symmetric,penden.env)
  assign("no",FALSE,penden.env)
  assign("max.iter",max.iter,penden.env)
  assign("base",base,penden.env)
  assign("adapt.grid",adapt.grid,penden.env)
  assign("plot.bsp",plot.bsp,penden.env)
  
  dd <- (2**d)+1 # Anzahl Knoten
  assign("dd",dd,penden.env)

  if(is.null(q)) q <- 1 # Grad des B spline
  assign("q",q,penden.env)

  ddb <- dd + (q-1) # Anzahl Basisfunktionen
  assign("ddb",ddb,penden.env)
 
  p <- get("p",penden.env)  #p Dimension der Kovariablen
  
  assign("D",D,penden.env) #max. Hierarchiestufe

  dimension <- c(rep(0,q+1),rep(1:d,2**(0:(d-1))))

  if(is.null(lambda)) lambda <- rep(10000,p)
  if(!is.null(lambda)) if(length(lambda)<p | length(lambda)>p) stop("length of lambda is wrong")
  assign("lambda",lambda,penden.env)

  assign("dimension",dimension,penden.env)
# Dimension gibt die Hierarchiestufe an, aus der der hierarchische B-spline berechnent wird
  
  ##################

  #D maximale Hierarchiestufe
  DIMENSION <- dimension
  Index.basis <- matrix(1:ddb)
  index.sparse <- DIMENSION <= D
  Index.basis.D <- matrix(Index.basis[index.sparse,])
  DIMENSION <- DIMENSION[index.sparse]
 
  for ( j in 2:p)
    {
      #print(j)
      DIMENSION.j <-  kronecker(matrix(1,ddb,1),DIMENSION) + kronecker( dimension, matrix(1, length(DIMENSION),1))
      Index.basis.plus.1 <- matrix(NA, dim(Index.basis.D)[1] * ddb , j)
      Index.basis.plus.1[,j] <- kronecker(matrix(1:ddb), matrix(1,dim(Index.basis.D)[1],1))
      Index.basis.plus.1[, 1:(j-1)] <-  kronecker(matrix(1, ddb,1),Index.basis.D)
      index.sparse <- DIMENSION.j <= D
      Index.basis.D <- Index.basis.plus.1[index.sparse,]
      DIMENSION <- DIMENSION.j[index.sparse]
    }
  DD <- dim(Index.basis.D)[1] # Dimension of sparse grid basis

  assign("DD",DD,penden.env) # DD Anzahl Koeffizienten
  assign("Index.basis.D",Index.basis.D,penden.env)

  ###################

  # Matrix zur Erstellung der marginalen Spline Koeffizienten
  # 
  j <- 1
  T.marg <- array(NA, dim=c(ddb,DD,p))
  
  for ( j in 1:p)
    {
      for ( l in 1:ddb)
        {
          T.marg[l,,j] <- (Index.basis.D[,j] == l)+0
        }
    }
  
  assign("T.marg",T.marg,penden.env)

  ####################

  tilde.Psi.d <-  array(NA, dim=c(get("n",penden.env),ddb,p))
 
  for (j in 1:p)
    {
      tilde.Psi.d[,,j] <-  hierarch.bs(get("Y",penden.env)[,j], d = d, plot.bsp = plot.bsp,typ=3,penden.env,int=FALSE)$B.tilde
    }

  assign("tilde.Psi.d",tilde.Psi.d,penden.env)

  assign("tilde.PSI.d.D",tilde.Psi.d[,Index.basis.D[,1],1],penden.env)

  for (j in 2:p)
    {
      assign("tilde.PSI.d.D",get("tilde.PSI.d.D",penden.env) * get("tilde.Psi.d",penden.env)[,Index.basis.D[,j],j],penden.env)
    }


  ####################

  #knots
  knots.start(penden.env)
 
  #####################

  #startwerte und gitter berechnen

  start.valgrid(penden.env)

  #####################################

  # Allgemein: Die Restriktionsmatrizen sind als Array aufgefasst:
  
  A <- array(NA, dim=c(get("ddb",penden.env),DD,p))

  for ( j in 1:p)
    {
      A[,,j] <- get("tilde.Psi.knots.d",penden.env) %*% T.marg[,,j]
    }

  assign("A.Restrict",A,penden.env)
  
  # Nun muss gelten A[,,j] %*% c = 1 fuer alle j

  #############################

  penalty.matrix(penden.env=penden.env)

  #############################

  liste <- matrix(0,1,3+DD+p)
  n.liste <- matrix(0,1,3+DD+p)

  lam <- coef <- c()
  for(i in 1:p) lam[i] <- paste("lambda.",i,sep="")
  for(j in 1:DD) coef[j] <- paste("b.",j,sep="")
 
  colnames(liste) <- c("pen.log.like","log.like","marg.log.like",lam,coef)
  
  #print("ck at the beginning")
   
  #print("fitted values at the beginning")
  help.str <- paste("d=",get("d",penden.env),"D=",get("D",penden.env),"lambda=",get("lambda",penden.env)[1],sep="")
  assign("help.str",help.str,penden.env)

  f.hat.val(penden.env,cal=TRUE)
    if(get("no",penden.env)) {
    assign("pen.log.like",0,penden.env)
    assign("log.like",0,penden.env)
    assign("AIC",0,penden.env)
    assign("BIC",0,penden.env)
    obj <- list(penden.env=penden.env)
    class(obj) <- "pencopula"
    return(obj)
  }

  #print(get("f.hat.val",penden.env))

  #print("marg log-likelihood at the beginning")
  pen.log.like(penden.env,cal=TRUE)
  #print(get("log.like",penden.env))
  Derv1(penden.env)
  Derv2(penden.env)
  #if(!fix.lambda) marg.likelihood(penden.env)
  #print(get("marg.log.like",penden.env))

  assign("i",i <- 1,penden.env)
  liste[i,1] <- get("pen.log.like",penden.env)
  liste[i,2] <- get("log.like",penden.env)
  #if(!fix.lambda) liste[i,3] <- get("marg.log.like",penden.env)
  liste[i,(4:(4+p-1))] <- get("lambda",penden.env)
  liste[i,((4+p):(4+p+DD-1))] <- get("ck.val",penden.env)

  assign("liste",liste,penden.env)
  
  assign("calc",TRUE,penden.env)
  
  Derv1(penden.env)
  Derv2(penden.env)
  if(new.weights(penden.env)=="fehler"){
    assign("pen.log.like",0,penden.env)
    assign("log.like",0,penden.env)
    assign("AIC",0,penden.env)
    assign("BIC",0,penden.env)
    obj <- list(penden.env=penden.env)
    class(obj) <- "pencopula"
    return(obj)
  }

  #if(!fix.lambda) {
  #  new.lambda(penden.env)
  #  penalty.matrix(penden.env,temp=TRUE)
  #  pen.log.like(penden.env,temp.lambda=TRUE)
  #  marg.likelihood(penden.env,temp=TRUE)
  #}
  #else assign("no.lambda",FALSE,penden.env)

  my.loop(penden.env)
  
  Derv1(penden.env)
  Derv2(penden.env)

  my.IC(penden.env)

  #obj <- list(penden.env=penden.env)
  class(penden.env) <- "pencopula"
  return(penden.env)
}
