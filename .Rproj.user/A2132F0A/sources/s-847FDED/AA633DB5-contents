# Aux 1 : guess_rank for r ------------------------------------------------
#' @keywords internal
guess_rank <- function(X,nnz){
  maxiter = 10000
  n = nrow(X)
  m = ncol(X)
  epsilon = nnz/sqrt(m*n)
  svdX = svd(X)
  S0 = svdX$d
  
  nsval0 = length(S0)
  S1 = S0[1:(nsval0-1)]-S0[2:nsval0]
  nsval1 = length(S1)
  if (nsval1>10){
    S1_ = S1/mean(S1[(nsval1-10):nsval1])
  } else {
    S1_ = S1/mean(S1[1:nsval1])
  }
  r1 = 0
  lam = 0.05
  
  itcounter = 0
  while (r1<=0){
    itcounter = itcounter+1
    cost = array(0,c(1,length(S1_)))
    for (idx in 1:length(S1_)){
      cost[idx] = lam*max(S1_[idx:length(S1_)]) + idx
    }
    v2 = min(cost)
    i2 = which(cost==v2)
    if (length(i2)==1){
      r1 = i2-1
    } else {
      r1 = max(i2)-1
    }
    lam = lam+0.05
    if (itcounter > maxiter){
      break
    }
  }
  
  if (itcounter<=maxiter){
    cost2 = array(0,c(1,(length(S0)-1)))
    for (idx in 1:(length(S0)-1)){
      cost2[idx] = (S0[idx+1]+sqrt(idx*epsilon)*S0[1]/epsilon)/S0[idx]
    }
    v2 = min(cost2)
    i2 = which(cost2==v2)
    if (length(i2)==1){
      r2 = i2
    } else {
      r2 = max(i2)
    }
    
    if (r1>r2){
      r = r1
    } else {
      r = r2
    }
    return(r)
  } else {
    r = min(nrow(X),ncol(X))
  }
}


# Aux 2 : compute the distortion ------------------------------------------
#' @keywords internal
aux_G <- function(X,m0,r){
  z = rowSums(X^2)/(2*m0*r)
  y = exp((z-1)^2) - 1
  idxfind = (z<1)
  y[idxfind] = 0
  out = sum(y)
  return(out)
}
#' @keywords internal
aux_F_t <- function(X,Y,S,M_E,E,m0,rho){
  n = nrow(X)
  r = ncol(X)
  
  out1 = (sum((((X%*%S%*%t(Y))-M_E)*E)^2))/2
  out2 = rho*aux_G(Y,m0,r)
  out3 = rho*aux_G(X,m0,r)
  out  = out1+out2+out3
  return(out)
}


# Aux 3 : compute the gradient --------------------------------------------
#' @keywords internal
aux_Gp <- function(X,m0,r){
  z = rowSums(X^2)/(2*m0*r)
  z = 2*exp((z-1)^2)/(z-1)
  idxfind = (z<0)
  z[idxfind] = 0
  
  out = (X*matrix(z,nrow=nrow(X),ncol=ncol(X),byrow=FALSE))/(m0*r)
}
#' @keywords internal
aux_gradF_t <- function(X,Y,S,M_E,E,m0,rho){
  n = nrow(X)
  r = ncol(X)
  m = nrow(Y)
  if (ncol(Y)!=r){
    stop("dimension error from aux_gradF_t")
  }
  
  XS  = (X%*%S)
  YS  = (Y%*%t(S))
  XSY = (XS%*%t(Y))
  
  Qx = ((t(X) %*% ((M_E-XSY)*E) %*% YS)/n)
  Qy = ((t(Y) %*% t((M_E-XSY)*E) %*% XS)/m)
  
  W = (((XSY-M_E)*E) %*% YS) + (X%*%Qx) + rho*aux_Gp(X,m0,r)
  Z = (t((XSY-M_E)*E) %*% XS)+ (Y%*%Qy) + rho*aux_Gp(Y,m0,r)
  
  resgrad = list()
  resgrad$W = W
  resgrad$Z = Z
  return(resgrad)
}


# Aux 4 : Sopt given X and Y ----------------------------------------------
#' @keywords internal
aux_getoptS <- function(X,Y,M_E,E){
  n = nrow(X)
  r = ncol(X)
  
  C = (t(X) %*% (M_E) %*% Y)
  C = matrix(as.vector(C))
  
  nnrow = ncol(X)*ncol(Y)
  A = matrix(NA,nrow=nnrow,ncol=(r^2))
  
  for (i in 1:r){
    for (j in 1:r){
      ind = ((j-1)*r+i)
      tmp = t(X) %*% (outer(X[,i],Y[,j])*E) %*% Y
      
      A[,ind] = as.vector(tmp)
    }
  }
  
  S = solve(A,C)
  out = matrix(S,nrow=r)
  return(out)
}

# Aux 5 : optimal line search ---------------------------------------------
#' @keywords internal
aux_getoptT <- function(X,W,Y,Z,S,M_E,E,m0,rho){
  norm2WZ = (norm(W,'f')^2)+(norm(Z,'f')^2)
  f = array(0,c(1,21))
  f[1] = aux_F_t(X,Y,S,M_E,E,m0,rho)
  t = -1e-1
  for (i in 1:20){
    f[i+1] = aux_F_t(X+t*W,Y+t*Z,S,M_E,E,m0,rho)
    if ((f[i+1]-f[i]) <= 0.5*t*norm2WZ){
      out = t
      break
    }
    t = t/2
  }
  out = t
  return(t)
}
