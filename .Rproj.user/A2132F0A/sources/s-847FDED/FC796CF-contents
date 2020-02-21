#' OptSpace : an algorithm for matrix reconstruction from a partially revealed set
#'
#' Let's assume an ideal matrix \eqn{M} with \eqn{(m\times n)} entries with rank \eqn{r} and
#' we are given a partially observed matrix \eqn{M\_E} which contains many missing entries.
#' Matrix reconstruction - or completion - is the task of filling in such entries.
#' OptSpace is an efficient algorithm that reconstructs \eqn{M} from \eqn{|E|=O(rn)} observed elements
#' with relative root mean square error (RMSE)
#' \deqn{RMSE \le C(\alpha)\sqrt{nr/|E|}}
#'
#' @param A an \eqn{(n\times m)} matrix whose missing entries should be flaged as NA.
#' @param ropt \code{NA} to guess the rank, or a positive integer as a pre-defined rank.
#' @param niter maximum number of iterations allowed.
#' @param tol stopping criterion for reconstruction in Frobenius norm.
#' @param showprogress a logical value; \code{TRUE} to show progress, \code{FALSE} otherwise.
#'
#' @return a named list containing
#' \describe{
#' \item{X}{an \eqn{(n \times r)} matrix as left singular vectors.}
#' \item{S}{an \eqn{(r \times r)} matrix as singular values.}
#' \item{Y}{an \eqn{(m \times r)} matrix as right singular vectors.}
#' \item{dist}{a vector containing reconstruction errors at each successive iteration.}
#' }
#'
#'
#' @references
#' \insertRef{keshavan_matrix_2010}{ROptSpace}
#'
#' @examples
#' ## Parameter Settings
#' n = 1000;
#' m = 100;
#' r = 3;
#' tolerance = 1e-7
#' eps = 10*r*log10(n)
#'
#' ## Generate a matrix with given data
#' U = matrix(rnorm(n*r),nrow=n)
#' V = matrix(rnorm(m*r),nrow=m)
#' Sig = diag(r)
#' M0 = U%*%Sig%*%t(V)
#'
#' ## Set some entries to be NA with probability eps/sqrt(m*n)
#' E = 1 - ceiling(matrix(rnorm(n*m),nrow=n) - eps/sqrt(m*n))
#' M_E = M0
#' M_E[(E==0)] = NA
#'
#' ## Create a noisy version
#' noiselevel = 0.1
#' M_E_noise  = M_E + matrix(rnorm(n*m),nrow=n)*noiselevel
#'
#' ## Use OptSpace for reconstruction
#' res1 = OptSpace(M_E,tol=tolerance)
#' res2 = OptSpace(M_E_noise,tol=tolerance)
#'
#' ## Compute errors for both cases using Frobenius norm
#' err_clean = norm(res1$X%*%res1$S%*%t(res1$Y)-M0,'f')/sqrt(m*n)
#' err_noise = norm(res2$X%*%res2$S%*%t(res2$Y)-M0,'f')/sqrt(m*n)
#'
#' ## print out the results
#' m1 = sprintf('RMSE without noise         : %e',err_clean)
#' m2 = sprintf('RMSE with noise of %.2f    : %e',noiselevel,err_noise)
#' print(m1)
#' print(m2)
#'
#' @export
OptSpace <- function(A,ropt=NA,niter=50,tol=1e-6,showprogress=TRUE){
  ## Preprocessing : A     : partially revelaed matrix
  if (!is.matrix(A)){
    stop("* OptSpace : an input A should be a matrix.")
  }
  if (any(is.infinite(A))){
    stop("* OptSpace : no infinite value in A is allowed.")
  }
  if (!any(is.na(A))){
    stop("* OptSpace : there is no unobserved values as NA.")
  }
  idxna = (is.na(A))
  M_E = array(0,c(nrow(A),ncol(A)))
  M_E[!idxna] = A[!idxna]
  
  ## Preprocessing : size information
  n = nrow(A)
  m = ncol(A)
  
  ## Preprocessing : other sparse-related concepts
  nnZ.E = sum(!idxna)
  E = array(0,c(nrow(A),ncol(A))); E[!idxna] = 1
  eps = nnZ.E/sqrt(m*n)
  
  ## Preprocessing : ropt  : implied rank
  if (is.na(ropt)){
    if (showprogress){
      print("* OptSpace: Guessing an implicit rank.")
    }
    r = min(max(round(guess_rank(M_E,nnZ.E)), 2), m-1)
    if (showprogress){
      print(paste0('* OptSpace: Guessing an implicit rank: Estimated rank : ',r))
    }
  } else {
    r = round(ropt)
    if ((!is.numeric(r))||(r<1)||(r>m)||(r>n)){
      stop("* OptSpace: ropt should be an integer in [1,min(nrow(A),ncol(A))].")
    }
  }
  ## Preprocessing : niter : maximum number of iterations
  if ((is.infinite(niter))||(niter<=1)||(!is.numeric(niter))){
    stop("* OptSpace: invalid niter number.")
  }
  niter = round(niter)
  
  m0 = 10000
  rho = 0
  
  ## Main Computation
  rescal_param = sqrt(nnZ.E*r/(norm(M_E,'f')^2))
  M_E = M_E*rescal_param
  
  # 1. Trimming
  if (showprogress){
    print("* OptSpace: Step 1: Trimming ...")
  }
  M_Et = M_E
  d  = colSums(E)
  d_ = mean(d)
  for (col in 1:m){
    if (sum(E[,col])>(2*d_)){
      listed = which(E[,col]>0)
      p = sample(1:length(listed),length(listed))
      M_Et[listed[p[ceiling(2*d_)]]:n,col] = 0
    }
  }
  
  d  = rowSums(E)
  d_ = mean(d)
  for (row in 1:n){
    if (sum(E[row,])>2*d_){
      listed = which(E[row,]>0)
      p = sample(1:length(listed),length(listed))
      M_Et[row,listed[p[ceiling(2*d_)]]:m] = 0
    }
  }
  
  # 2. SVD
  if (showprogress){
    print("* OptSpace: Step 2: SVD ...")
  }
  svdEt = svd(M_Et)
  X0 = svdEt$u[,1:r]
  S0 = diag(svdEt$d[1:r])
  Y0 = svdEt$v[,1:r]
  
  # 3. Initial Guess
  if (showprogress){
    print("* OptSpace: Step 3: Initial Guess ...")
  }
  X0 = X0*sqrt(n)
  Y0 = Y0*sqrt(m)
  S0 = S0/eps
  
  # 4. Gradient Descent
  if (showprogress){
    print("* OptSpace: Step 4: Gradient Descent ...")
  }
  X = X0
  Y = Y0
  S = aux_getoptS(X,Y,M_E,E)
  
  # initialize
  dist = array(0,c(1,(niter+1)))
  dist[1] = norm((M_E - (X%*%S%*%t(Y)))*E,'f')/sqrt(nnZ.E)
  for (i in 1:niter){
    # compute the gradient
    tmpgrad = aux_gradF_t(X,Y,S,M_E,E,m0,rho)
    W = tmpgrad$W
    Z = tmpgrad$Z
    
    # line search for the optimum jump length
    t = aux_getoptT(X,W,Y,Z,S,M_E,E,m0,rho)
    X = X+t*W;
    Y = Y+t*Z;
    S = aux_getoptS(X,Y,M_E,E)
    
    # compute the distortion
    dist[i+1] = norm(((M_E - X%*%S%*%t(Y))*E),'f')/sqrt(nnZ.E)
    if (showprogress){
      pmsg=sprintf('* OptSpace: Step 4: Iteration %d: distortion: %e',i,dist[i+1])
    }
    
    if (dist[i+1]<tol){
      dist = dist[1:(i+1)]
      break
    }
  }
  S = S/rescal_param
  
  
  # Return Results
  out = list()
  out$X = X
  out$S = S
  out$Y = Y
  out$dist = dist
  if (showprogress){
    print('* OptSpace: estimation finished.')
  }
  return(out)
}