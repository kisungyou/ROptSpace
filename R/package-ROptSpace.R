#' Matrix Reconstruction from a Few Entries
#'
#' Matrix reconstruction, also known as matrix completion, 
#' is the task of inferring missing entries of a partially observed matrix. 
#' This package provides a method called OptSpace, 
#' which was proposed by Keshavan, R.H., Oh, S., and Montanari, A. (2009) 
#' <doi:10.1109/ISIT.2009.5205567> for a case under low-rank assumption.
#' 
#' @noRd
#' @name ROptSpace
#' @aliases package-ROptSpace
#' @importFrom utils packageVersion
#' @import stats Rdpack
#' @importFrom Rcpp evalCpp
#' @useDynLib ROptSpace, .registration=TRUE
NULL
