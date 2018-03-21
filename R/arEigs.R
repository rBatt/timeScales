
#' Extend Diagonal
#' 
#' Place 1's on diagonal, or on superdiagonal 'offset' rows above diagonal
#' 
#' @param size scalar, size of matrix (nrow, nrow=ncol)
#' @param offset scalar, -1 puts 1's one row below diagonal, 1 puts 1's one row above diagonal, -2 two rows below, etc
#' 
#' @return a matrix with 1's on the 'offset' of the diagonal, 0's everywhere else
#' @export
diagExtend <- function(size, offset=-1){ 
	M <- matrix(0, size, size)
	M[row(M)+offset == col(M)] <- 1
	return(M)
}

#' Eigenvalues
#' 
#' @param b numeric vector of AR coefficients
#' 
#' @return eigenvalues
#' 
#' @seealso \code{\link{diagExtend}}, \code{\link{ar}}
#' @export
arEigs <- function(b){
	nb <- length(b)
	B <- diagExtend(nb,-1)
	B[1,1:nb] <- b
	eigs <- eigen(B)$values
	# eigB <- max(abs(eigs))
	return(eigs)
}