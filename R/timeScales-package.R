#' timeScales: A pacakge for assessing statistical early warning signals of critical transitions at multiple time scales
#' 
#' The package uses lake experiment data to 1) test whether autocorrelation, a common early warning statistic, is sensitive to the time scale of the analysis/ observation frequency; 2) corroborate these findings using an ecological simulation; 3) propose statistics that might be time-scale-agnostic, assess whether they depend on time scales using simple statistical and more complicated ecological simulations, and then apply these statistics to the field experiment data.  
#'   
#' Some functions in this package may be useful for other applications. However, this package, and the functions and data herein, are designed as a transparent and documented record of analyses performed in the preparation of a manuscript being prepared for submission to a peer-reviewed journal. Therefore, the primary objective of this package is not for general use to assess time-scale-sensitivity of early warning statistics of arbitrary time series, but instead to test the hypotheses in the manuscript. That said, an effort was made to make many of the functions relatively general, and others may find it easy to adapt these statistics and code to their own needs.  
#' 
#' @importFrom graphics image
#' @importFrom stats ar coef density frequency
#' @importFrom utils write.csv
#'   
#' @docType package
#' @name timeScales
#' 

NULL