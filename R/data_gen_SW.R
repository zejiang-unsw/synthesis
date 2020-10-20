#' Generate predictor and response data: sinusoidal model
#'
#' @param nobs    The data length to be generated.
#' @param fp      The frequencies in the generated response. Default fp=50.
#' @param fd      A vector of frequencies for potential predictors. Default fd=50.
#' @param sd      The noise level in the predictor.
#'
#' @return A list of 3 elements: a vector of response (x), a matrix of potential predictors (dp) with each column containing one potential predictor, and a vector of true predictor numbers.
#'
#' @references Shumway, R. H., & Stoffer, D. S. (2011). Characteristics of Time Series. In D. S. Stoffer (Ed.), Time series analysis and its applications (pp. 8-14). New York : Springer.
#'
data.gen.SW<-function(nobs=512,fp=50,fd=50,sd=1)
{

  t <- 1:nobs

  index <- which(fd %in% fp)
  ndim <- length(fd)

  dp<-matrix(0,nobs,ndim)
  for(i in 1:ndim){
    eps<-rnorm(nobs,mean=0,sd=sd)
    dp[,i]<- 2*cos(2*pi*t/fd[i] + 0.6*pi)

  }

  x <- 2*cos(2*pi*t/fp + 0.6*pi) + rnorm(nobs,mean=0,sd=sd)

  data_generated<-list(x=x,
                       dp=dp,
                       true.cpy=index)

  return(data_generated)
}

#' Generate predictor and response data: sinusoidal model
#'
#' @param nobs    The data length to be generated.
#' @param fp      The frequencies in the generated response.
#' @param fd      A vector of frequencies for potential predictors. fd = c(3,5,10,15,25,30,55,70,95) used in the WRR paper.
#' @param sd.x    The noise level in the predictor.
#' @param sd.y    The noise level in the response.
#'
#' @return A list of 3 elements: a vector of response (x), a matrix of potential predictors (dp) with each column containing one potential predictor, and a vector of true predictor numbers.
#'
data.gen.SW1<-function(nobs=512,fp=25,fd,sd.x=0.1,sd.y=0.1)
{

  t <- seq(0,1,length.out = nobs)

  index <- which(fd %in% fp)
  ndim <- length(fd)

  dp<-matrix(0,nobs,ndim)
  for(i in 1:ndim){
    eps<-rnorm(nobs,mean=0,sd=1)
    dp[,i]<- sin(2*pi*fd[i]*t) + sd.x*eps

  }

  x <- rowSums(sapply(1:length(index), function(i) sin(2*pi*fd[index[i]]*t))) + rnorm(nobs,0,sd.y)

  data_generated<-list(x=x,
                       dp=dp,
                       true.cpy=index)

  return(data_generated)
}
