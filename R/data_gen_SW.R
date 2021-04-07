#' Generate predictor and response data: sinusoidal model
#'
#' @param nobs    The data length to be generated.
#' @param freq    The frequencies in the generated response. Default freq=50.
#' @param A       The amplitude of the sinusoidal series
#' @param phi     The phase of the sinusoidal series
#' @param mu      The mean of Gaussian noise in the variable.
#' @param sd      The standard deviation of Gaussian noise in the variable.
#'
#' @export
#' @return A list of time and x.
#'
#' @references Shumway, R. H., & Stoffer, D. S. (2011). Characteristics of Time Series. In D. S. Stoffer (Ed.), Time series analysis and its applications (pp. 8-14). New York : Springer.
#'
#' @examples
#' sample=500
#' sw <- data.gen.SW(nobs=sample, freq = 25, A = 2, phi = 0.6*pi, mu=0, sd = 0.1)
#' plot(sw$t,sw$x, type='o', ylab='Cosines', xlab="t")

data.gen.SW<-function(nobs=500,freq=50,A=2,phi=pi,mu=0,sd=1)
{

  t <- seq(0,1,length.out = nobs)

  x <- A*cos(2*pi*freq*t + phi) + rnorm(nobs,mean=mu,sd=sd)

  return(list(t=t, x=x))
}

