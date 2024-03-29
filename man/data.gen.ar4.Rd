% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_gen_AR.R
\name{data.gen.ar4}
\alias{data.gen.ar4}
\title{Generate predictor and response data from AR4 model.}
\usage{
data.gen.ar4(nobs, ndim = 9)
}
\arguments{
\item{nobs}{The data length to be generated.}

\item{ndim}{The number of potential predictors (default is 9).}
}
\value{
A list of 2 elements: a vector of response (x), and a matrix of potential predictors (dp) with each column containing one potential predictor.
}
\description{
Generate predictor and response data from AR4 model.
}
\examples{
# AR4 model from paper with total 9 dimensions
data.ar4<-data.gen.ar4(500)
plot.ts(cbind(data.ar4$x,data.ar4$dp))

}
