library(testthat)
library(synthesis)


nobs=1000
data.tar<-data.gen.tar(nobs,phi1=c(0.6,-0.1),phi2 = c(0.2,0.6))
data.tar2<-data.gen.tar2(nobs)

plot.ts(cbind(data.tar$x,data.tar2$x))
