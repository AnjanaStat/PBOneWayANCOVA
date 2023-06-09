% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCI1.R
\name{SCI1}
\alias{SCI1}
\title{Find simultaneous confidence intervaks for the successive pair-wise differences of treatment effects in one-way ANCOVA model with single covariate}
\usage{
SCI1(Y, X, k, alpha)
}
\arguments{
\item{Y}{a real dataset}

\item{X}{a real dataset}

\item{k}{a positive integer}

\item{alpha}{real number between 0 and 1 called significance level}
}
\value{
numeric values
}
\description{
More detailed description
}
\examples{
k=4;N=c(20,30,20,50);S=c(1,1,1,2);alpha=0.05
g=NULL
for(i in 1:k)
{
 g[[i]]=rnorm(N[i],2,sqrt(S[i]))
}
Y=g
G=NULL
for(i in 1:k)
{
G[[i]]=rnorm(N[i],3,sqrt(S[i]))
}
X=G
SCI1(Y,X,4,0.05)
}
