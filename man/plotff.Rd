\name{plotff}
\alias{plotff}

\title{
This function plots objects generated by choose.fform or fform
}
\description{
plots the output from fform() and choose.fform()}
\usage{
plotff(input)
}
\arguments{
  \item{input}{
A PJ-objected produced by fform or choose.fform
}
}
\author{
Toke Emil Panduro and Cathrine Ulla Jensen
}
\examples{

library(PanJen) 
data("hvidovre")

form<-formula(lprice ~brick+roof_tile+roof_cemen+rebuild90+rebuild00+y7+y8+y9) 
PanJenArea<-fform(data=hvidovre,variable="area",base_form=form) 
plotff(PanJenArea) 
}
