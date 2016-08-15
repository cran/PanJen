\name{choose.fform}
\alias{choose.fform}

\title{Compare a number of user-specified  transformations with a semiparametric smoothing and a model without the variable
}
\description{
This function runs a regression for each user-defined transformations plus  one without the variable and one with a smoothing of the variable. 
The printed output is a table with the transformations sorted according to their explanatory power measured by AIC. The table shows both AIC and BIC  for each regression, where BIC penalize for the number of explanatory variables
}
\usage{
choose.fform(data,variable,forms_to_test,distribution=gaussian,smoothing_splines=20)
}
\arguments{
  \item{data}{
A data.frame
}
  \item{variable}{
A character-string with the name of the variable to test
}

  \item{forms_to_test}{
A formula object with the regressions without the variable that is tested and each transformation to test seperated by "|") Please see the example for an elaboration
}
  \item{distribution}{
Assumed distribution, see mgcv-vignette for an elaboration. Default is gaussian
}
  \item{smoothing_splines}{
Number of splines, see mgcv-vignette or Wood(2006) for an elaboration, Default is 20
}
}



\value{
 \item{rank.table }{The printed table. The transformations are ranked according to their explanatory power measured by AIC. The table shows both the AIC and BIC value, where the BIC value penalizes for the number of variables}
 \item{models }{A list of estimated models}
 \item{dataset }{A dataframe with the dataset}
 \item{fforms }{The formula provided by the user}
}

\references{
Simon Wood, Generalized Additive Models: an introduction with R.,hapman and Hall/CRC,  2006 
}
\author{
Toke Emil Panduro & Cathrine Ulla Jensen
}
\examples{
## Test a linear specification (x), a log (I(log(x))) and a square I(x^2))
library(PanJen) 
data("hvidovre")
form<-formula(lprice ~brick+roof_tile+roof_cemen+rebuild90+rebuild00+y7+y8+y9|x|I(log(x))|I(x^2))
PanJenAreaC<-choose.fform(hvidovre,"area",form) 
}