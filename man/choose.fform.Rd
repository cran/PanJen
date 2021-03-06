\name{choose.fform}
\alias{choose.fform}

\title{Compare a number of user-specified  transformations with a semiparametric smoothing and a model without the variable
}

\description{
PanJen is built on the idea that the choice of a functional form can be extrapolated from model fit measures. The function provides a  ranking of different transformations according to their Baysian Information Criterion (BIC). The BIC provides a relative goodness-of-fit measure while accounting for the complexity of the model.
The function provides BIC for each transformation, a model without the variable in question and a socalled smoothing of the variable. The models are estimated a Generalized additive model (GAM).A GAM is a special case of the Generalized Linear Model (GLM), where it is possible to include one or more so called smoothing functions. A smoothing function is a non-parametric way to include a continuous independent variable in a parametric model and thus make it semi-parametric. Please see Wood (2006) for an elaboration. 
The printed output is a table with the transformations sorted according to their explanatory power measured by AIC. The table shows both AIC and BIC  for each regression, where BIC penalize for the number of explanatory variables
}

\usage{
choose.fform(data,base_form,variable,functionList, distribution=gaussian)
}
\arguments{
  \item{data}{
A data.frame
}
  \item{variable}{
A character-string with the name of the variable to test
}
  \item{base_form}{
A formula-object with the regressions without the variable that is tested
}
  \item{functionList}{
A list of transformations. Please see the example for an elaboration
}
  \item{distribution}{
Assumed distribution, see mgcv-vignette for an elaboration. Default is gaussian
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
form<-formula(lprice ~brick+roof_tile+roof_cemen+rebuild90+rebuild00+y7+y8+y9)
fxlist= list(
  linear = function(x) x,
  sqr = function(x) x^2,
  log=function(x) log(x)
)
  
PanJenAreaC<-choose.fform(data=hvidovre,variable="area",base_form=form, functionList=fxlist) 
}