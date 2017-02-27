
\name{fform}
\alias{fform}

\title{Compare a number of user-specified  transformations with a semiparametric smoothing and a model without the variable
}

\description{
This function runs 9 regressions, one for each of the seven predefined transformations, one without the variable and one with a smoothing of the variable. 
The printed output is a table with the transformations sorted according to their explanatory power measured by AIC. The table shows both AIC and BIC  for each regression, where BIC penalize for the number of explanatory variables.
}

\usage{
fform(data,variable,base_form, distribution=gaussian)
}
\arguments{
  \item{data}{
A data.frame
}
  \item{variable}{
A character-string with the name of the variable to test
}
  \item{base_form}{
A formula object with the regressions without the variable that is tested
}
  \item{distribution}{
Assumed distribution, see mgcv-vignette for an elaboration. Default is gaussian
}

}



\value{
 \item{rank.table }{The printed table. The transformations are ranked according to their explanatory power measured by AIC. The table shows both the AIC and BIC value, where the BIC value penalizes for the number of variables}
 \item{models }{A list of estimated models}
 \item{dataset }{A dataframe with the dataset}
}

\references{
Simon Wood, Generalized Additive Models: an introduction with R.,hapman and Hall/CRC,  2006 
}
\author{
Toke Emil Panduro & Cathrine Ulla Jensen
}

\examples{
library(PanJen) 
data("hvidovre")

form<-formula(lprice ~brick+roof_tile+roof_cemen+rebuild70+rebuild80+rebuild90+rebuild00+y7+y8+y9) 
PanJenArea<-fform(data=hvidovre,variable="area",base_form=form) 
}

