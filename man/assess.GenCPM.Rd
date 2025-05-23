% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assess.GenCPM.R
\name{assess.GenCPM}
\alias{assess.GenCPM}
\title{Given a summary performance measures for the GenCPM model}
\usage{
assess.GenCPM(object, model = "linear", edge = "separate")
}
\arguments{
\item{object}{Returned GenCPM object from \code{.GenCPM} or \code{.regularized.GenCPM} functions.}

\item{model}{A character string representing one of the built-in regression models.
"linear" for \code{linear.GenCPM} and \code{linear.regularized.GenCPM};
"logistic" for \code{logit.GenCPM} and \code{logit.regularized.GenCPM};
"multinom" for \code{multinom.GenCPM} and \code{multinom.regularized.GenCPM};
and  "cox" for \code{cox.GenCPM} and \code{cox.regularized.GenCPM}.
The default is "linear".}

\item{edge}{Usage of edges to fit models, and it should be decided by the edge usage in the "object" input.
"separate" for fitting two separate models using positive edges
and negative edges respectively, and "combined" for fitting only one model
use all edges selected. The default is "separate".}
}
\value{
A list contains metrics assessing the model performance (MSE, AUC, C-index, etc), predicted response and actual response.
}
\description{
Given a summary performance measures for the GenCPM model
}
