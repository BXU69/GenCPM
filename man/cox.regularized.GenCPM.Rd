% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cox.regularized.GenCPM.R
\name{cox.regularized.GenCPM}
\alias{cox.regularized.GenCPM}
\title{Penalized Version of Cox Model Prediction}
\usage{
cox.regularized.GenCPM(
  connectome,
  x = NULL,
  time,
  status,
  cv = "leave-one-out",
  k = dim(connectome)[3],
  thresh = 0.01,
  edge = "separate",
  type = "lasso",
  lambda = NULL,
  alpha = NULL,
  seed = 1220
)
}
\arguments{
\item{connectome}{an array indicating the connectivity between M1 edges and over N subjects. The dimension should be \code{M1*M1*N}.}

\item{x}{non-image covariates matrix, of \code{n (obs)* p (vars)}.}

\item{time}{the follow-up time for all individuals.}

\item{status}{the status indicator, normally 0=alive and 1=event.}

\item{cv}{a character indicating the method of cross-validation. The default method is "leave-one-out" cross-validation.}

\item{k}{a parameter used to set the number of folds for k-fold cross-validation.}

\item{thresh}{the value of the threshold for selecting significantly related edges. The default value is .01.}

\item{edge}{a character indicating the model is fitted with either positive and negative edges respectively or combined edges together. The default is "separate".}

\item{type}{type of penalty. “lasso” for LASSO, “ridge” for ridge, and “EN” for elastic net. The default is “lasso”.}

\item{lambda}{a user-specified lambda sequence or the optimal one automatically searched by cv.glmnet.}

\item{alpha}{the elastic net mixing parameter, ranging from 0 to 1. The default is 0.95.}

\item{seed}{the value used to set seed for random sampling in the process of cross-validation. The default value is 1220.}
}
\value{
A list contains positive edges, negative edges, predicted survival, actual survival outcomes, and the optimal lambda.
}
\description{
Penalized Version of Cox Model Prediction
}
