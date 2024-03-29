% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelAssessment.R
\name{modelAssessment}
\alias{modelAssessment}
\title{Model assessment for LDA based on the maximum abundance for each taxon.}
\usage{
modelAssessment(
  ps,
  stan.fit,
  iterUse = 1000,
  ASVsIndexToPlot = c(1, 3, 10:14, 19:26, 36, 51:53, 148)
)
}
\arguments{
\item{ps}{phyloseq object.}

\item{stan.fit}{An instance of stanfit.}

\item{iterUse}{integer. number of iterations used in each chain after subtracting warm-up samples.}

\item{ASVsIndexToPlot}{An integer vector. Use to select ASVs to show the goodness of fit in histograms.}
}
\value{
A ggplot2 object. Histogram of data generated from the posterior estimates and the observed data.
}
\description{
Model assessment for LDA based on the maximum abundance for each taxon.
}
