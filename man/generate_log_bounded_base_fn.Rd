% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/e_detectors.R
\name{generate_log_bounded_base_fn}
\alias{generate_log_bounded_base_fn}
\title{Generate log of baseline function for bounded RVs}
\usage{
generate_log_bounded_base_fn(
  lambda,
  m = 0.5,
  bound_lower = 0,
  bound_upper = 1,
  is_flipped = FALSE
)
}
\arguments{
\item{lambda}{Lambda parameter of the target baseline function.}

\item{m}{Mean parameter of the target baseline function. It must be strictly larger than \code{bound_lower}.}

\item{bound_lower}{Lower bound of observations. Default is \code{0}.}

\item{bound_upper}{Upper bound of observations. Default is \code{1}.}

\item{is_flipped}{A Boolean to indicate whether the model should take a flipped input or not. If the input \eqn{X} is in \eqn{[0,1]} then the flipped input is defined by \eqn{1-X}.}
}
\value{
A function compute log of baseline process given an observation.
}
\description{
Generate logarithm of the baseline function for bounded RVs in \eqn{[0,1]} given lambda and mean parameters
}
