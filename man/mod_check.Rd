% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod_check.R
\name{mod_check}
\alias{mod_check}
\title{Checking the main parameters
\code{\link{struc_preparation2}}}
\usage{
mod_check(
  param.weed = NA,
  area = NA,
  duration = NA,
  repetitions = NA,
  crop_list = NA,
  max_vec_length = NA
)
}
\arguments{
\item{param.weed}{A \code{data.frame} with population dynamic parameters with or without stochasticity. The structure of \code{param.weed} is essential (see details). The easiest way to create the \code{data.frame} is to adopt an example (\link{param.ECHCG}).}

\item{area}{number of area units. \code{positive numeric}.}

\item{duration}{maximum number of simulation loops in the simulation. \code{positive integer}.}

\item{repetitions}{number of repetitions of the simulation. \code{positive integer}}

\item{crop_list}{crop rotation. \code{character vector}, elements must fit to the names in the \code{data.frame} weed.}

\item{max_vec_length}{used internally, a technical term, defining the maximum length of vectors to be used.}
}
\description{
Before the preparation (\code{\link{struc_preparation2}}) it is necessary to check, whether the necessary information is provided to principally conduct a simulation run. Furthermore the crop.list is brought to the correct form.
}
\details{
If no genetics are required \code{af} has to be set to \code{NA}. The value of \code{dom} is adjusted to the length of \code{af}. The \code{dom} is cut or the first given value is repeated to fit the number of given alleles in \code{af}. The variables af and \code{dom} are corrected and \code{n_loci} is created, which is 0 when no genetics is included.
}
\examples{
mod_check(param.weed=param.GALAP, area=100, duration=2,
     repetitions=1, crop_list=c("corn"), max_vec_length=1000000)
}
