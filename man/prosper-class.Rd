% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S4-Klasse-Methode-PROSPER.R
\docType{class}
\name{prosper-class}
\alias{prosper-class}
\title{Class for PROSPER results}
\description{
This \code{class} simplifies the analysis of PROSPER results, by providing a fixed structure with some basic summary and plot methods.
}
\section{Slots}{

\describe{
\item{\code{Call}}{\code{character}, the model call of the simulation.}

\item{\code{simstruc}}{\code{numeric vector} with 2 numbers, the first defines the number of repetitions for the simulation, the second defines the number of simulation cycles in each repetition of the simulation.}

\item{\code{loci}}{\code{numeric}, the number of loci used, 0 when no genetic is included.}

\item{\code{simData}}{\code{data.table}, all simulated data is stored here. The first column is treated as the number of the simulation run, the second column the number of the repetition, the third as the simulation cycle/year, the forth as the type of data (genotype based or only year based), the fifth as the genotype (\code{character}) (NA when data is year based). All following columns are up to specific the model design.}
}}

\keyword{internal}
