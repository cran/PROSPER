% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S4-Klasse-Methode-PROSPER.R
\name{prosper}
\alias{prosper}
\title{Creation of a PROSPER object}
\usage{
prosper(
  Call = match.call(definition = sys.function(sys.parent(n = 2)), call =
    sys.call(sys.parent(n = 2))),
  simstruc = numeric(),
  loci = numeric(),
  simData = as.data.table()
)
}
\arguments{
\item{Call}{\code{character}, the model call of the simulation.}

\item{simstruc}{\code{numeric vector} with two numbers, the first defines the number of repetitions for the simulation, the second defines the number of simulation cycles in each repetition of the simulation.}

\item{loci}{\code{numeric}, the number of loci used (up to 4), 0 when no genetic is included.}

\item{simData}{\code{data.table}, all simulated data is stored here. The first column is treated as the number of the simulation run, the second column the number of the repetition, the third as the simulation cycle/year, the forth as the type of data (genotype based or only year based), the fifth as the genotype (\code{character}) (NA when data is year based). All following columns are up to specific the model design.}
}
\value{
An object of class "prosper".
}
\description{
The creation method for PROSPER objects.
}
\details{
The genotypes are described by their relevant loci (up to 4). Each locus can have 0, 1 or 2 alleles in a diploid genome. So a genotype for two loci is coded like that: 00 (no relevant alleles), 21 (2 alleles at the first locus, and 1 at the second). For this coding \code{loci} must be "2".
}
\examples{
\donttest{
require(data.table)
struc_preparation2(Rmx=10, af=c(0.01,0.8), epis=0, dom=c(1,0.3))
simdata <- as.data.table(dfgenotype)
   simdata[,"repetition":=1]
   simdata[,"simcycle":=1]
   setcolorder(simdata, c(c("repetition", "simcycle"),
   colnames(simdata)[!(colnames(simdata) \%in\%
   c("simulation_run", "repetition", "simcycle"))])) 
prosper(simstruc = c(5,10), loci=2, simData = simdata)  }

}
\keyword{internal}
