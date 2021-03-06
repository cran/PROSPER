% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pop_germc.R
\name{pop_germc}
\alias{pop_germc}
\title{Germination}
\usage{
pop_germc(init_sb, germ, max_vec_length = 1e+07)
}
\arguments{
\item{init_sb}{column name of the initial seed bank in the data.frame dfgenotype that is delivered by struc_preparation. \code{character}.}

\item{germ}{germination probabilities for the different cohorts. See details.}

\item{max_vec_length}{used internally, a technical term, defining the maximum length of vectors to be used.}
}
\value{
Columns are added to \code{dfgenotype}: "germ_dorm" contains the numbers of each genotype that remain dormant, "germ1" to "germX" contain the numbers for each of X cohorts.
}
\description{
\code{pop_germc} describes germination as random event for individual seeds. The function considers different cohorts and dormancy.
}
\details{
Each individual has a chance to germinate or to stay dormant. In case of germination, it emerges in one of the cohorts. The distribution of individual seeds to cohorts or dormancy is random. The function uses the columns \code{init_SB} of the data.frame \code{dfgenotype} as input. The output values are the numbers of seedlings of each genotype and each cohort.

\code{germ} must be given as a numeric vector or - in case of multiple columns in init_sb - data.frame or matrix. In that case the row number must fit to the length of \code{init_sb}. If the columns in init_sb represent cohorts, the rows of \code{germ} give the germination probabilities for these specific cohorts. The sum of one row of \code{germ} shall be \eqn{\ge 0} and \eqn{\le 1}. The difference of 1 and the sum of on row \code{germ} is the probability of dormancy.
}
\examples{
struc_preparation2(Rmx=10, af=c(0.01,0.8), epis=0, dom=1)
ls()
gen_freq( af=c(0.01,0.8), n_seeds=10000)
#Distribute the individuals to three cohorts with the germination 
#probabilities 0.2, 0.4 and 0.4.
pop_germc( init_sb="initialSB", germ=c(0.2,0.4,0.4))
rm(mf, dfgenotype, xprobab) 
}
\seealso{
\code{\link{pop_reprod}} \code{\link{pop_step}} \code{\link{quanti}}
}
