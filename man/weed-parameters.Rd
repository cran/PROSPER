% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weed-parameters.R
\docType{data}
\name{weed-parameters}
\alias{weed-parameters}
\title{Parameterisations for population dynamic models - Examples}
\format{
The provided data are given in tables with 10 columns.
\describe{
\item{weed}{the weed species.} 
\item{crop}{the crop for which the parameters apply, \code{character}.}
\item{variable}{the name of the variable, simulation step, \code{character}.}
\item{name}{a naming of parameter, like a variable definition, \code{character}.}
\item{mean}{the mean of the parameters defined with "name", \code{numeric}.}
\item{se}{the standard error of the parameters defined with "name", \code{numeric}.}
\item{model}{a (statistical) model to describe a specific simulation step. The model must be constructed of defined parameters with "name", \code{character}.}
\item{range_low}{a model might have a minimum level of use, then this value can be defined here, \code{numeric}.}
\item{range_up}{a model might have a maximum level of use, then this value can be defined here, \code{numeric}.}
\item{source}{if the value is taken from sources not owned by the author they can be named here.}
}
}
\description{
PROSPER entails full parameterized models. The parameterisation for these models are described here.
}
\details{
\enumerate{
\item \code{param.LOLRI} entails the data used for the introduction of PERTH by Renton et al. (2011) using the parameters of Pannell et al. (2004). With \code{prosper.LOLRI} one of the PERTH models can be reproduced.
\item \code{param.ECHCG} contains parameters for population dynamic of \emph{Echinochloa crus-galli} (Redwitz et al., 2016).
\item \code{param.GALAP} contains parameters for population dynamic off \emph{Galium aparine} taken from field experiments (Daedlow, 2015). 
}
}
\references{
\enumerate{
\item Bagavathiannan, MV, Norsworthy, JK (2013) Postdispersal loss of important arable weed seeds in the midsouthern United States. Weed Science, 61(4), 570-579.
\item Bosnic, AC, Swanton, C (1997) Influence of Barnyardgrass (Echinochloa crus-galli) Time of Emergence and Density on Corn (Zea mays). Weed Science, 45(2), 276-282.
\item Clay, SA, Kleinjan, J, Clay, DE, Forcella, Frank, Batchelor, W (2005) Growth and fecundity of several weed species in corn and soybean. Agronomy journal, 97(1), 294-302.
\item Daedlow, D. (2015): About the contribution of seed predation on weed demography. (Doctoral dissertation). University of Rostock, Faculty of Agriculture and Environmental Sciences. GBV Gemeinsamer Verbundkatalog (Accession No. 839752644).
\item Ogg, AG, Dawson, JH (1984) Time of emergence of eight weed species. Weed Science, 32(3), 327-335.
\item Pannell, DJ, Stewart, V, Bennett, A, Monjardino, M, Schmidt, C, Powles, SB (2004) RIM:a bioeconomic model for integrated weedmanagement of Lolium rigidum in Western Australia. Agricultural Systems, 79, 305-325.
\item Redwitz, C von, Pannwitt, H, Gerowitt, B (2016): About the interplay of sensitive and resistant biotypes in weed populations - simulation exercises for Echinochloa crus-galli in maize crops. Julius-Kuehn-Archiv 452, 93-99.
\item Renton, M.; Diggle, A.; Manalil, S. & Powles, S. (2011): Does cutting herbicide rates threaten the sustainability of weed management in cropping systems? Journal of Theoretical Biology, 283, 14-27.
}
}
\seealso{
\code{\link{prosper-models}} \code{\link{quanti}}
}
