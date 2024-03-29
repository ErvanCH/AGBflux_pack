% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AGBfluxes_functions.R
\name{loss.flux}
\alias{loss.flux}
\title{Unbiased loss flux}
\usage{
loss.flux(A0, A1, S1, area, time)
}
\arguments{
\item{A0}{biomass at initial census.}

\item{A1}{biomass of alive trees at final census.}

\item{S1}{initial biomass of individuals that survived to time t}

\item{time}{cenusus interval in year}
}
\value{
absolute recruitment flux in % per year
}
\description{
Compute unbiased recruitment rate (i.e. account for unmeasured recruitment)
}
\author{
Helene Muller-Laudau (hmullerlandau@gmail.com)
}
