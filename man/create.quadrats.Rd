% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AGBfluxes_functions.R
\name{create.quadrats}
\alias{create.quadrats}
\title{Create quadrats}
\usage{
create.quadrats(census, grid_size, x = "gx", y = "gy")
}
\arguments{
\item{census}{a data.frame where trees have relative X-Y coordinates.}

\item{grid_size}{the size of the grid (in meter)}

\item{x}{the identifier of X coordinates (i.e. 'gx')}

\item{y}{the identifier of Y coordinates (i.e. 'gy')}
}
\value{
add three columns to the data.frame: quadrat's number, centroids X and Y.
}
\description{
Creat a grid where all trees are allocated to a given quadrat of size (= grid size).
}
\author{
Ervan Rutishauser (ervan.rutishauser@infoflora.ch)
}
