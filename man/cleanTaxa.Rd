% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{cleanTaxa}
\alias{cleanTaxa}
\title{Clean taxonomic names}
\usage{
cleanTaxa(x, column = 1)
}
\arguments{
\item{x}{a \code{data.frame}.}

\item{column}{the column (index or name) of \code{x} containing the taxonomic names.}
}
\value{
A \code{data.frame}.
}
\description{
This function applies a serie of text-processing steps to clean, simplify and normalize taxonomic names.
}
