\encoding{UTF-8}
\name{squares}
\alias{squares}
\title{Auxiliary Function to Plot a Square}
\description{An auxiliary function to plot squares, used by the \code{\link{plot.factorplot}} function}
\usage{
squares(ll, width = 1, col)
}
\arguments{
  \item{ll}{The (\code{x},\code{y}) coordinate of the lower-left corder of the square}
  \item{width}{a scalar indicating how wide the squares should be}
  \item{col}{a color with which the square will be filled in}
}
\details{
This is a function called by \code{\link{plot.factorplot}} and not intended to be directly used by the user; however, it is possible that this could be of more general use as a utility.  The function is simply a wrapper to \code{polygon} that obviates the need to specify all (\code{x},\code{y}) coordinates for the polygon. 
}
\value{
\item{square}{A square is printed on the graph, but nothing else is returned}
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)
}
