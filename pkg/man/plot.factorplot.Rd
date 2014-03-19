\name{plot.factorplot}
\alias{plot.factorplot}
\title{Plot method for objects of class factorplot}
\description{Creates a plot akin to an upper-triangular levelplot (though using \code{plot} rather than \code{levelplot}) where the coloring of the squares represents significance and text inside the squares represents the pairwise difference and its correspopnding standard error. 
}
\usage{
\method{plot}{factorplot}(x, ..., abbrev.char=10, polycol=NULL, 
	textcol=NULL, trans=NULL, print.sig.leg=TRUE, print.square.leg=TRUE, 
	scale.text=1, space.text=1, print.est = TRUE, print.se=TRUE)
}
\arguments{
  \item{x}{An object of class factorplot, produced by \code{\link{factorplot}}.
}
  \item{abbrev.char}{The number of characters that should be used to abbreviate the levels of the factor.  Set to a large value for unabbreviated names.}
  \item{polycol}{A vector of three colors indicating the colors of polygons when the difference is significant negative, insignificant, and significant positive, in that order.  Defaults to c(\sQuote{gray80}, \sQuote{white}, \sQuote{gray40}).}
  \item{textcol}{A vector of three colors indicating the text color for polygons that are significant negative, insignificant, and significant positive, in that order.  Defaults to c(\sQuote{black}, \sQuote{black}, \sQuote{white})}
  \item{trans}{A character string representing the post-hypothesis-testing transformation to be performed on the estimates.  For example, if the estimates provided to the \code{factorplot} command are log-floating absolute risks, you could use the transformation \sQuote{exp}.  The transformation is performed through a call to \code{do.call}}
  \item{print.sig.leg}{logical indicating whether the legend identifying the meaning of the different colors should be included.}
  \item{print.square.leg}{logical indicating whether the legend identifying the meaning of the numbers in each square should be included.}
  \item{scale.text}{optional scale factor to be applied to text, numbers bigger than 1 make text bigger than default and numbers smaller than 1 do the opposite}
  \item{space.text}{optional text spacing factor, numbers bigger than 1 push text toward the extent of the boxes and numbers smaller than one bring text in toward the center}
  \item{print.est}{logical argument indicating whether the estimates should be printed in the boxes}
  \item{print.se}{logical argument indicating whether the standard errors should be printed in the boxes}
  \item{...}{Other arguments to be passed to plot, currently not implemented}
}
\value{
\item{a graph}{}}
\author{Dave Armstrong (Department of Political Science, UW-Milwaukee)
}
\seealso{\code{\link{factorplot}}}
\examples{
est1 <- log(c(1.00,2.12,1.44,1.31,1.44,
     1.46,0.90))
var1 <- c(0.242,0.096,0.156,0.140,
     0.380,0.484,0.375)^2
names(est1) <- c(
     "Normal & superficial gastritis", 
     "Chronic gastritis", 
     "Chronic atrophic gastritits", 
     "Intestinal metaplasia I", 
     "Intestinal metaplasia II", 
     "Intestinal metaplasia III", 
     "Dysplasia")

plummer_fp1 <- factorplot(est1, var=var1, resdf=Inf)
plot(plummer_fp1, trans="exp", abbrev.char = 100)
}