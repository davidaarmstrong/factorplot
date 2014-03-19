\name{print.factorplot}
\alias{print.factorplot}
\title{Print method for objects of class factorplot}
\description{Prints the output from an object of class \code{\link{factorplot}}.  By default, the function prints all pairwise differences along with standard errors and p-values (optionally adjusted for multiple testing). Optionally, it can print only significant differences.}
\usage{
\method{print}{factorplot}(x, ..., digits = 3, sig = FALSE, trans=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link{factorplot}}. 
}
  \item{digits}{The number of digits to print in each column
}
\item{trans}{A character string representing the post-hypothesis-testing transformation to be performed on the estimates.  For example, if the estimates provided to the \code{factorplot} command are log-floating absolute risks, you could use the transformation \sQuote{exp}.  The transformation is performed through a call to \code{do.call}}
  \item{sig}{Logical indicating whether only significant differences should be printed.}
  \item{...}{Other arguments passed to print, currently not implemented}
}
\author{Dave Armstrong (Department of Political Science, UW-Milwaukee)
}
\seealso{\code{\link{factorplot}}
}
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
print(plummer_fp1, trans="exp")}
