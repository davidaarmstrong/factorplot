\name{summary.factorplot}
\alias{summary.factorplot}
\title{Summary method for objects of class factorplot
}
\description{Summarizes the number of significant positive and negative differences for objects of class \code{\link{factorplot}}
}
\usage{
	\method{summary}{factorplot}(object, \dots)
}
\arguments{
  \item{object}{An object of class \code{\link{factorplot}}}
  \item{\dots}{Other arguments passed to summary, currently not implemented}
}
\author{Dave Armstrong (Department of Political Science, UW-Milwaukee)
}
\seealso{\code{\link{factorplot}}}
\examples{
x <- as.factor(round(runif(1000, .5,5.5)))
levels(x) <- paste("lab", 1:20, sep="")
X <- model.matrix(~x)
b <- rnorm(ncol(X),0,4)
Y.hat <- X \%*\% b 
Y <- Y.hat  + rnorm(1000)
mod <- lm(Y ~ x)
fp <- factorplot(mod, factor.variable="x", pval=0.05, order="alph")
summary(fp)
}
