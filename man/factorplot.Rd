\encoding{UTF-8}
\name{factorplot}
\alias{factorplot}
\alias{factorplot.lm}
\alias{factorplot.glm}
\alias{factorplot.glht}
\alias{factorplot.summary.glht}
\alias{factorplot.multinom}
\alias{factorplot.default}

\title{Calculate Pairwise Differences}
\description{
This function calculates all pairwise difference from the input data.  The input data can be the result of a GLM (produced with \code{\link{glm}}), a multinomial logit model (produced with \code{multinom} from the \pkg{nnet} package), a general linear hypothesis test (produced with \code{\link{glht}} from the \pkg{multcomp} package) or any vector of values and a corresponding variance-covariance matrix. 
}
\usage{
\method{factorplot}{glm}(obj, adjust.method="none", order="natural", 
    factor.variable=NULL, pval=.05, two.sided=TRUE, \dots)
\method{factorplot}{lm}(obj, adjust.method="none", order="natural", 
    factor.variable=NULL, pval=.05, two.sided=TRUE, \dots)
\method{factorplot}{glht}(obj, adjust.method="none", pval=.05, \dots)
\method{factorplot}{summary.glht}(obj, \dots)
\method{factorplot}{multinom}(obj, adjust.method="none", order="natural", 
    variable, pval = .05, two.sided=TRUE, \dots)
\method{factorplot}{default}(obj, adjust.method="none", order="natural", 
    var, resdf, pval=0.05, two.sided=TRUE, \dots)
}
\arguments{
  \item{obj}{An object of class \code{glm} or \code{lm}, \code{glht}, \code{summary.glht}, \code{multinom} or a vector of values (of class \code{numeric}) for which pairwise differences will be calculated.}
  \item{factor.variable}{String containing the name of the factor for which pairwise coefficient differences will be calculated (if a \code{glm} or \code{lm} class object is passed to the function)}
  \item{variable}{String containing the name of the column of the model matrix for which pairwise differences will be calculated if a \code{multinom} class object is passed to the function}
  \item{var}{Variance-covariance matrix to be used if \code{obj} is a numeric vector.  This could also be a vector of quasi/floating variances from which a diagonal variance-covariance matrix will be produced}
  \item{resdf}{Residual degrees of freedom used as the degrees of freedom for the t-distribution from which p-values will be generated if \code{obj} is a numeric vector}
  \item{pval}{The (uncorrected) Type I error probability required, default = 0.05}
  \item{two.sided}{Logical argument indicating whether the hypothesis test should be against a two-sided alternative if TRUE (default) or a one-sided alternative if FALSE}
  \item{order}{One of \sQuote{natural}, \sQuote{alph}, or \sQuote{size} indicating how the levels of the factor should be ordered for presentation.  The \sQuote{natural} option (the default) leaves the levels as they are in the factor contrasts.  \sQuote{alph} sorts the levels alphabetically and \sQuote{size} sorts the levels by size of coefficient.}
  \item{adjust.method}{For objects of class \code{multinom} and \code{numeric} - one of the methods allowed by \code{\link{p.adjust}} in \pkg{stats} - \sQuote{holm}, \sQuote{hochberg}, \sQuote{hommel}, \sQuote{bonferroni}, \sQuote{BH}, \sQuote{BY}, \sQuote{fdr}, \sQuote{none}.  See help for the \code{\link{p.adjust}} for more information on these different adjustment methods.  For objects of class \code{glm}, \code{lm} or \code{glht}, additional arguments of \sQuote{single-step}, \sQuote{Shaffer}, \sQuote{Westfall} and \sQuote{free} are possible.  See \code{\link{glht}} from the \pkg{multcomp} package for details. }	
  \item{\dots}{Additional arguments to be passed to \code{\link{summary.glht}}, including, but not limited to \code{level} and \code{alternative}.}
}
\details{
This function calculates pairwise differences that can be passed to a novel plotting method that does not suffer from some of the same problems as floating/quasi confidence intervals and is easier to apprehend immediately than a compact letter display.  

While the factorplot function and its print and summary methods work equally well regardless of the number of levels in the \code{factor.variable}, the plot function automatically scales the resulting graph to the appropriate size, but will be less useful as the number of contrasts gets large (e.g., > 30).  If more than one factor covariate is present and the \code{factor.variable} option is NULL, the function generates a text-based menu in the R GUI that will allow the users to pick the term for which they want to calculate the results.    
}
\value{
\item{b.diff}{An upper-triangular matrix of pairwise differences between row and column levels of the factor}
\item{b.sd}{An upper-triangular matrix of standard errors of the pairwise differences represented in b.diff}
\item{pval}{An upper-triangular matrix of uncorrected (one-sided) p-values corresponding to the entries of b.diff}
\item{p}{The p-value specified in the command}
}
\references{
Easton, D.F., J. Peto and G.A.G. Babiker. 1991. Floating absolute risk: An alternative to relative risk in survival and case control analysis avoiding an arbitrary reference group.  \emph{Statistics in Medicine} \bold{10}: 1025--1035.\cr
Firth, David and Renee X. de Menzes.  2004.  Quasi-variances.  \emph{Biometrika} \bold{91.1}: 65--80.\cr
Plummer, M. 2004. Improved estimates of floating absolute risk.  \emph{Statistics in Medicine} \bold{23}: 93--104.\cr}
\author{Dave Armstrong (Department of Political Science, UW-Milwaukee)}
\examples{
## for lm/glm
x <- as.factor(round(runif(1000, .5,5.5)))
levels(x) <- paste("lab", 1:20, sep="")
X <- model.matrix(~x)
Y <- X \%*\% rnorm(ncol(X),0,4) + rnorm(1000)
mod <- lm(Y ~ x)
fp <- factorplot(mod, factor.variable="x",  pval = 0.05, order="alph")

## for glht
library(multcomp)
mod.glht <- glht(mod, linfct = mcp('x' = 'Tukey'))
fp2 <- factorplot(mod.glht, adjust.method='single-step')

## for vector of values
b <- c(0, mod$coef[-1])
v <- rbind(0, cbind(0, vcov(mod)[-1,-1]))
names(b) <- colnames(v) <- rownames(v) <- mod$xlevels[["x"]]
fp3 <- factorplot(b, var=v, resdf=mod$df.residual)

## for multinomial logit
data(france)
library(nnet)
multi.mod <- multinom(vote ~ retnat + lrself + male + age, data=france)
fp4 <- factorplot(multi.mod, variable="lrself")
}

