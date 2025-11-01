globalVariables(c("column", "Significance", "difference", "se"))

#' factorplot
#' 
#' Factorplot is a way to summarize and plot information from categorical
#' predictors from linear models and GLMs.  It creates all simple contrasts and
#' analytical standard errors for those contrasts.
#' 
#' \tabular{ll}{ Package: \tab factorplot\cr Type: \tab Package\cr Version:
#' \tab 1.2.3\cr Date: \tab 2024-05-15\cr License: \tab GPL (>=2)\cr LazyLoad:
#' \tab yes\cr } After a linear model or GLM has been estimated, the factorplot
#' command creates all pairwise differences among the levels (including the
#' reference category) of the indicated factor as well as their associated
#' standard errors to facilitate hypothesis testing directly.  The print method
#' prints the pairwise difference, standard error, p-value and
#' Bonferroni-corrected p-value.  The summary method prints the number of
#' significant positive/negative pairwise differences.  The plot method makes
#' something akin to an upper-triangular levelplot that indicates whether
#' differences are positive/negative and statistically significant.
#' 
#' @name factorplot-package
#' @aliases factorplot-package factorplot-package
#' @author Dave Armstrong Maintainer: Dave Armstrong <dave@@quantoid.net>
#' @references Armstrong, David A., II. 2013. factorplot: Improving
#' Presentation of Simple Contrasts in Generalized Linear Models.  \emph{The R
#' Journal} \bold{5(2)}: 4--15.
#' @keywords package
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' Example data for factorplot function
#' 
#' A subset of data from the 1994 Eurobarometer for France
#' 
#' 
#' @name france
#' @docType data
#' @format A data frame with 542 observations on the following 5 variables.
#' \describe{ \item{lrself}{respondent's left-right self-placement on a
#' 1(left)-10(right) scale} \item{male}{a dummy variable coded 1 for
#' males and 0 for females} \item{age}{respondent's age}
#' \item{vote}{a factor indicating vote choice with levels PCF, PS,
#' Green, RPR and UDF} \item{retnat}{a factor indicating the
#' respondent's retrospective national economic evaluation with levels Better,
#' Same and Worse} }
#' @references Reif, Karlheinz and Eric Marlier.  1997. \emph{Euro-barometer
#' 42.0: The First Year of the New European Union, November-December 1994}.
#' Inter-university Consortium for Political and Social Research (ICPSR)
#' [distributor].
#' @keywords datasets
NULL



