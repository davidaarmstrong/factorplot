#' Calculate Pairwise Differences
#' 
#' This function calculates all pairwise difference from the input data.  The
#' input data can be the result of a GLM (produced with \code{\link[stats]{glm}}), a
#' multinomial logit model (produced with \code{multinom} from the \pkg{nnet}
#' package), a general linear hypothesis test (produced with \code{\link[multcomp]{glht}}
#' from the \pkg{multcomp} package), an object of class \code{eff} from the
#' \code{effects} package or any vector of values and a corresponding
#' variance-covariance matrix.
#' 
#' This function calculates pairwise differences that can be passed to a novel
#' plotting method that does not suffer from some of the same problems as
#' floating/quasi confidence intervals and is easier to apprehend immediately
#' than a compact letter display.
#' 
#' While the factorplot function and its print and summary methods work equally
#' well regardless of the number of levels in the \code{factor.variable}, the
#' plot function automatically scales the resulting graph to the appropriate
#' size, but will be less useful as the number of contrasts gets large (e.g., >
#' 30).  If more than one factor covariate is present and the
#' \code{factor.variable} option is NULL, the function generates a text-based
#' menu in the R GUI that will allow the users to pick the term for which they
#' want to calculate the results.
#' 
#' @aliases factorplot factorplot.lm factorplot.glm factorplot.glht
#' factorplot.summary.glht factorplot.multinom factorplot.eff
#' factorplot.default factorplot.sims
#' @param obj An object of class \code{glm} or \code{lm}, \code{glht},
#' \code{summary.glht}, \code{multinom} or a vector of values (of class
#' \code{numeric}) for which pairwise differences will be calculated.
#' Alternatively, an object of class \code{sims} which gives an Nsim x
#' Nstimulus matrix of predictions from which differences will be calculated.
#' @param factor.variable String containing the name of the factor for which
#' pairwise coefficient differences will be calculated (if a \code{glm} or
#' \code{lm} class object is passed to the function)
#' @param variable String containing the name of the column of the model matrix
#' for which pairwise differences will be calculated if a \code{multinom} class
#' object is passed to the function
#' @param var Variance-covariance matrix to be used if \code{obj} is a numeric
#' vector.  This could also be a vector of quasi/floating variances from which
#' a diagonal variance-covariance matrix will be produced
#' @param resdf Residual degrees of freedom used as the degrees of freedom for
#' the t-distribution from which p-values will be generated if \code{obj} is a
#' numeric vector
#' @param pval The (uncorrected) Type I error probability required, default =
#' 0.05
#' @param two.sided Logical argument indicating whether the hypothesis test
#' should be against a two-sided alternative if TRUE (default) or a one-sided
#' alternative if FALSE
#' @param order One of \sQuote{natural}, \sQuote{alph}, or \sQuote{size}
#' indicating how the levels of the factor should be ordered for presentation.
#' The \sQuote{natural} option (the default) leaves the levels as they are in
#' the factor contrasts.  \sQuote{alph} sorts the levels alphabetically and
#' \sQuote{size} sorts the levels by size of coefficient.
#' @param adjust.method For objects of class \code{multinom} and \code{numeric}
#' - one of the methods allowed by \code{\link[stats]{p.adjust}} in \pkg{stats} -
#' \sQuote{holm}, \sQuote{hochberg}, \sQuote{hommel}, \sQuote{bonferroni},
#' \sQuote{BH}, \sQuote{BY}, \sQuote{fdr}, \sQuote{none}.  See help for the
#' \code{\link[stats]{p.adjust}} for more information on these different adjustment
#' methods.  For objects of class \code{glm}, \code{lm} or \code{glht},
#' additional arguments of \sQuote{single-step}, \sQuote{Shaffer},
#' \sQuote{Westfall} and \sQuote{free} are possible.  See \code{\link[multcomp]{glht}}
#' from the \pkg{multcomp} package for details.
#' @param ordby For objects of class \code{eff} with interactions, \code{ordby}
#' is a string indicating the variable by which the plot should be ordered.
#' @param \dots Additional arguments to be passed to
#' \code{\link[multcomp]{summary.glht}}, including, but not limited to \code{level} and
#' \code{alternative}.
#' @return \item{b.diff}{An upper-triangular matrix of pairwise differences
#' between row and column levels of the factor} \item{b.sd}{An upper-triangular
#' matrix of standard errors of the pairwise differences represented in b.diff}
#' \item{pval}{An upper-triangular matrix of uncorrected (one-sided) p-values
#' corresponding to the entries of b.diff} \item{est}{The values of the estimates used in the calculations.}
#' \item{p}{The p-value specified in
#' the command}
#' @author Dave Armstrong
#' @references Easton, D.F., J. Peto and G.A.G. Babiker. 1991. Floating
#' absolute risk: An alternative to relative risk in survival and case control
#' analysis avoiding an arbitrary reference group.  \emph{Statistics in
#' Medicine} \bold{10}: 1025--1035.\cr Firth, David and Renee X. de Menzes.
#' 2004.  Quasi-variances.  \emph{Biometrika} \bold{91.1}: 65--80.\cr Plummer,
#' M. 2004. Improved estimates of floating absolute risk.  \emph{Statistics in
#' Medicine} \bold{23}: 93--104.\cr
#' @examples
#' 
#' ## for lm/glm
#' x <- as.factor(round(runif(1000, .5,5.5)))
#' levels(x) <- paste("lab", 1:20, sep="")
#' X <- model.matrix(~x)
#' Y <- X %*% rnorm(ncol(X),0,4) + rnorm(1000)
#' mod <- lm(Y ~ x)
#' fp <- factorplot(mod, factor.variable="x",  pval = 0.05, order="alph")
#' 
#' ## for glht
#' library(multcomp)
#' mod.glht <- glht(mod, linfct = mcp('x' = 'Tukey'))
#' fp2 <- factorplot(mod.glht, adjust.method='single-step')
#' 
#' ## for vector of values
#' b <- c(0, mod$coef[-1])
#' v <- rbind(0, cbind(0, vcov(mod)[-1,-1]))
#' names(b) <- colnames(v) <- rownames(v) <- mod$xlevels[["x"]]
#' fp3 <- factorplot(b, var=v, resdf=mod$df.residual)
#' 
#' ## for multinomial logit
#' data(france)
#' library(nnet)
#' multi.mod <- multinom(vote ~ retnat + lrself + male + age, data=france)
#' fp4 <- factorplot(multi.mod, variable="lrself")
#' 
#' @export factorplot
#' @importFrom graphics axis legend par plot polygon strwidth text
#' @importFrom stats coef contrasts model.matrix na.omit 
#' p.adjust pt terms vcov sd
#' @importFrom utils combn menu
#' @import multcomp
factorplot <- function(obj, adjust.method="none",  ...){
	UseMethod("factorplot")
}


#' @method factorplot glm
#' @rdname factorplot
#' @export
factorplot.glm <-function(obj, adjust.method="none", order="natural", factor.variable=NULL, pval=0.05, two.sided=TRUE, ...){
	tmp.classes <- attr(terms(obj), "dataClasses")
	tmp.classes <- tmp.classes[tmp.classes == "factor"]
	tmp.levs <- NULL
	for(i in 1:length(tmp.classes)){
	    tmp.levs <- c(tmp.levs, length(levels(obj$model[[names(tmp.classes)[i]]])))
	}
	tmp.class <- tmp.classes[tmp.levs > 2]
	if(is.null(factor.variable)){
	{if(length(tmp.classes) > 1){
	    myvar <- names(tmp.classes)[menu(names(tmp.classes))]
	}
	else{ 
	    myvar <- names(tmp.classes[1])
	}}
	}
	else{
	    myvar <- factor.variable
	}
	varind <- which(attr(terms(obj), "term.labels") == myvar)
	bcols <- which(attr(model.matrix(obj), "assign") == varind)
	ref <- which(apply(contrasts(obj$model[[myvar]]), 1, sum) == 0)
	b <- c(0, coef(obj)[bcols])
	names(b) <- c(levels(obj$model[[myvar]])[ref], levels(obj$model[[myvar]])[-ref])
	v <- vcov(obj)[bcols, bcols]
	v <- cbind(0, v)
	v <- rbind(0, v)
	colnames(v) <- rownames(v) <- names(b)
	levs <- obj$xlevels[[myvar]]
	if(!(order %in% c("alph", "natural", "size")))stop("Order must be one of 'size', 'natural', 'alph'")
	tmp.ord <- switch(order, 
		alph = order(names(b)), 
		size = order(b), 
		natural = 1:length(levs))
	b <- b[tmp.ord]
	v <- v[tmp.ord, tmp.ord]
	cmbn <- t(combn(length(b), 2))
	diffs <- matrix(0, nrow=length(b), ncol=nrow(cmbn))
	diffs[cbind(cmbn[,1], 1:ncol(diffs))] <- -1
	diffs[cbind(cmbn[,2], 1:ncol(diffs))] <- 1

	b.diff <- b.sd <- matrix(NA, ncol=length(b), nrow=length(b))
	b.diff[cmbn] <- b %*% diffs
	b.sd[cmbn] <- sqrt(diag(t(diffs) %*% v %*% diffs))
	colnames(b.diff) <- rownames(b.diff) <- colnames(b.sd) <- rownames(b.sd) <- names(b)
	b.diff <- b.diff[-nrow(b.diff),-1, drop=FALSE]
	b.diff <- -b.diff
	b.sd <- b.sd[-nrow(b.sd),-1, drop=FALSE]
	b.t <- b.diff/b.sd
	rns <- rownames(b.t)
	cns <- colnames(b.t)
	rns.split <- strsplit(rns, split=" ")
	rns.out <- sapply(rns.split, paste, collapse="_")
	cns.split <- strsplit(cns, split=" ")
	cns.out <- sapply(cns.split, paste, collapse="_")
	rownames(b.t) <- rownames(b.diff) <- rownames(b.sd) <- rns.out
	colnames(b.t) <- colnames(b.diff) <- colnames(b.sd) <- cns.out
	b.p <- (2^(as.numeric(two.sided)))*pt(abs(b.t), obj$df.residual,lower.tail=FALSE)
	b.bp <- array(p.adjust(b.p, method=adjust.method), dim=dim(b.p))
	ret <- list(b.diff=b.diff, b.sd=b.sd, pval = b.bp, est = b, p=pval)
	class(ret) <- c("factorplot", "list")
	ret
}

#' @method factorplot lm
#' @rdname factorplot
#' @export
factorplot.lm <-function(obj, adjust.method="none", order="natural", factor.variable=NULL, pval=0.05, two.sided=TRUE, ...){
	tmp.classes <- attr(terms(obj), "dataClasses")
	tmp.classes <- tmp.classes[tmp.classes == "factor"]
	tmp.levs <- NULL
	for(i in 1:length(tmp.classes)){
	    tmp.levs <- c(tmp.levs, length(levels(obj$model[[names(tmp.classes)[i]]])))
	}
	tmp.class <- tmp.classes[tmp.levs > 2]
	if(is.null(factor.variable)){
	{if(length(tmp.classes) > 1){
	    myvar <- names(tmp.classes)[menu(names(tmp.classes))]
	}
	else{ 
	    myvar <- names(tmp.classes[1])
	}}
	}
	else{
	    myvar <- factor.variable
	}
	varind <- which(attr(terms(obj), "term.labels") == myvar)
	bcols <- which(attr(model.matrix(obj), "assign") == varind)
	ref <- which(apply(contrasts(obj$model[[myvar]]), 1, sum) == 0)
	b <- c(0, coef(obj)[bcols])
	names(b) <- c(levels(obj$model[[myvar]])[ref], levels(obj$model[[myvar]])[-ref])
	v <- vcov(obj)[bcols, bcols]
	v <- cbind(0, v)
	v <- rbind(0, v)
	colnames(v) <- rownames(v) <- names(b)
	levs <- obj$xlevels[[myvar]]
	if(!(order %in% c("alph", "natural", "size")))stop("Order must be one of 'size', 'natural', 'alph'")
	tmp.ord <- switch(order, 
		alph = order(names(b)), 
		size = order(b), 
		natural = 1:length(levs))
	b <- b[tmp.ord]
	v <- v[tmp.ord, tmp.ord]
	cmbn <- t(combn(length(b), 2))
	diffs <- matrix(0, nrow=length(b), ncol=nrow(cmbn))
	diffs[cbind(cmbn[,1], 1:ncol(diffs))] <- -1
	diffs[cbind(cmbn[,2], 1:ncol(diffs))] <- 1

	b.diff <- b.sd <- matrix(NA, ncol=length(b), nrow=length(b))
	b.diff[cmbn] <- b %*% diffs
	b.sd[cmbn] <- sqrt(diag(t(diffs) %*% v %*% diffs))
	colnames(b.diff) <- rownames(b.diff) <- colnames(b.sd) <- rownames(b.sd) <- names(b)
	b.diff <- b.diff[-nrow(b.diff),-1, drop=FALSE]
	b.diff <- -b.diff
	b.sd <- b.sd[-nrow(b.sd),-1, drop=FALSE]
	b.t <- b.diff/b.sd
	rns <- rownames(b.t)
	cns <- colnames(b.t)
	rns.split <- strsplit(rns, split=" ")
	rns.out <- sapply(rns.split, paste, collapse="_")
	cns.split <- strsplit(cns, split=" ")
	cns.out <- sapply(cns.split, paste, collapse="_")
	rownames(b.t) <- rownames(b.diff) <- rownames(b.sd) <- rns.out
	colnames(b.t) <- colnames(b.diff) <- colnames(b.sd) <- cns.out
	b.p <- (2^(as.numeric(two.sided)))*pt(abs(b.t), obj$df.residual,lower.tail=FALSE)
	b.bp <- array(p.adjust(b.p, method=adjust.method), dim=dim(b.p))
	ret <- list(b.diff=b.diff, b.sd=b.sd, pval = b.bp, est=b,  p=pval)
	class(ret) <- c("factorplot", "list")
	ret
}

#' @method factorplot summary.glht
#' @rdname factorplot
#' @export
factorplot.summary.glht <-function(obj, ...){
	otherargs <- list(...)
	if("pval" %in% names(otherargs)){pval <- otherargs$pval}
	else{pval <- .05}
	proc.names <- do.call(rbind, strsplit(names(obj$test$coef), split=" - "))
	levs <- unique(c(proc.names[,c(2,1)]))
	b.diff <- b.sd <- b.p <- array(NA, dim=c(length(levs), length(levs)))
	colnames(b.diff) <- rownames(b.diff) <- colnames(b.sd) <- rownames(b.sd) <- colnames(b.p) <- rownames(b.p) <- levs
	proc.names <- do.call(rbind, strsplit(names(obj$test$coef), split=" - "))
	rc <- apply(proc.names, c(1,2), function(x)which(colnames(b.diff) == x))[,c(2,1)]
	b.diff[rc] <- -obj$test$coef
	b.sd[rc] <- obj$test$sigma
	b.p[rc] <- obj$test$pvalues
	b.diff <- b.diff[-nrow(b.diff),-1, drop=FALSE]
	b.sd <- b.sd[-nrow(b.sd), -1]
	b.p <- b.p[-nrow(b.p), -1]
	rns <- rownames(b.diff)
	cns <- colnames(b.diff)
	rns.split <- strsplit(rns, split=" ")
	rns.out <- sapply(rns.split, paste, collapse="_")
	cns.split <- strsplit(cns, split=" ")
	cns.out <- sapply(cns.split, paste, collapse="_")
	rownames(b.p) <- rownames(b.diff) <- rownames(b.sd) <- rns.out
	colnames(b.p) <- colnames(b.diff) <- colnames(b.sd) <- cns.out
	ret <- list(b.diff=b.diff, b.sd=b.sd, pval = b.p, p=pval)
	class(ret) <- c("factorplot", "list")
	ret
}

#' @method factorplot glht
#' @rdname factorplot
#' @export

factorplot.glht <-function(obj, adjust.method="none", pval=.05, ...){
	s.glht.obj <- summary(obj, test=adjusted(adjust.method), ...)
	ret <- factorplot(s.glht.obj)
	class(ret) <- "factorplot"
	ret
}

#' @method factorplot sims
#' @rdname factorplot
#' @export
factorplot.sims <- function(obj, adjust.method="none", order="natural", pval=.05,...){
	cmbn <- t(combn(ncol(obj), 2))
	b <- colMeans(obj)
	diffs <- matrix(0, nrow=ncol(obj), ncol=nrow(cmbn))
	diffs[cbind(cmbn[,1], 1:ncol(diffs))] <- -1
	diffs[cbind(cmbn[,2], 1:ncol(diffs))] <- 1
	sim.diffs <- obj %*% diffs
	tmp.diff <- colMeans(sim.diffs)
	tmp.sd <- apply(sim.diffs, 2, sd)
	tmp.p <- apply(sim.diffs, 2, function(x)mean(x > 0))
	tmp.p <- ifelse(tmp.p > .5, 1-tmp.p, tmp.p)
	b.diff <- b.sd <- b.p <- matrix(NA, ncol=ncol(obj), nrow=ncol(obj))
	b.diff[cmbn] <- tmp.diff
	b.sd[cmbn] <- tmp.sd
	b.p[cmbn] <- tmp.p
	colnames(b.diff) <- rownames(b.diff) <- colnames(b.sd) <- rownames(b.sd) <- colnames(b.p) <- rownames(b.p) <- colnames(obj)
	b.diff <- b.diff[-nrow(b.diff),-1, drop=FALSE]
	b.diff <- -b.diff
	b.sd <- b.sd[-nrow(b.sd),-1, drop=FALSE]
	b.p <- b.p[-nrow(b.p), -1]
	b.bp <- array(p.adjust(b.p, method=adjust.method), dim=dim(b.p))
	ret <- list(b.diff=b.diff, b.sd=b.sd, pval = b.bp, est = b, p = pval)
	class(ret) <- c("factorplot", "list")
	ret
}

#' @method factorplot default
#' @rdname factorplot
#' @export
factorplot.default <-function(obj, adjust.method="none", order="natural", var, resdf=Inf, pval=0.05, two.sided=TRUE, ...){
	b <- obj
	if(!is.matrix(var)){
		v <- diag(var)
	} else{
		v <- var
	}
	if(is.null(names(b))){
		names(b) <- as.character(1:length(b))
	}
	levs <- colnames(v) <- rownames(v) <- names(b)
	if(!(order %in% c("alph", "natural", "size")))stop("Order must be one of 'size', 'natural', 'alph'")
	tmp.ord <- switch(order, 
		alph = order(names(b)), 
		size = order(b), 
		natural = 1:length(levs))
	b <- b[tmp.ord]
	v <- v[tmp.ord,tmp.ord]
	cmbn <- t(combn(length(b), 2))
	diffs <- matrix(0, nrow=length(b), ncol=nrow(cmbn))
	diffs[cbind(cmbn[,1], 1:ncol(diffs))] <- -1
	diffs[cbind(cmbn[,2], 1:ncol(diffs))] <- 1

	b.diff <- b.sd <- matrix(NA, ncol=length(b), nrow=length(b))
	b.diff[cmbn] <- b %*% diffs
	b.sd[cmbn] <- sqrt(diag(t(diffs) %*% v %*% diffs))
	colnames(b.diff) <- rownames(b.diff) <- colnames(b.sd) <- rownames(b.sd) <- names(b)
	b.diff <- b.diff[-nrow(b.diff),-1, drop=FALSE]
	b.diff <- -b.diff
	b.sd <- b.sd[-nrow(b.sd),-1, drop=FALSE]

	b.t <- b.diff/b.sd
	rns <- rownames(b.t)
	cns <- colnames(b.t)
	rns.split <- strsplit(rns, split=" ")
	rns.out <- sapply(rns.split, paste, collapse="_")
	cns.split <- strsplit(cns, split=" ")
	cns.out <- sapply(cns.split, paste, collapse="_")
	rownames(b.t) <- rownames(b.diff) <- rownames(b.sd) <- rns.out
	colnames(b.t) <- colnames(b.diff) <- colnames(b.sd) <- cns.out
	b.p <- (2^(as.numeric(two.sided)))*pt(abs(b.t), resdf,lower.tail=FALSE)
	b.bp <- array(p.adjust(b.p, method=adjust.method), dim=dim(b.p))
	ret <- list(b.diff=b.diff, b.sd=b.sd, pval = b.bp, est = b, p = pval)
	class(ret) <- c("factorplot", "list")
	ret
}

#' @method factorplot eff
#' @rdname factorplot
#' @export
factorplot.eff <-function(obj, adjust.method="none", order="natural", pval=0.05, two.sided=TRUE, ordby = NULL,...){
	vars <- strsplit(obj$term, split="*", fixed=T)[[1]]
	b <- obj$fit
	v <- vcov(obj)
	if(ncol(obj$x) > 1){
		n <- apply(obj$x[,vars], 1, paste, collapse=":")
	}
	else{
		n <- as.character(obj$x)
	}
	names(b) <- n
	colnames(v) <- rownames(v) <- NULL
	if(!is.null(ordby)){
		if(!(ordby %in% vars))stop("Variable specifed in ordby not part of effect term")
		ord <- order(obj$x[,ordby])
		b <- b[ord]
		v <- v[ord, ord]
	}
	resdf <- nrow(obj$data)-ncol(obj$model.matrix)
	ret <- factorplot.default(b, var=v, adjust.method=adjust.method, order=order, resdf=resdf, pval=pval, two.sided=two.sided, ...)
	return(ret)
}


#' @method factorplot multinom
#' @rdname factorplot
#' @export
factorplot.multinom <- function(obj, adjust.method="none", order="natural", variable, pval = .05, two.sided=TRUE, ...){
  order <- match.arg(order)
  v <- vcov(obj)
	b <- c(t(coef(obj)))
	names(b) <- nb <- rownames(v)
  if(variable %in% obj$coefnames){
    inds <- grep(paste0("^", variable), gsub("^[^\\:]+\\:\\s*", "", names(b)))
  	v <- v[inds,inds]
  	b <- b[inds]
  	b <- c(0,b)
  	v <- rbind(0, cbind(0, v))
  	levs <- names(b) <- obj$lev
  } else if(variable %in% attr(terms(obj), "term.labels")){
      inds <- grep(paste0("^[^\\:]+\\:", variable), names(b))
      v <- v[inds,inds]
      b <- b[inds]
      b <- c(0,b)
      v <- rbind(0, cbind(0, v))
      levs <- names(b) <- c("Reference", nb[inds])
  }
		resdf <- nrow(obj$residuals) - length(c(coef(obj)))
	tmp.ord <- switch(order, 
		alph = order(names(b)), 
		size = order(b), 
		natural = 1:length(levs))
	b <- b[tmp.ord]
	v <- v[tmp.ord,tmp.ord]
	cmbn <- t(combn(length(b), 2))
	diffs <- matrix(0, nrow=length(b), ncol=nrow(cmbn))
	diffs[cbind(cmbn[,1], 1:ncol(diffs))] <- -1
	diffs[cbind(cmbn[,2], 1:ncol(diffs))] <- 1

	b.diff <- b.sd <- matrix(NA, ncol=length(b), nrow=length(b))
	b.diff[cmbn] <- b %*% diffs
	b.sd[cmbn] <- sqrt(diag(t(diffs) %*% v %*% diffs))
	colnames(b.diff) <- rownames(b.diff) <- colnames(b.sd) <- rownames(b.sd) <- names(b)
	b.diff <- b.diff[-nrow(b.diff),-1, drop=FALSE]
	b.diff <- - b.diff
	b.sd <- b.sd[-nrow(b.sd),-1, drop=FALSE]
	b.t <- b.diff/b.sd
	rns <- rownames(b.t)
	cns <- colnames(b.t)
	rns.split <- strsplit(rns, split=" ")
	rns.out <- sapply(rns.split, paste, collapse="_")
	cns.split <- strsplit(cns, split=" ")
	cns.out <- sapply(cns.split, paste, collapse="_")
	rownames(b.t) <- rownames(b.diff) <- rownames(b.sd) <- rns.out
	colnames(b.t) <- colnames(b.diff) <- colnames(b.sd) <- cns.out
	b.p <- (2^(as.numeric(two.sided)))*pt(abs(b.t), resdf,lower.tail=FALSE)
	b.bp <- array(p.adjust(b.p, method=adjust.method), dim=dim(b.p))
	ret <- list(b.diff=b.diff, b.sd=b.sd, pval = b.bp,  est = b, p = pval)
	class(ret) <- c("factorplot", "list")
	ret
}





#' Auxiliary Function to Plot a Square
#' 
#' An auxiliary function to plot squares, used by the
#' \code{\link[factorplot]{plot.factorplot}} function
#' 
#' This is a function called by \code{\link[factorplot]{plot.factorplot}} and not intended
#' to be directly used by the user; however, it is possible that this could be
#' of more general use as a utility.  The function is simply a wrapper to
#' \code{polygon} that obviates the need to specify all (\code{x},\code{y})
#' coordinates for the polygon.
#' 
#' @param ll The (\code{x},\code{y}) coordinate of the lower-left corder of the
#' square
#' @param width a scalar indicating how wide the squares should be
#' @param col a color with which the square will be filled in
#' @return \item{square}{A square is printed on the graph, but nothing else is
#' returned}
#' @author Dave Armstrong
#' @export squares
squares <- function(ll, width=1,col){ 
    poly.x <- c(ll[1], ll[1]+width, ll[1]+width, ll[1], ll[1])
    poly.y <- c(ll[2], ll[2], ll[2]+width, ll[2]+width, ll[2])
    polygon(poly.x, poly.y, col=col, border="gray85", lwd=.25)
}

#' Convert factorplot output to data frame
#' 
#' Converts the output from factorplot to a data frame that is convenient for custom plotting.  
#' The plot method for factorplot objects from the factorplot package will produce a plot that is 
#' lightly customisable.  However,  for more control over the appearance of the plot, and to plot estimates
#' on the diagnoal of the display, returning the data and making a plot is easier.  
#' 
#' @param obj An object of class `factorplot` produced by the `factorplot()` function from the `factorplot` package.
#' @param type Indicates whether you want the resulting plot to have differences on the upper triangle and p-values on the lower triangle (if `"both_tri"`) or 
#' both p-values and differences to be plotted in the upper triangle.  
#' @param ... Other arguments, currently ignored.
#' @returns A data frame with columns
#'  - `row`: The name of the parameter in the row
#'  - `col`: The name of the parameter in the column
#'  - `estimate`: The estimate for the parameter (on the diagonal)
#'  - `difference`: Pairwise differences between parameters. 
#'  - `p_value`: The p-value of the pairwise difference. 
#'  
#' @export
#' @examples
#' library(factorplot)
#' library(ggplot2)
#' data(chickwts)
#' mod <- lm(weight ~ feed, data=chickwts)
#' fp <- factorplot(mod, factor.variable = "feed", order="size")
#' chick_df <- fp_to_df(fp, type="upper_tri")
#' ggplot(chick_df, aes(x=row, y=column)) + 
#'   geom_tile(aes(fill=difference), color="black") + 
#'   geom_text(aes(label = ifelse(p_value < .05, "*", "")), color="white", size=10) + 
#'   scale_fill_viridis_c(option = "D", na.value = "transparent") + 
#'   theme_classic() + 
#'   geom_text(data=chick_df, 
#'             aes(x=row, y=column, label=round(estimate, 2))) + 
#'   labs(fill="Difference", x="", y="")
fp_to_df <- function(obj, type=c("upper_tri", "both_tri"), ...){
  typ <- match.arg(type)
  if(!"est" %in% names(obj)){
    stop("The fp_to_matrix() function does not work for factorplot calculated on glht or summary.glht objects.\n")
  }
  param_names <- c(rownames(obj$b.diff), colnames(obj$b.diff)[ncol(obj$b.diff)])
  mat <- mat_se <- matrix(NA, nrow=length(param_names), ncol = length(param_names))
  colnames(mat) <- rownames(mat) <- colnames(mat_se) <- rownames(mat_se) <- param_names
  if(typ == "both_tri"){
    mat <- mat_se <- matrix(NA, nrow=length(param_names), ncol = length(param_names))
    colnames(mat) <- rownames(mat) <- colnames(mat_se) <- rownames(mat_se) <- param_names
    mat[lower.tri(mat)] <- t(obj$b.diff)[lower.tri(obj$b.diff, diag = TRUE)]
    mat[upper.tri(mat)] <- obj$pval[upper.tri(obj$pval, diag=TRUE)]
    mat_se[upper.tri(mat)] <- obj$b.sd[upper.tri(obj$b.sd, diag=TRUE)]
    tmpl <- array(dim = dim(mat))
    diag(tmpl) <- "estimate"
    tmpl[lower.tri(tmpl)] <- "difference"
    tmpl[upper.tri(tmpl)] <- "pval"
    rownames(tmpl) <- colnames(tmpl) <- rownames(mat)
    diag(mat) <- obj$est
    mat <- mat[, ncol(mat):1]
    mat_se <- mat_se[, ncol(mat_se):1]
    tmpl <- tmpl[,ncol(tmpl):1]
    out <- as.data.frame(as.table(mat))
    out_se <- as.data.frame(as.table(mat_se))
    tmpl <- as.data.frame(as.table(tmpl))
    names(out) <- c("row", "column", "value")
    names(out_se) <- c("row", "column", "se")
    names(tmpl) <- c("row", "column", "type")
    out <- merge(out, tmpl)
    out <- merge(out, out_se)
    out$p_value <- ifelse(out$type == "pval", out$value, NA)
    out$difference <- ifelse(out$type == "difference", out$value, NA)
    out$estimate <- ifelse(out$type == "estimate", out$value, NA)
  }else{
    est_mat <- diff_mat <- pv_mat <- se_mat <- mat
    pv_mat[lower.tri(pv_mat)] <- t(obj$pval)[lower.tri(obj$pval, diag=TRUE)]
    se_mat[lower.tri(se_mat)] <- t(obj$b.sd)[lower.tri(obj$b.sd, diag=TRUE)]
    diff_mat[lower.tri(diff_mat)] <- t(obj$b.diff)[lower.tri(obj$b.diff, diag=TRUE)]
    diag(est_mat) <- obj$est
    pv_mat <- pv_mat[, ncol(pv_mat):1]
    se_mat <- se_mat[, ncol(se_mat):1]
    diff_mat <- diff_mat[, ncol(pv_mat):1]
    est_mat <- est_mat[, ncol(pv_mat):1]
    est_dat <- as.data.frame(as.table(est_mat))
    diff_dat <- as.data.frame(as.table(diff_mat))
    pv_dat <- as.data.frame(as.table(pv_mat))
    se_dat <- as.data.frame(as.table(se_mat))
    names(est_dat) <- c("row", "column", "estimate")
    names(diff_dat) <- c("row", "column", "difference")
    names(pv_dat) <- c("row", "column", "p_value")
    names(se_dat) <- c("row", "column", "se")
    out <- merge(diff_dat, pv_dat)
    out <- merge(out, est_dat)
    out <- merge(out, se_dat)
    out <- out[-which(rowSums(is.na(out)) == 4), ]
  }
  out
}


#' Plot method for objects of class factorplot
#' 
#' Creates a plot akin to an upper-triangular levelplot (though using
#' \code{plot} rather than \code{levelplot}) where the coloring of the squares
#' represents significance and text inside the squares represents the pairwise
#' difference and its correspopnding standard error.
#' 
#' 
#' @param x An object of class factorplot, produced by
#' \code{\link[factorplot]{factorplot}}.
#' @param abbrev.char The number of characters that should be used to
#' abbreviate the levels of the factor.  Set to a large value for unabbreviated
#' names.
#' @param print.est logical argument indicating whether the estimates should be
#' printed in the boxes
#' @param print.se logical argument indicating whether the standard errors
#' should be printed in the boxes
#' @param text.nudge Scalar giving the value the estimate text will be moved in the positive direction on y and the 
#' standard error text will be moved in the negative direction on y.
#' @param text.args A list of other arguments to be passed to `geom_text()`, such as `size` or `color`.
#' @param remove.caption Logical indicating whether the caption should be removed from the plot.
#' @param ... Other arguments to be passed to plot, currently not implemented
#' @return \item{a graph}{For m categories, the plot returns an m-1 x m-1 matrix
#' where the nexus of the row and column values represent the pairwise differences
#' between the row and column values along with the standard error of the difference
#' on the linear scale (unless a transformation is performed).}
#' @export
#' @author Dave Armstrong
#' @seealso \code{\link[factorplot]{factorplot}}
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual geom_text position_nudge theme_bw theme element_line element_text labs unit
#' @importFrom ggtext element_textbox_simple
#' @examples
#' 
#' est1 <- log(c(1.00,2.12,1.44,1.31,1.44,
#'      1.46,0.90))
#' var1 <- c(0.242,0.096,0.156,0.140,
#'      0.380,0.484,0.375)^2
#' names(est1) <- c(
#'      "Normal & superficial gastritis", 
#'      "Chronic gastritis", 
#'      "Chronic atrophic gastritits", 
#'      "Intestinal metaplasia I", 
#'      "Intestinal metaplasia II", 
#'      "Intestinal metaplasia III", 
#'      "Dysplasia")
#' 
#' plummer_fp1 <- factorplot(est1, var=var1, resdf=Inf)
#' plot(plummer_fp1, trans="exp", abbrev.char = 100)
#' 
#' @method plot factorplot
#' @export
plot.factorplot <- function(x, 
                            ..., 
                            abbrev.char=100, 
                            print.est = TRUE, 
                            print.se = TRUE, 
                            text.nudge = .1, 
                            text.args = list(), 
                            remove.caption = FALSE
                            ){
  d <- fp_to_df(x, type = "upper_tri")
  d$Significance <- NA
  sn_inds <- which(d$difference < 0 & d$p_value < .05)
  ns_inds <- which(d$p_value > .05)
  sp_inds <- which(d$difference > 0 & d$p_value < .05)
  if(length(sn_inds) > 0)d$Significance[sn_inds] <- "Significant Negative"
  if(length(sp_inds) > 0)d$Significance[sp_inds] <- "Significant Positive"
  if(length(ns_inds) > 0)d$Significance[ns_inds] <- "Not Significant"
  d$Significance  <- factor(d$Significance, 
                            levels = c("Significant Negative", "Not Significant", "Significant Positive"))
  levels(d$row) <- abbreviate(levels(d$row), abbrev.char)
  levels(d$column) <- abbreviate(levels(d$column), abbrev.char)
  g <- ggplot(mapping = aes(x=row, y=column, fill=Significance)) + 
    geom_tile(data = subset(d, !is.na(Significance)), 
              color="black", linewidth=.05) + 
    scale_fill_manual(values=c("Significant Negative" = "gray80", 
                                "Not Significant" = "white", 
                                "Significant Positive" = "gray30"), 
                       na.value="transparent") 
    if(print.est){
      est_text_args <- text.args
      est_text_args$fontface = "bold"
      est_text_args$position = position_nudge(y=text.nudge)
      est_text_args$data = subset(d, !is.na(difference))
      est_text_args$mapping = aes(label = sprintf("%.2f", difference))
      g <- g + do.call(geom_text, est_text_args)
    }
    if(print.se){
      se_text_args <- text.args
      se_text_args$fontface = "italic"
      se_text_args$position = position_nudge(y=-text.nudge)
      se_text_args$data = subset(d, !is.na(se))
      se_text_args$mapping = aes(label = sprintf("(%.2f)", se))
      g <- g + do.call(geom_text, se_text_args)
    }
  if(print.est & print.se){
    cap <- "Bold entries are row estimate minus column estimate, italic entries are the standard errors of the differences"
  }
  if(print.est & !print.se){
    cap <- "Text entries are row estimate minus column estimate"
  }
  if(!print.est & print.se){
    cap <- "Text entries are the standard errors of the differences"
  }
  if((!print.est & !print.se) | remove.caption){
    cap <- ""
  }
  g <- g + 
    labs(x="", y="", fill="Significance", 
         caption = cap) + 
    theme_bw() + 
    theme(panel.grid = element_line(linewidth=.2), 
        axis.text.x = element_text(angle=45, hjust=1), 
        legend.position="top", 
        plot.caption = element_textbox_simple(width = unit(1, "npc"),  # Wraps to plot width
                                              halign = 0)) 
  g
  }

#' Print method for objects of class factorplot
#' 
#' Prints the output from an object of class \code{\link[factorplot]{factorplot}}.  By
#' default, the function prints all pairwise differences along with standard
#' errors and p-values (optionally adjusted for multiple testing). Optionally,
#' it can print only significant differences.
#' 
#' 
#' @param x An object of class \code{\link[factorplot]{factorplot}}.
#' @param digits The number of digits to print in each column
#' @param trans A character string representing the post-hypothesis-testing
#' transformation to be performed on the estimates.  For example, if the
#' estimates provided to the \code{factorplot} command are log-floating
#' absolute risks, you could use the transformation \sQuote{exp}.  The
#' transformation is performed through a call to \code{do.call}
#' @param sig Logical indicating whether only significant differences should be
#' printed.
#' @param ... Other arguments passed to print, currently not implemented
#' @return \item{Printed output}{The printed output shows the difference between 
#' all pairs of stimuli (i.e., levels of the factor) along with their standard errors 
#' and (optionally adjusted) p-values.  If a transformation is implemented, the difference
#' is transformed accordingly, but the standard errors and other values are on the 
#' linear scale. }
#' @export
#' @author Dave Armstrong
#' @seealso \code{\link[factorplot]{factorplot}}
#' @examples
#' 
#' est1 <- log(c(1.00,2.12,1.44,1.31,1.44,
#'      1.46,0.90))
#' var1 <- c(0.242,0.096,0.156,0.140,
#'      0.380,0.484,0.375)^2
#' names(est1) <- c(
#'      "Normal & superficial gastritis", 
#'      "Chronic gastritis", 
#'      "Chronic atrophic gastritits", 
#'      "Intestinal metaplasia I", 
#'      "Intestinal metaplasia II", 
#'      "Intestinal metaplasia III", 
#'      "Dysplasia")
#' plummer_fp1 <- factorplot(est1, var=var1, resdf=Inf)
#' print(plummer_fp1, trans="exp")
#' @method print factorplot
#' @export
print.factorplot <- function(x, ..., digits=3, sig=FALSE, trans=NULL){
	eg <- expand.grid(rownames(x$b.diff), colnames(x$b.diff))
	mc <- apply(eg, 2, function(x)max(nchar(x)))
	strnames <- paste(sprintf("%*s", mc[1], eg[,1]), sprintf("%*s", mc[1], eg[,2]), 
		sep = " - ")
	if(!is.null(trans)){
		x$b.diff <- do.call(trans, list(x$b.diff))
	}
	tmp <- cbind(c(x$b.diff), c(x$b.sd), c(x$pval))
	tmp <- round(tmp, 3)
	rownames(tmp) <- strnames
	colnames(tmp) <- c("Difference", "SE", "p.val")
	tmp <- as.data.frame(tmp)
	tmp <- na.omit(tmp)
	if(sig == TRUE){
		if(any(tmp$p.val > ifelse("p" %in% names(x), x$p, .05))){
		tmp <- tmp[-which(tmp$p.val > ifelse("p" %in% names(x), x$p, .05)), ]
		}
	}
	t1c <- do.call(rbind, strsplit(as.character(tmp[,1]), split=".", fixed=T))
	t1mc <- apply(t1c, 2, function(x)max(nchar(x)))
	t2c <- do.call(rbind, strsplit(as.character(tmp[,2]), split=".", fixed=T))
	t2mc <- apply(t2c, 2, function(x)max(nchar(x)))

	tmp[,1] <- sprintf(paste("%*.", digits, "f", sep=""),t1mc[1], tmp[,1])
	tmp[,2] <- sprintf(paste("%*.", digits, "f", sep=""),t2mc[1], tmp[,2])
	tmp[,3] <- sprintf(paste("%1.", digits, "f", sep=""),tmp[,3])
	tmp
}



#' Summary method for objects of class factorplot
#' 
#' Summarizes the number of significant positive and negative differences for
#' objects of class \code{\link[factorplot]{factorplot}}
#' 
#' 
#' @param object An object of class \code{\link[factorplot]{factorplot}}
#' @param \dots Other arguments passed to summary, currently not implemented
#' @return \item{Printed Output}{The printed output summarises the number of 
#' stimuli that are significantly higher or lower and not significantly different
#' from each other.}
#' @export
#' @author Dave Armstrong
#' @seealso \code{\link[factorplot]{factorplot}}
#' @examples
#' 
#' x <- as.factor(round(runif(1000, .5,5.5)))
#' levels(x) <- paste("lab", 1:20, sep="")
#' X <- model.matrix(~x)
#' b <- rnorm(ncol(X),0,4)
#' Y.hat <- X %*% b 
#' Y <- Y.hat  + rnorm(1000)
#' mod <- lm(Y ~ x)
#' fp <- factorplot(mod, factor.variable="x", pval=0.05, order="alph")
#' summary(fp)
#' 
#' @method summary factorplot
#' @export
summary.factorplot <- function(object, ...){
	tmp <- object$b.diff
	tmp <- cbind(NA, tmp)
	tmp <- rbind(tmp, NA)
	rownames(tmp)[nrow(tmp)] <- colnames(tmp)[ncol(tmp)]
	colnames(tmp)[1] <- rownames(tmp)[1]
	tmp[lower.tri(tmp)] <- -t(tmp)[lower.tri(t(tmp))]
	tmp.p <- object$pval
	tmp.p <- cbind(NA, tmp.p)
	tmp.p <- rbind(tmp.p, NA)
	tmp.p[lower.tri(tmp.p)] <- t(tmp.p)[lower.tri(t(tmp.p))]
	rownames(tmp.p)[nrow(tmp.p)] <- colnames(tmp.p)[ncol(tmp.p)]
	colnames(tmp.p)[1] <- rownames(tmp.p)[1]
	diag(tmp.p) <- 1
	tmp1 <- tmp
	tmp1[which(tmp.p > object$p, arr.ind=T)] <- 0
	sig.plus <- apply(tmp1, 1, function(object)sum(object > 0))
	sig.minus <- apply(tmp1, 1, function(object)sum(object < 0))
	insig <- (nrow(tmp) -1) - (sig.plus + sig.minus)
    out <- cbind(sig.plus, sig.minus, insig)
    colnames(out) <- c("sig+", "sig-", "insig")
    out
}

