<p align="center">
  <img width="200" height="200" src="https://quantoid.net/files/images/fpsticker.png">
</p>

factorplot
==========

factorplot is an R package that helps visualize pairwise comparisons.  It is particularly useful as a post-estimation technique following a (G)LM, multinomial logistic regression or any other multiple comparison procedure done with [multcomp](https://CRAN.R-project.org/package=multcomp).  A more thorough discussion of the method and what it produces are [here](https://journal.r-project.org/archive/2013-2/armstrong.pdf):

> David A. Armstrong II.  2013.  "factorplot: Improving Presentation of Simple Contrasts in Generalized Linear Models."  The R Journal. 5(2): 4--15. 

factorplot can be installed from this repository with: 

	library(devtools)
	install_github("factorplot", username="davidaarmstrong")
