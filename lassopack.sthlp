{smcl}
{* *! version 1.0.12  31july2020}{...}
{hline}
{cmd:help lassopack}{right: lassopack v1.4.1}
{hline}

{title:Package}

{p2colset 5 16 18 2}{...}
{p2col:{hi:LASSOPACK}}{p_end}
{p2colreset}{...}

{title:Overview}

{pstd}
LASSOPACK is a suite of programs for penalized regression methods suitable 
for the high-dimensional setting where the number of predictors p may be 
large and possibly greater than the number of observations.

{pstd}
The package consists of six main programs: 

{pstd}
{helpb lasso2} implements lasso, square-root lasso, elastic net, ridge regression,
adaptive lasso and post-estimation OLS. 
The lasso (Least Absolute Shrinkage and Selection Operator, Tibshirani 1996), 
the square-root-lasso (Belloni et al. 2011) and the adaptive lasso (Zou 2006) 
are regularization methods that use L1 norm penalization to achieve sparse solutions: 
of the full set of p predictors, typically most will have coefficients set to 
Ridge regression (Hoerl & Kennard 1970) relies on L2 norm penalization; 
the elastic net (Zou & Hastie 2005) uses a mix of L1 and L2 penalization. 

{pstd}
{helpb cvlasso} supports K-fold cross-validation and rolling cross-validation 
for cross-section, panel and time-series data. 

{pstd}
{helpb rlasso} implements theory-driven penalization for the lasso and square-root 
lasso for cross-section and panel data. 
rlasso uses the theory-driven penalization methodology of 
Belloni et al. (2012, 2013, 2014, 2016) for the lasso and square-root lasso.

{pstd}
{helpb lassologit}, {helpb cvlassologit} and {helpb rlassologit} are the 
corresponding programs for logistic lasso regression. 

{pstd}
For more information, please see our website
{browse "https://statalasso.github.io/"}, 
the help files and our paper below.

{title:Citation of lassopack}

{pstd}{opt lassopack} is not an official Stata package. It is a free contribution
to the research community, like a paper. Please cite it as such: {p_end}

{phang}Ahrens, A., Hansen, C.B., Schaffer, M.E. 2018 (updated 2020).
LASSOPACK: Stata module for lasso, square-root lasso, elastic net, ridge, adaptive lasso estimation and cross-validation
{browse "http://ideas.repec.org/c/boc/bocode/s458458.html"}{p_end}

{phang}Ahrens, A., Hansen, C.B. and M.E. Schaffer. 2020.
lassopack: model selection and prediction with regularized regression in Stata.
{it:The Stata Journal}, 20(1):176-235.
{browse "https://journals.sagepub.com/doi/abs/10.1177/1536867X20909697"}.
Working paper version: {browse "https://arxiv.org/abs/1901.05397"}.{p_end}

{title:Authors}

	Achim Ahrens, Public Policy Group, ETH Zurich, Switzerland
	achim.ahrens@gess.ethz.ch
	
	Christian B. Hansen, University of Chicago, USA
	Christian.Hansen@chicagobooth.edu

	Mark E. Schaffer, Heriot-Watt University, UK
	m.e.schaffer@hw.ac.uk

