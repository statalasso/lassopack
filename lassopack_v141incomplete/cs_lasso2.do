* certification script for 
* lassopack package 1.4.X 18aug2020, aa/ms
* parts of the script use R's glmnet and Matlab code "SqrtLassoIterative.m".

set more off
cscript "lasso2" adofile lasso2 lasso2_p lassoutils
clear all
capture log close
set rmsg on
program drop _all
log using cs_lasso2,replace
about
which lasso2
which lasso2_p
which lassoutils

* data source
//global prostate prostate.data
global prostate https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data

* simple ridge regression program
cap program drop estridge
program define estridge, rclass
	syntax varlist , Lambda(real) [NOCONStant]
	local yvar		: word 1 of `varlist'
	local xvars		: list varlist - yvar
	qui putmata X=(`xvars') y=(`yvar'), replace
	mata: n=rows(y)
	if ("`noconstant'"=="") {
		mata: X=X:-mean(X)
		mata: y=y:-mean(y)
	}
	mata: p=cols(X)
	mata: beta=lusolve(X'X+(`lambda')/2*I(p),X'y)
	tempname bhat
	mata: st_matrix("`bhat'",beta')
	mat list `bhat'
	return matrix bhat = `bhat'
end

cap program drop comparemat
program define comparemat , rclass
	syntax anything [, tol(real 10e-3)] 
	local A		: word 1 of `0'
	local B		: word 2 of `0'
	tempname Amat Bmat
	mat `Amat' = `A'
	mat `Bmat' = `B'
	local diff=mreldif(`Amat',`Bmat')
	di as text "mreldif=`diff'. tolerance = `tol'"
	mat list `Amat'
	mat list `Bmat'
	return scalar mreldif = `diff'
	assert `diff'<`tol'
end

* program to compare two vectors using col names
cap program drop comparevec
program define comparevec , rclass
	syntax anything [, tol(real 10e-3)] 
	local A		: word 1 of `0'
	local B		: word 2 of `0'
	tempname Amat Bmat
	mat `Amat' = `A'
	mat `Bmat' = `B'
	local Anames: colnames `Amat' 
	local Bnames: colnames `Bmat'
	local maxdiff = 0
	local num = 0
	foreach var of local Anames {
		local aix = colnumb(`Amat',"`var'")
		local bix = colnumb(`Bmat',"`var'")
		//di `aix'
		//di `bix'
		local thisdiff=reldif(el(`Amat',1,`aix'),el(`Bmat',1,`bix'))
		if `thisdiff'>`maxdiff' {
			local diff = `thisdiff'
		}
		local num=`num'+1
	}
	di as text "Max rel dif = `maxdiff'. tolerance = `tol'"
	mat list `Amat'
	mat list `Bmat'
	return scalar maxdiff = `maxdiff'
	assert `maxdiff'<`tol'
end


********************************************************************************
*** replicate glmnet														 ***
********************************************************************************

sysuse auto, clear
drop if rep78==.

global model price mpg-foreign

// # the following R code was run using ‘glmnet’ version 4.0-2
/*
library("glmnet")
library("haven")
library("tidyverse")

auto <- read_dta("http://www.stata-press.com/data/r9/auto.dta")

auto <- auto %>% drop_na()
n <- nrow(auto)

price <- auto$price

X <- auto[,c("mpg","rep78","headroom","trunk",
             "weight","length","turn",
             "displacement","gear_ratio","foreign")]
X <- X %>% 
  mutate(foreign = as.integer(foreign)) %>%
  as.matrix()
*/

// single lambda, lasso
/*
> r<-glmnet(X,price,alpha=1,lambda=1000,thresh=1e-15)
> t(coef(r))
1 x 11 sparse Matrix of class "dgCMatrix"
   [[ suppressing 11 column names ‘(Intercept)’, ‘mpg’, ‘rep78’ ... ]]
                                             
s0 4336.84 . . . . 0.3819041 . . 3.289185 . .
*/
mat G = 0.3819041, 3.289185, 4336.84
lasso2 $model, lambda(1000) lglmnet
assert mreldif(e(b),G) <1e-5

// single lambda, ridge
/*
> r<-glmnet(X,price,alpha=0,lambda=1000,thresh=1e-15)
> t(coef(r))
1 x 11 sparse Matrix of class "dgCMatrix"
   [[ suppressing 11 column names ‘(Intercept)’, ‘mpg’, ‘rep78’ ... ]]
                                                
s0 3425.776 -57.0062 296.9702 -392.4219 23.70664
                                                 
s0 0.9938447 5.562674 -24.18045 8.855065 -536.018
           
s0 1726.658
*/
mat G = -57.0062, 296.9702, -392.4219, 23.70664, 0.9938447, 5.562674,	///
			-24.18045, 8.855065, -536.018, 1726.65, 3425.776
lasso2 $model, alpha(0) lambda(1000) lglmnet
assert mreldif(e(b),G) <1e-5

// single lambda, elastic net
/*
> r<-glmnet(X,price,alpha=0.5,lambda=1000,thresh=1e-15)
> t(coef(r))
1 x 11 sparse Matrix of class "dgCMatrix"
   [[ suppressing 11 column names ‘(Intercept)’, ‘mpg’, ‘rep78’ ... ]]
                                                   
s0 2883.129 -5.096381 . . . 0.692915 . . 5.942777 .
          
s0 308.225
*/
mat G = -5.096381, 0.692915, 5.942777, 308.225, 2883.129
lasso2 $model, alpha(0.5) lambda(1000) lglmnet
assert mreldif(e(b),G) <1e-5

// lasso, lambda grid
/*
> r<-glmnet(X,price,alpha=1,nlambda=5,thres=1e-15)
> r$lambda
[1] 1584.1894208  158.4189421   15.8418942    1.5841894
[5]    0.1584189
> r$dev.ratio
[1] 0.0000000 0.5272556 0.5970093 0.5988338 0.5988521
*/
mat L =  1584.1894208, 158.4189421, 15.8418942, 1.5841894, 0.1584189
mat L = L'
mat D = 0.0000000, 0.5272556, 0.5970093, 0.5988338, 0.5988521
mat D = D'
lasso2 $model, lglmnet lcount(5)
assert mreldif(e(lambdamat0),L) < 1e-5
assert mreldif(e(rsq),D) < 1e-5

// ridge, lambda grid
/*
> r<-glmnet(X,price,alpha=0,nlambda=5,thres=1e-15)
> r$lambda
[1] 1584189.4208  158418.9421   15841.8942    1584.1894
[5]     158.4189
> r$dev.ratio
[1] 2.777378e-36 4.347508e-02 2.060125e-01 4.453574e-01
[5] 5.744210e-01
*/
mat L =  1584189.4208, 158418.9421, 15841.8942, 1584.1894, 158.4189
mat L = L'
// first R-sq of glmnet appears to be wrong - see below
* mat D = 2.777378e-36, 4.347508e-02, 2.060125e-01, 4.453574e-01, 5.744210e-01
* mat D = D'
lasso2 $model, alpha(0) lglmnet lcount(5) long
assert mreldif(e(lambdamat0),L) < 1e-5
* assert mreldif(e(rsq),D) < 1e-5

// single lambda, ridge - see above
/*
> r<-glmnet(X,price,alpha=0,lambda=1584189.4208,thresh=1e-15)
> r$dev.ratio
[1] 0.004941719
*/
lasso2 $model, alpha(0) lglmnet lambda(1584189.4208)
assert reldif(e(r2),0.004941719) < 1e-5

// elastic net, grid of 5
/*
> r<-glmnet(X,price,alpha=0.5,nlambda=5,thresh=1e-15)
> r$lambda
[1] 3168.3788416  316.8378842   31.6837884    3.1683788
[5]    0.3168379
> r$dev.ratio
[1] 0.0000000 0.5088361 0.5942304 0.5987942 0.5988517
*/
mat L = 3168.3788416, 316.8378842, 31.6837884, 3.1683788, 0.3168379
mat L = L'
mat D = 0.0000000, 0.5088361, 0.5942304, 0.5987942, 0.5988517
mat D = D'
lasso2 $model, alpha(0.5) lglmnet lcount(5) long
assert mreldif(e(lambdamat0),L) < 1e-5
assert mreldif(e(rsq),D) < 1e-5

// single lambda, lasso, nocons
/*
> r<-glmnet(X,price,alpha=1,lambda=1000,intercept=FALSE,thresh=1e-15)
> t(coef(r))
1 x 11 sparse Matrix of class "dgCMatrix"
   [[ suppressing 11 column names ‘(Intercept)’, ‘mpg’, ‘rep78’ ... ]]
                                     
s0 . . . . . 0.379965 26.1441 . . . .
*/
mat G = 0.379965, 26.1441
lasso2 $model, lambda(1000) lglmnet nocons
assert mreldif(e(b),G) <1e-5

// single lambda, lasso, no standardisation
/*
> r<-glmnet(X,price,alpha=1,lambda=1000,standardize=FALSE,thresh=1e-15)
> t(coef(r))
1 x 11 sparse Matrix of class "dgCMatrix"
   [[ suppressing 11 column names ‘(Intercept)’, ‘mpg’, ‘rep78’ ... ]]
                                                
s0 10119.16 . . . . 3.549222 -62.23642 -106.0307
               
s0 6.079465 . .
*/
mat G = 3.549222, -62.23642, -106.0307, 6.079465, 10119.16
lasso2 $model, lambda(1000) lglmnet unitloadings
assert mreldif(e(b),G) <1e-5

// single lambda, lasso, no standardisation, noconstant
/*
> r<-glmnet(X,price,alpha=1,lambda=1000,standardize=FALSE,intercept=FALSE,thresh=1e-15)
> t(coef(r))
1 x 11 sparse Matrix of class "dgCMatrix"
   [[ suppressing 11 column names ‘(Intercept)’, ‘mpg’, ‘rep78’ ... ]]
                                                    
s0 . 6.830147 . . . 1.85882 -2.334519 . 3.886964 . .
*/
mat G = 6.830147, 1.85882, -2.334519, 3.886964
lasso2 $model, lambda(1000) lglmnet unitloadings nocons
// note looser tolerance
assert mreldif(e(b),G) <1e-4

// single lambda, lasso, mpg and foreign unpenalized
/*
> p.fac = rep(1, 10)
> p.fac[c(1, 10)] = 0
> r<-glmnet(X,price,alpha=1,lambda=500,penalty.factor=p.fac, thresh=1e-15)
> t(coef(r))
1 x 11 sparse Matrix of class "dgCMatrix"
   [[ suppressing 11 column names ‘(Intercept)’, ‘mpg’, ‘rep78’ ... ]]
                                                     
s0 9140.732 -225.1226 . . . . . . 6.065554 . 1962.096
*/
mat G = -225.1226, 6.065554, 1962.096, 9140.732
lasso2 $model, lambda(500) lglmnet notpen(mpg foreign)
assert mreldif(e(b),G) <1e-5

// single lambda, ridge, mpg and foreign unpenalized
/*
> p.fac = rep(1, 10)
> p.fac[c(1, 10)] = 0
> r<-glmnet(X,price,alpha=0,lambda=500,penalty.factor=p.fac, thresh=1e-15)
> t(coef(r))
1 x 11 sparse Matrix of class "dgCMatrix"
   [[ suppressing 11 column names ‘(Intercept)’, ‘mpg’, ‘rep78’ ... ]]
                                                        
s0 4420.934 -65.4987 167.7168 -446.8571 22.46099 1.26494
                                               
s0 2.118738 -25.1642 10.62293 -910.106 3185.862
*/
mat G = -65.4987, 167.7168, -446.8571, 22.46099, 1.26494, 			///
		2.118738, -25.1642, 10.62293, -910.106, 3185.862, 4420.934
lasso2 $model, lambda(500) alpha(0) lglmnet notpen(mpg foreign)
assert mreldif(e(b),G) <1e-5

// single lambda, elastic net, mpg and foreign unpenalized
/*
> p.fac = rep(1, 10)
> p.fac[c(1, 10)] = 0
> r<-glmnet(X,price,alpha=0.5,lambda=500,penalty.factor=p.fac, thresh=1e-15)
> t(coef(r))
1 x 11 sparse Matrix of class "dgCMatrix"
   [[ suppressing 11 column names ‘(Intercept)’, ‘mpg’, ‘rep78’ ... ]]
                                                             
s0 4089.962 -132.7721 . . . 0.7644731 . . 8.896648 . 2639.545
*/
mat G = -132.7721, 0.7644731, 8.896648, 2639.545, 4089.962
lasso2 $model, lambda(500) alpha(0.5) lglmnet notpen(mpg foreign)
assert mreldif(e(b),G) <1e-5

// single lambda, lasso, misc penalty loadings
/*
> p.fac = rep(1, 10)
> p.fac[c(1, 10)] = 4
> p.fac[c(2,9)] = 3
> p.fac[c(3,8)] = 2
> r<-glmnet(X,price,alpha=1,lambda=400,penalty.factor=p.fac, thresh=1e-15)
> t(coef(r))
1 x 11 sparse Matrix of class "dgCMatrix"
   [[ suppressing 11 column names ‘(Intercept)’, ‘mpg’, ‘rep78’ ... ]]
                                                     
s0 1254.559 . . . . 2.136852 . -42.89822 . . 393.3358
*/
mat G =  2.136852, -42.89822, 393.3358, 1254.559
mat psi = 4, 3, 2, 1, 1, 1, 1, 2, 3, 4
lasso2 $model, lambda(400) lglmnet ploadings(psi)
assert mreldif(e(b),G) <1e-5

// Elastic net, single lambda, misc penalty loadings
/*
> p.fac = rep(1, 10)
> p.fac[c(1, 10)] = 4
> p.fac[c(2,9)] = 3
> p.fac[c(3,8)] = 2
> r<-glmnet(X,price,alpha=0.5,lambda=400,penalty.factor=p.fac,thresh=1e-15)
> t(coef(r))
1 x 11 sparse Matrix of class "dgCMatrix"
   [[ suppressing 11 column names ‘(Intercept)’, ‘mpg’, ‘rep78’ ... ]]
                                                                
s0 1232.917 . 91.49722 -256.4055 . 2.29446 . -78.49864 5.61806 .
           
s0 1398.059
*/
mat G = 91.49722, -256.4055, 2.29446, -78.49864, 5.61806, 1398.059, 1232.917
mat psi = 4, 3, 2, 1, 1, 1, 1, 2, 3, 4
lasso2 $model, lambda(400) alpha(0.5) lglmnet ploadings(psi)
assert mreldif(e(b),G) <1e-5

// Elastic net, single lambda, misc penalty loadings, no standardization
/*
> p.fac = rep(1, 10)
> p.fac[c(1, 10)] = 4
> p.fac[c(2,9)] = 3
> p.fac[c(3,8)] = 2
> r<-glmnet(X,price,alpha=0.5,lambda=400,penalty.factor=p.fac,standardize=FALSE,thresh=1e-15)
> t(coef(r))
1 x 11 sparse Matrix of class "dgCMatrix"
   [[ suppressing 11 column names ‘(Intercept)’, ‘mpg’, ‘rep78’ ... ]]
                                                                   
s0 16039.13 -58.51055 322.0282 -218.2771 21.3079 4.093964 -72.31631
                         
s0 -247.7182 8.463531 . .
*/
mat G = -58.51055, 322.0282, -218.2771, 21.3079, 4.093964, -72.31631,	///
		-247.7182, 8.463531, 16039.13 
mat psi = 4, 3, 2, 1, 1, 1, 1, 2, 3, 4
lasso2 $model, lambda(400) alpha(0.5) lglmnet ploadings(psi) unitloadings
// note looser tolerance
assert mreldif(e(b),G) <1e-4

********************************************************************************
*** lglmnet option - consistency  with lasso2 default parameterization       ***
********************************************************************************

* uses auto dataset
* can change sd(y) and lambdas to suit
* note lglmnet implies prestd

// any lambda, any alpha
sysuse auto, clear
gen double y = price/10
local glmnetlambda_a = 100
local glmnetlambda_b = 10
drop if rep78==.
qui sum y
local sd = r(sd) * 1/sqrt( r(N)/(r(N)-1)   )
local glmnetalpha = 0.5
local L1lambda_a = `glmnetlambda_a' * 2 * 69
local L2lambda_a = `L1lambda_a' / `sd'
local L1lambda_b = `glmnetlambda_b' * 2 * 69
local L2lambda_b = `L1lambda_b' / `sd'
di "glmnetlambda_a=`glmnetlambda_a' L1lambda_a=`L1lambda_a' L2lambda_a=`L2lambda_a'"
di "glmnetlambda_b=`glmnetlambda_b' L1lambda_b=`L1lambda_b' L2lambda_b=`L2lambda_b'"
// Lasso / alpha=1
lasso2 y mpg-foreign, lambda(`glmnetlambda_a') lglmnet
local pmse = e(pmse)
storedresults save glmnet e()
lasso2 y mpg-foreign, lambda(`L1lambda_a') prestd
storedresults compare glmnet e(), tol(1e-8)							///
	exclude(macros: lasso2opt scalar: niter lambda pmse)
assert reldif(2*`pmse',e(pmse))<1e-8
// lambda list
lasso2 y mpg-foreign, lambda(`glmnetlambda_a' `glmnetlambda_b') lglmnet
storedresults save glmnet e()
lasso2 y mpg-foreign, lambda(`L1lambda_a' `L1lambda_b') prestd
storedresults compare glmnet e(), tol(1e-8)							///
	exclude(macros: lasso2opt										///
	scalar: niter lmax lmax0 lmin lmin0 laic laicc lbic lebic		///
	matrix: lambdamat lambdamat0)
// Enet
local alpha=(`glmnetalpha'*`sd')/( (1-`glmnetalpha') + `glmnetalpha'*`sd' )
local lambda_a=(1-`glmnetalpha')*`L2lambda_a' + `glmnetalpha'*`L1lambda_a'
local lambda_b=(1-`glmnetalpha')*`L2lambda_b' + `glmnetalpha'*`L1lambda_b'
di "alpha=`alpha' lambda_a=`lambda_a' lambda_b=`lambda_b'"
lasso2 y mpg-foreign, alpha(`glmnetalpha') lambda(`glmnetlambda_a') lglmnet
local pmse = e(pmse)
storedresults save glmnet e()
lasso2 y mpg-foreign, alpha(`alpha') lambda(`lambda_a') prestd
storedresults compare glmnet e(), tol(1e-8)							///
	exclude(macros: lasso2opt scalar: niter alpha lambda pmse)
di 2*`pmse'
di e(pmse)
assert reldif(2*`pmse',e(pmse))<1e-8
// lambda list
lasso2 y mpg-foreign, alpha(`glmnetalpha') lambda(`glmnetlambda_a' `glmnetlambda_b') lglmnet
storedresults save glmnet e()
lasso2 y mpg-foreign, alpha(`alpha') lambda(`lambda_a' `lambda_b') prestd
storedresults compare glmnet e(), tol(1e-8)							///
	exclude(macros: lasso2opt										///
	scalar: niter lmax lmax0 lmin lmin0 laic laicc lbic lebic alpha	///
	matrix: lambdamat lambdamat0)
// Ridge / alpha=0
lasso2 y mpg-foreign, alpha(0) lambda(`glmnetlambda_a') lglmnet
local pmse = e(pmse)
storedresults save glmnet e()
lasso2 y mpg-foreign, alpha(0) lambda(`L2lambda_a') prestd
storedresults compare glmnet e(), tol(1e-8)							///
	exclude(macros: lasso2opt scalar: niter lambda pmse)
di 2*`pmse'
di e(pmse)
assert reldif(2*`pmse',e(pmse))<1e-8
// lambda list
lasso2 y mpg-foreign, alpha(0) lambda(`glmnetlambda_a' `glmnetlambda_b') lglmnet
storedresults save glmnet e()
lasso2 y mpg-foreign, alpha(0) lambda(`L2lambda_a' `L2lambda_b') prestd
storedresults compare glmnet e(), tol(1e-8)							///
	exclude(macros: lasso2opt										///
	scalar: niter lmax lmax0 lmin lmin0 laic laicc lbic lebic		///
	matrix: lambdamat lambdamat0)


********************************************************************************
*** Validate ploadings(.) option                                             ***
********************************************************************************

* Default lasso2 behavior is to standardize, either on the fly or in advance.
* Confirm ploadings(.) work correctly by providing SDs of X.
* Behavior should be identical to standardization on-the-fly.

sysuse auto, clear
drop if rep78==.

cap mat drop psi_sd
foreach var of varlist mpg-foreign {
	qui sum `var'
	mat psi_sd = nullmat(psi_sd) , r(sd) * sqrt( (r(N)-1)/r(N)   )
}

// lasso
lasso2 price mpg-foreign, lambda(1000) alpha(1)
storedresults save nopload e()
lasso2 price mpg-foreign, lambda(1000) alpha(1) pload(psi_sd)
storedresults compare nopload e(), tol(1e-8)						///
	exclude(macros: lasso2opt)

// ridge
lasso2 price mpg-foreign, lambda(1000) alpha(0)
storedresults save nopload e()
lasso2 price mpg-foreign, lambda(1000) alpha(0) pload(psi_sd)
storedresults compare nopload e(), tol(1e-8)						///
	exclude(macros: lasso2opt)

// elastic net
lasso2 price mpg-foreign, lambda(1000) alpha(0.5)
storedresults save nopload e()
lasso2 price mpg-foreign, lambda(1000) alpha(0.5) pload(psi_sd)
storedresults compare nopload e(), tol(1e-8)						///
	exclude(macros: lasso2opt)


********************************************************************************
*** Validate adaptive lasso option                                           ***
********************************************************************************

sysuse auto, clear
drop if rep78==.
cap mat drop sdvec
// standardized variables used to get OLS with std coefs
foreach var of varlist price mpg-foreign {
	cap drop `var'_sd
	qui sum `var', meanonly
	qui gen double `var'_sd = `var'-r(mean)
	qui sum `var'
	qui replace `var'_sd = `var'_sd * 1/r(sd) * sqrt( r(N)/(r(N)-1) )
	mat sdvec = nullmat(sdvec), r(sd)/sqrt( r(N)/(r(N)-1) )
}
mat ysd = sdvec[1,1]
mat xsd = sdvec[1,2..11]

// replicate adaptive lasso

// use standardized coefficients + prestd
qui reg price_sd mpg_sd-foreign_sd

mata: st_matrix("psi",1:/abs(st_matrix("e(b)")))
mat psi = psi[1,1..10]
mat psi2 = J(1,10,1)

// lasso
lasso2 price mpg-foreign, lambda(10000) adaptive prestd
mat b=e(b)
lasso2 price mpg-foreign, lambda(10000) ploadings(psi) unitloadings prestd
assert mreldif(b,e(b)) < 1e-7
// elastic net
lasso2 price mpg-foreign, lambda(10000) adaptive prestd alpha(0.5)
mat b=e(b)
lasso2 price mpg-foreign, lambda(10000) ploadings(psi) ploadings2(psi2) unitloadings prestd alpha(0.5)
assert mreldif(b,e(b)) < 1e-7

// use unstandardized coefficients
qui reg price mpg-foreign

mata: st_matrix("psi",1:/abs(st_matrix("e(b)")))
// note that psi needs to be rescaled by sd(y) in order for lambda to be the same
mat psi = psi[1,1..10] * ysd[1,1]
mat psi2 = xsd

// lasso
lasso2 price mpg-foreign, lambda(10000) adaptive
mat b=e(b)
lasso2 price mpg-foreign, lambda(10000) ploadings(psi)
assert mreldif(b,e(b)) < 1e-7
// elastic net
lasso2 price mpg-foreign, lambda(10000) adaptive alpha(0.5)
mat b=e(b)
lasso2 price mpg-foreign, lambda(10000) ploadings(psi) ploadings2(psi2) alpha(0.5)
assert mreldif(b,e(b)) < 1e-7


* lglmnet version
* lglmnet standardizes automatically unless overridden by unitloadings.

// replicate adaptive lasso with lglmnet and no standardization
// unitloadings overrides standardization

// use unstandardized coefficients
qui reg price mpg-foreign

mata: st_matrix("psi",1:/abs(st_matrix("e(b)")))
mat psi = psi[1,1..10]
mat psi2 = J(1,10,1)

// lasso
lasso2 price mpg-foreign, lambda(1000) adaptive lglmnet unitloadings
mat b=e(b)
lasso2 price mpg-foreign, lambda(1000) ploadings(psi) unitloadings lglmnet
assert mreldif(b,e(b)) < 1e-7
// elastic net
lasso2 price mpg-foreign, lambda(1000) adaptive lglmnet unitloadings alpha(0.5)
mat b=e(b)
lasso2 price mpg-foreign, lambda(1000) ploadings(psi) ploadings2(psi2) unitloadings lglmnet alpha(0.5)
assert mreldif(b,e(b)) < 1e-7

// replicate adaptive lasso with lglmnet and standardization
// lglmnet standardizes by default unless overridden by unitloadings

// use standardized coefficients
// note dep var doesn't have to be standardized
qui reg price mpg_sd-foreign_sd

mata: st_matrix("psi",1:/abs(st_matrix("e(b)")))
mat psi = psi[1,1..10]
mat psi2 = J(1,10,1)

// lasso
lasso2 price mpg-foreign, lambda(1000) adaptive lglmnet
mat b=e(b)
lasso2 price mpg-foreign, lambda(1000) ploadings(psi) lglmnet
assert mreldif(b,e(b)) < 1e-7
// elastic net
lasso2 price mpg-foreign, lambda(1000) adaptive lglmnet alpha(0.5)
mat b=e(b)
lasso2 price mpg-foreign, lambda(1000) ploadings(psi) ploadings2(psi2) lglmnet alpha(0.5)
assert mreldif(b,e(b)) < 1e-7


********************************************************************************
*** replicate sqrt-lasso Matlab program										 ***
********************************************************************************

* load example data
insheet using "$prostate", tab clear
global model lpsa lcavol lweight age lbph svi lcp gleason pgg45

// uses the Matlab code "SqrtLassoIterative.m" (available on request)

lasso2 $model, sqrt l(40) unitload
mat a=e(betaAll)
/*
ans =

    0.3627
         0
         0
         0
         0
         0
         0
    0.0103
    1.7383
*/
mat b = (0.3627,0,0,0,0,0,0,0.0103,1.7383)
comparemat a b

lasso2 $model, sqrt l(10) unitload
mat a=e(betaAll)
/*
ans =

    0.5771
    0.1965
   -0.0092
    0.0773
    0.0685
         0
         0
    0.0063
    1.3946
*/
mat b = (0.5771,0.1965,-0.0092,0.0773,0.0685,0,0,0.0063,1.3946)
comparemat a b

lasso2 $model, sqrt l(1) unitload
mat a=e(betaAll)
/*
ans =

    0.5610
    0.5774
   -0.0196
    0.0950
    0.6700
   -0.0766
    0.0088
    0.0049
    0.5259
*/
mat b = (0.5610, 0.5774,-0.0196,0.0950,0.6700,-0.0766,0.0088,0.0049,0.5259)
comparemat a b

********************************************************************************
*** validation using Stata's elasticnet										 ***
********************************************************************************

// ridge
cap noi elasticnet linear $model, alphas(0) grid(2, min(0.25))
lassoselect alpha = 0 lambda = 0.25
lassocoef, display(coef, penalized)
mat b=e(b)
// Stata lambda = 2N*lambda
global L=2*e(N)*0.25
lasso2 $model, lambda($L) alpha(0)
assert mreldif(b,e(b))<1e-7

// lasso
cap noi elasticnet linear $model, alphas(1) grid(2, min(0.25))
lassoselect alpha = 1 lambda = 0.25
lassocoef, display(coef, penalized)
mat b=e(b)
// Stata lambda = 2N*lambda
global L=2*e(N)*0.25
lasso2 $model, lambda($L)
assert mreldif(b,e(b))<1e-7

// elastic net
cap noi elasticnet linear $model, alphas(0.5) grid(2, min(0.25))
lassoselect alpha = 0.5 lambda = 0.25
lassocoef, display(coef, penalized)
mat b=e(b)
// Stata lambda = 2N*lambda
global L=2*e(N)*0.25
lasso2 $model, lambda($L) alpha(0.5)
assert mreldif(b,e(b))<1e-7


********************************************************************************
*** norecover option														 ***
********************************************************************************

// partial() with constant
lasso2 $model, partial(age) l(50 20 10) 
mat A = e(betas)
mat A = A[2,1..9]
lasso2 $model, l(20) partial(age) postall
mat B = e(b)
comparemat A B

lasso2 $model, partial(age) l(50 20 10) nor 
mat A = e(betas)
mat A = A[2,1..7]
lasso2 $model, l(20) partial(age) nor postall
mat B = e(b)
comparemat A B

// partial() with constant, unitloadings
lasso2 $model, partial(age) l(50 20 10) unitl
mat A = e(betas)
mat A = A[2,1..9]
lasso2 $model, l(20) partial(age) postall unitl
mat B = e(b)
comparemat A B

lasso2 $model, partial(age) l(50 20 10) nor  unitl
mat A = e(betas)
mat A = A[2,1..7]
lasso2 $model, l(20) partial(age) nor postall unitl
mat B = e(b)
comparemat A B

// no partial() w/ constant, unitloadings
lasso2 $model, l(50 20 10) unitl
mat A = e(betas)
mat A = A[2,1..9]
lasso2 $model, l(20) postall unitl
mat B = e(b)
comparemat A B

lasso2 $model, l(50 20 10) nor  unitl
mat A = e(betas)
mat A = A[2,1..9]
lasso2 $model, l(20) nor postall unitl
mat B = e(b)
comparemat A B

// no partial() w/o constant, unit loadings
lasso2 $model, l(50 20 10) unitl nocons
mat A = e(betas)
mat A = A[2,1..8]
lasso2 $model, l(20) postall unitl nocons
mat B = e(b)
comparemat A B

lasso2 $model, l(50 20 10) nor unitl nocons
mat A = e(betas)
mat A = A[2,1..8]
lasso2 $model, l(20) nor postall unitl nocons
mat B = e(b)
comparemat A B


********************************************************************************
*** options																	 ***
********************************************************************************

cap lasso2 $model, alpha(0) sqrt
if _rc != 198 {
	exit 1
} 
*
// should say that lcount/lmax/lminr are being ignored
lasso2 $model, lambda(10) lcount(10)
lasso2 $model, lambda(10) lmax(100)
lasso2 $model, lambda(10) lminr(0.01)

// plotting only supported for lambda list
lasso2 $model, lambda(10) plotpath(lambda)

// incompatible options wrt penalty loadings
cap lasso2 $model, ploadings(abc) adaptive
if _rc != 198 {
	exit 1
} 
*
cap lasso2 $model, ploadings(abc) adatheta(3)
if _rc != 198 {
	exit 1
}
*

// var may not appear in partial() and notpen()
cap lasso2 $model, partial(age svi lcp) notpen(age svi)
if _rc != 198 {
	exit 1
} 
*

// controls the output and content of e(b)
lasso2 $model, l(20) displayall
lasso2 $model, l(20) postall
mat list e(b)

lasso2 $model, l(20) displayall postall
mat list e(b)

********************************************************************************
*** verify results are the same for scalar lambda vs lambda list			 ***
********************************************************************************

global lambdalist 150 130 100 80 60 30 10 5 3 1

* lasso
lasso2 $model, l($lambdalist)
mat A = e(betas)
local j=1
foreach i of numlist $lambdalist {
	mat a = A[`j',1..9]
	lasso2 $model, l(`i')
	mat b = e(betaAll)
	comparemat a b
	local j=`j'+1
}
*

* lasso (w/o constant)
lasso2 $model, l($lambdalist) nocons
mat A = e(betas)
local j=1
foreach i of numlist $lambdalist {
	mat a = A[`j',1..8]
	lasso2 $model, l(`i') nocons
	mat b = e(betaAll)
	comparemat a b
	local j=`j'+1
}
*

* post-lasso
lasso2 $model, l($lambdalist) ols
mat A = e(betas)
local j=1
foreach i of numlist $lambdalist {
	mat a = A[`j',1..9]
	lasso2 $model, l(`i')  
	mat b = e(betaAllOLS)
	comparemat a b
	local j=`j'+1
}
*

// replaced lambda=40 with lambda=60
// lambda=40 in path used beta=zeros from lambda=100
// as initial beta; would end up at poor solution
global sqrtlambdalist 100 60 20 10 5 1
* sqrt-lasso
lasso2 $model, l($sqrtlambdalist) sqrt   
mat A = e(betas)
local j=1
foreach i of numlist $sqrtlambdalist {
	mat a = A[`j',1..8]
	di `i'
	lasso2 $model, l(`i') sqrt 
	mat b = e(betaAll)
	mat b = b[1,1..8]
	comparemat a b
	local j=`j'+1
}
*

* post-sqrt-lasso ols
lasso2 $model, l($sqrtlambdalist) sqrt ols
mat A = e(betas)
local j=1
foreach i of numlist $sqrtlambdalist {
	di "this lambda: `i'"
	mat a = A[`j',1..9]
	lasso2 $model, l(`i') sqrt ols
	mat b = e(betaAllOLS)
	comparemat a b
	local j=`j'+1
}
*

* ridge
lasso2 $model, l($lambdalist) alpha(0)
mat A = e(betas)
local j=1
foreach i of numlist $lambdalist {
	mat a = A[`j',1..9]
	lasso2 $model, l(`i') alpha(0)
	mat b = e(betaAll)
	comparemat a b
	local j=`j'+1
}
*

* ols ridge
lasso2 $model, l($lambdalist) alpha(0) ols
mat A = e(betas)
local j=1
foreach i of numlist $lambdalist {
	mat a = A[`j',1..9]
	lasso2 $model, l(`i') alpha(0) ols
	mat b = e(betaAllOLS)
	comparemat a b
	local j=`j'+1
}
*
	
* elastic net
foreach ai of numlist 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 {
lasso2 $model, l($lambdalist) alpha(`ai')
mat A = e(betas)
local j=1
foreach i of numlist $lambdalist {
	mat a = A[`j',1..9]
	lasso2 $model, l(`i') alpha(`ai')
	mat b = e(betaAll)
	comparemat a b
	local j=`j'+1
}
}
*

* elastic net with ols
foreach ai of numlist 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 {
lasso2 $model, l($lambdalist) alpha(`ai') ols
mat A = e(betas)
local j=1
foreach i of numlist $lambdalist {
	mat a = A[`j',1..9]
	lasso2 $model, l(`i') alpha(`ai') ols
	mat b = e(betaAllOLS)
	comparemat a b
	local j=`j'+1
}
}
*


********************************************************************************
*** verify adapative weights												 ***
********************************************************************************

/*
** prestd should match no prestd (no prestd = ada weights incorporate scaling)
** ada theta = 1
lasso2 $model, lambda(10) adaptive
mat A = e(b)
lasso2 $model, lambda(10) adaptive prestd
mat B = e(b)
comparemat A B
** ada theta = 2
lasso2 $model, lambda(1) adaptive adatheta(2)
mat A = e(b)
lasso2 $model, lambda(1) adaptive adatheta(2) prestd
mat B = e(b)
comparemat A B


** lasso with ada theta = 1
lasso2 $model, adaptive verb
mat psi = e(Psi)

reg $model
mat bols = e(b)
// need to remove depvar from scaling of weights
qui sum `e(depvar)'
mat bols = bols * 1/r(sd) * sqrt(r(N)/(r(N)-1))

mat checkpsi = J(1,8,.)
forvalues i=1(1)8 {
	mat checkpsi[1,`i'] = abs(1/bols[1,`i'])
}
comparemat psi checkpsi


** lasso with ada theta = 2
lasso2 $model, adaptive verb adatheta(2)
mat psi = e(Psi)

reg $model
mat bols = e(b)
// need to remove depvar from scaling of weights
qui sum `e(depvar)'
mat bols = bols * 1/r(sd) * sqrt(r(N)/(r(N)-1))
global depvar `e(depvar)'
global indepvars	: list global(model) - global(depvar)
mat checkpsi = J(1,8,.)
// need to remove theta from scaling of weights
forvalues i=1(1)8 {
	mat checkpsi[1,`i'] = abs(1/bols[1,`i'])^2
	local v		: word `i' of $indepvars
	qui sum `v'
	mat checkpsi[1,`i'] = checkpsi[1,`i'] * (1/r(sd) * sqrt(r(N)/(r(N)-1)))
}
comparemat psi checkpsi

// use of adaloadings option

** prestd should match no prestd (no prestd = ada weights incorporate scaling)
** ada theta = 1
lasso2 $model , l(10) alph(0)
mat b = e(betaAll)
lasso2 $model, lambda(10) adaloadings(b) adaptive
mat A = e(b)
lasso2 $model, lambda(10) adaloadings(b) adaptive prestd
mat B = e(b)
comparemat A B
** ada theta = 2
lasso2 $model , l(10) alph(0)
mat b = e(betaAll)
lasso2 $model, lambda(1) adaloadings(b) adaptive adatheta(2)
mat A = e(b)
lasso2 $model, lambda(1) adaloadings(b) adaptive adatheta(2) prestd
mat B = e(b)
comparemat A B

// adaloadings with theta=1
lasso2 $model , l(10) alph(0)
mat b = e(betaAll)
lasso2 $model, adaptive adal(b) adat(1)
mat psi = e(Psi)
mat checkpsi = J(1,8,.)
forvalues i=1(1)8 {
	mat checkpsi[1,`i'] = abs(1/b[1,`i'])
}
comparemat psi checkpsi

// adaloadings with theta=2
lasso2 $model , l(10) alph(0)
mat b = e(betaAll)
lasso2 $model, adaptive adal(b) adat(2)
mat psi = e(Psi)
global depvar `e(depvar)'
global indepvars	: list global(model) - global(depvar)
mat checkpsi = J(1,8,.)
// need to remove theta from scaling of weights
forvalues i=1(1)8 {
	mat checkpsi[1,`i'] = abs(1/b[1,`i'])^2
	local v		: word `i' of $indepvars
	qui sum `v'
	mat checkpsi[1,`i'] = checkpsi[1,`i'] * (1/r(sd) * sqrt(r(N)/(r(N)-1)))
}
comparemat psi checkpsi
*/

********************************************************************************
*** pre-estimation standardisation vs std on the fly   				 		 ***
********************************************************************************

// lasso
// standardisation using penalty loadings (default)
lasso2 $model, l(10) 
mat A = e(beta)
// pre-estimation standardisation of data
lasso2 $model, l(10) prestd  
mat B = e(beta)
comparemat A B , tol(10e-6)

// lasso [nocons]
// standardisation using penalty loadings (default)
lasso2 $model, l(10)  nocons
mat A = e(beta)
// pre-estimation standardisation of data
lasso2 $model, l(10) prestd nocons 
mat B = e(beta)
comparemat A B , tol(10e-6)

// sqrt lasso
// standardisation using penalty loadings (default)
lasso2 $model, l(10) sqrt 
mat A = e(beta)
// pre-estimation standardisation of data
lasso2 $model, l(10) sqrt prestd
mat B = e(beta)
comparemat A B , tol(10e-6)

// sqrt lasso [nocons]
// standardisation using penalty loadings (default)
lasso2 $model, l(10) sqrt nocons
mat A = e(beta)
// pre-estimation standardisation of data
lasso2 $model, l(10) sqrt prestd nocons
mat B = e(beta)
comparemat A B , tol(10e-6)

// elastic net
foreach lam of numlist 1 10 50 150 160 {
foreach ai of numlist 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 {
	// in original units
	// standardisation using penalty loadings (default)
	lasso2 $model, l(`lam')  alpha(`ai')
	mat A = e(beta)

	// pre-estimation standardisation of data
	lasso2 $model, l(`lam') prestd alpha(`ai')
	mat B = e(beta)
	comparemat A B , tol(10e-6)
}
}
*

// elastic net [nocons]
foreach lam of numlist 1 10 50 150 160 {
foreach ai of numlist 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 {
	// in original units
	// standardisation using penalty loadings (default)
	lasso2 $model, l(`lam')  alpha(`ai') nocons
	mat A = e(beta)

	// pre-estimation standardisation of data
	lasso2 $model, l(`lam') prestd alpha(`ai') nocons
	mat B = e(beta)
	comparemat A B , tol(10e-6)
}
}
*

********************************************************************************
*** verify ridge regression results											 ***
********************************************************************************

lasso2 $model, l(150) alpha(0) unitload
mat A = e(beta)
mat A = A[1,1..8] // excl intercept

estridge $model, l(150)
mat B = r(bhat)

comparemat A B , tol(10e-6)


********************************************************************************
*** verify post-estimation OLS results										 ***
********************************************************************************

* lasso
foreach i of numlist 0.5 1 4 10 15 50 150 {
	lasso2 $model, l(`i') ols
	mat A = e(b)
	reg lpsa `e(selected)'
	mat B = e(b)
	comparemat A B
}
*

* lasso (w/o constant)
foreach i of numlist 0.5 1 4 10 15 50 150 {
	lasso2 $model, l(`i') ols nocons
	mat A = e(b)
	reg lpsa `e(selected)', nocons
	mat B = e(b)
	comparemat A B
}
*

* sqrt lasso
foreach i of numlist 0.5 1 4 10 15 50 150 {
	lasso2 $model, l(`i') ols sqrt
	mat A = e(b)
	reg lpsa `e(selected)'
	mat B = e(b)
	comparemat A B
}
*

* elastic net
foreach i of numlist 0.5 1 4 10 15 50 150 {
	lasso2 $model, l(`i') ols alpha(.5)
	mat A = e(b)
	reg lpsa `e(selected)'
	mat B = e(b)
	comparemat A B
}
*

* ridge
foreach i of numlist 0.5 1 4 10 15 50 150 {
	lasso2 $model, l(`i') ols alpha(0)
	mat A = e(b)
	reg lpsa `e(selected)'
	mat B = e(b)
	comparemat A B
}
*

********************************************************************************
*** partial() vs notpen()													 ***
********************************************************************************

lasso2 $model, partial(lcp) l(50)
mat A = e(b) 
lasso2 $model, notpen(lcp) l(50)
mat B = e(b)
comparevec A B  

lasso2 $model, partial(lcp) l(50) sqrt
mat A = e(b) 
lasso2 $model, notpen(lcp) l(50) sqrt
mat B = e(b)
comparevec A B  

lasso2 $model, partial(lcp) l(50) alpha(0.5)
mat A = e(b) 
lasso2 $model, notpen(lcp) l(50) alpha(0.5)
mat B = e(b)
comparevec A B  

lasso2 $model, partial(lcp) l(50) alpha(0)
mat A = e(b) 
lasso2 $model, notpen(lcp) l(50) alpha(0)
mat B = e(b)
comparevec A B  

lasso2 $model, lambda(10) partial(age) notpen(svi)
mat A = e(b) 
lasso2 $model, lambda(10) partial(svi) notpen(age)
mat B = e(b)
comparevec A B  

********************************************************************************
*** penalty loadings vs notpen (see help file)						         ***
********************************************************************************

lasso2 $model, l(10) notpen(lcavol) unitloadings
mat A = e(b)

mat myloadings = (0,1,1,1,1,1,1,1)
lasso2 $model, l(10) ploadings(myloadings)
mat B = e(b)

comparemat A B


********************************************************************************
*** ic option to control display of output 							***
********************************************************************************

di as red "should display EBIC (the default):"
lasso2 $model  
sleep 1000
di as red "should display AIC:"
lasso2 $model , ic(aic)
sleep 1000
di as red "should display AICc:"
lasso2 $model , ic(aicc)
sleep 1000
di as red "should display BIC:"
lasso2 $model , ic(bic)
sleep 1000
di as red "should display EBIC:"
lasso2 $model , ic(ebic)


********************************************************************************
*** degrees of freedom calculation          						   ***
********************************************************************************

// replicate dof w/o constant and no standardisation [OK]
lasso2 $model ,  alpha(0) l(20 .1)  long unitload    nocons nopath
mat D= e(dof)
mat list e(dof)

putmata y=(lpsa) X=(lcavol lweight  age lbph  svi  lcp gleason pgg45  ), replace
mata: df=trace(X*invsym((X'X):+20/2*I(8))*X') 
mata: df

mata: st_local("df",strofreal(df))
assert reldif(el(D,1,1),`df')<10^-6

// standardisation w/o constant [OK]
lasso2 $model ,  alpha(0)   l(20 .1)  long      nocons nopath
mat D1 = e(dof)
mat list e(dof)
lasso2 $model ,  alpha(0)   l(20 .1)  long prestd   nocons nopath
mat D2 = e(dof)
mat list e(dof)
comparemat D1 D2

putmata y=(lpsa) X=(lcavol lweight  age lbph  svi  lcp gleason pgg45  ), replace
mata: s = sqrt(mean((X:-mean(X)):^2))
mata: ssq = s:^2
mata: Xs=X:/s 
mata: df1=trace(X*invsym((X'X):+20/2*diag(ssq))*X')  // "on the fly" standardisation
mata: df2=trace(Xs*invsym((Xs'Xs):+20/2*I(8))*Xs')   // pre-standardisation
mata: df1,df2
mata: st_local("df1",strofreal(df1))
mata: st_local("df2",strofreal(df2))
assert reldif(el(D1,1,1),`df1')<10^-6
assert reldif(el(D1,1,1),`df2')<10^-6
assert reldif(`df1',`df2')<10^-6

// dof with constant and no standardisation [OK]
lasso2 $model ,  alpha(0) l(20 .1)  long unitload  nopath   
mat list e(dof)
mat D = e(dof)

putmata y=(lpsa) X=(lcavol lweight  age lbph  svi  lcp gleason pgg45  ), replace
mata: Xone=(X,J(97,1,1))
mata: Psicons = I(9)
mata: Psicons[9,9]=0
mata: Xdm = X :- mean(X)
mata: trace(X*invsym((X'X):+20/2*I(8))*X') 			// this is wrong (ignores constant)
mata: trace(Xdm*invsym((Xdm'Xdm):+20/2*I(8))*Xdm')   // this is correct but missing the constant
mata: trace(Xone*invsym((Xone'Xone):+20/2*Psicons)*Xone')  // this is correct (incl constant)
mata: df1=trace(Xone*invsym((Xone'Xone):+20/2*Psicons)*Xone')  // this should be correct
mata: df2=trace(Xdm*invsym((Xdm'Xdm):+20/2*I(8))*Xdm')+1  // this is correct after +1 for constant
mata: df1,df2
mata: st_local("df1",strofreal(df1))
mata: st_local("df2",strofreal(df2))
assert reldif(el(D,1,1),`df1')<10^-6
assert reldif(el(D,1,1),`df2')<10^-6
assert reldif(`df1',`df2')<10^-6

// standardisation w/  constant [OK]
lasso2 $model ,  alpha(0) l(20 10 1 .1)  long nopath
mat list e(dof)
mat D1 = e(dof)
lasso2 $model ,  alpha(0) l(20 10 1 .1)  long prestd nopath  
mat list e(dof)
mat D2 = e(dof)
comparemat D1 D2

putmata y=(lpsa) X=(lcavol lweight  age lbph  svi  lcp gleason pgg45), replace
mata: s = sqrt(mean((X:-mean(X)):^2))
mata: ssq = s:^2
mata: Xs=(X:-mean(X)):/s 
mata: Xone=(X,J(97,1,1))
mata: df1=trace(Xdm*invsym((Xdm'Xdm):+20/2*diag(ssq))*Xdm') +1
mata: df2=trace(Xs*invsym((Xs'Xs):+20/2*I(8))*Xs') +1
mata: df1,df2
mata: st_local("df1",strofreal(df1))
mata: st_local("df2",strofreal(df2))
assert reldif(el(D1,1,1),`df1')<10^-6
assert reldif(el(D1,1,1),`df2')<10^-6
assert reldif(`df1',`df2')<10^-6

********************************************************************************
*** lic option 						    							***
********************************************************************************

* check that right lambda is used 
foreach ic of newlist ebic aic aicc bic {
	lasso2 $model 
	local optlambda=e(l`ic') 
	lasso2, lic(`ic') postres
	local thislambda=e(lambda)
	assert reldif(`optlambda',`thislambda')<10^-8
}
*

********************************************************************************
*** predicted values (see help file)										 ***
********************************************************************************

// xbhat1 is generated by re-estimating the model for lambda=10.  The noisily 
// option triggers the display of the
// estimation results.  xbhat2 is generated by linear approximation using the 
// two beta estimates closest to
//    lambda=10.

* load example data
insheet using "$prostate", tab clear
global model lpsa lcavol lweight age lbph svi lcp gleason pgg45

lasso2 $model
cap drop xbhat1
predict double xbhat1, xb l(10) noisily
cap drop xbhat2
predict double xbhat2, xb l(10) approx

//    The model is estimated explicitly using lambda=100.  If lasso2 is 
//called with a scalar lambda value, the
//   subsequent predict command requires no lambda() option.
lasso2 $model, lambda(10)
cap drop xbhat3
predict double xbhat3, xb

//    All three methods yield the same results.  However note that the linear 
// approximation is only exact for the lasso
//   which is piecewise linear.
assert (xbhat1-xbhat2<10e-8) & (xbhat3-xbhat2<10e-8) 

//It is also possible to obtain predicted values by referencing a specific
// lambda ID using the lid() option.
lasso2 $model
cap drop xbhat4
predict double xbhat4, xb lid(21)
cap drop xbhat5
predict double xbhat5, xb l(25.45473900468241)
assert (xbhat4-xbhat5<10e-8)


********************************************************************************
*** misc options/syntax checks		                                         ***
********************************************************************************

// Support for inrange(.) and similar [if] expressions:
lasso2 $model if inrange(age,50,70)


********************************************************************************
*** plotting                                                                 ***
********************************************************************************

lasso2 $model

lasso2, plotpath(lambda) plotlabel plotopt(legend(off))

lasso2, plotpath(lnlambda) plotlabel plotopt(legend(off))

lasso2, plotpath(norm) plotlabel plotopt(legend(off))

lasso2, plotpath(norm) plotlabel plotopt(legend(off)) plotvar(lcavol)

********************************************************************************
*** validate RSS / r-squared    				                             ***
********************************************************************************

// loop over three mehods. "nopath" corresponds to default standardisation on the fly.
// "nopath" is just a placeholder that doesn't affect calculations.
foreach method in prestd unitl nopath {
foreach a of numlist 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 {
	lasso2 $model, alpha(`a') l(10)  `method'
	cap drop xb
	cap drop r
	predict double xb, xb
	predict double r, r

	sum lpsa
	di r(sd)^2*98
	local tss=r(sd)^2*98
	sum xb
	di r(sd)^2*98
	local ess=r(sd)^2*98
	sum r
	di r(sd)^2*98
	local rss=r(sd)^2*98
	 
	di "r-squared"
	local rsq = 1-`rss'/`tss'
	di `rsq'

	lasso2 $model, alpha(`a') l(10 0.1) `method'
	mat list e(rsq)
	mat RSQ = e(rsq)
	di el(RSQ,1,1)
	di reldif(`rsq',el(RSQ,1,1))
	assert reldif(`rsq',el(RSQ,1,1))<10^(-3)
}
}
*

********************************************************************************
*** validate EBIC default gamma     				                         ***
********************************************************************************

webuse air2, clear

lasso2 air L(1/24).air

local myebicgamma = 1-log(e(N))/(2*log(e(p)))

di `myebicgamma'
di e(ebicgamma)
assert reldif(`myebicgamma',e(ebicgamma))<10^(-3)

********************************************************************************
*** panel example: validate within transformation                            ***
********************************************************************************

use "http://fmwww.bc.edu/ec-p/data/macro/abdata.dta", clear

lasso2 ys l(0/3).k l(0/3).n, fe l(10)  
mat A = e(b)

lasso2 ys l(0/3).k l(0/3).n ibn.id, partial(ibn.id) l(10) nor
mat B = e(b)

comparemat A B

// noftools option
lasso2 ys l(0/3).k l(0/3).n, fe l(10)
savedresults save ftools e()
cap noi assert "`e(noftools)'"==""  // will be error if ftools not installed
lasso2 ys l(0/3).k l(0/3).n, fe l(10) noftools
assert "`e(noftools)'"=="noftools"
savedresults comp ftools e(), exclude(macros: lasso2opt)

********************************************************************************
***  check if partial() works with fe				                         ***
***  and followed by lic()													 ***
********************************************************************************

clear 
use https://www.stata-press.com/data/r16/nlswork

lasso2 ln_w grade age c.age#c.age ttl_exp c.ttl_exp#c.ttl_exp tenure ///
		c.tenure#c.tenure 2.race not_smsa south i.year, fe  
ereturn list 

lasso2 ln_w grade age c.age#c.age ttl_exp c.ttl_exp#c.ttl_exp tenure ///
		c.tenure#c.tenure 2.race not_smsa south i.year, fe partial(i.year) 
lasso2, lic(ebic)

lasso2 ln_w grade age c.age#c.age ttl_exp c.ttl_exp#c.ttl_exp tenure ///
		c.tenure#c.tenure 2.race not_smsa south i.year, fe partial(i.year) ///
		lic(ebic)

********************************************************************************
***  check residuals with fe												 ***
********************************************************************************
		
clear
use https://www.stata-press.com/data/r16/nlswork

foreach opt in xb u e ue xbu {
	cap drop `opt'hat_lasso2
	cap drop `opt'hat_xtreg
}
cap drop esample

**replace ln_w = . if year == 80

lasso2 ln_w grade age c.age#c.age ttl_exp c.ttl_exp#c.ttl_exp tenure ///
		c.tenure#c.tenure 2.race not_smsa south , fe
gen byte esample=e(sample)
lasso2, lic(ebic) postres ols
mat bl2 = e(b)
local selected = e(selected)
di "`selected'"
//ereturn list

// confirm post-lasso LSDV matches xtreg,fe
foreach opt in xb u e ue xbu {
	predict double `opt'hat_lasso2, ols `opt'
}
xtreg ln_w `selected' if esample, fe
foreach opt in xb u e ue xbu {
	predict double `opt'hat_xtreg, `opt'
}
foreach opt in xb u e ue xbu {
	di "checking option `opt'
	assert reldif(`opt'hat_lasso2, `opt'hat_xtreg) < 1e-6
}


** alt approach from 2nd half
clear
use https://www.stata-press.com/data/r16/nlswork

replace ln_w = . if year == 80

lasso2 ln_w grade age c.age#c.age ttl_exp c.ttl_exp#c.ttl_exp tenure ///
		c.tenure#c.tenure 2.race not_smsa south , fe   
lasso2, lic(ebic) postres ols	
		
local sel = e(selected)
di "`sel'"

predict double uehat  , ue noi  
predict double ehat  , e noi  
predict double xbhat  , xb noi  
predict double xbuhat  , xbu noi  
predict double uhat  , u noi  

xtreg ln_w `sel' if e(sample), fe 
mat bxtreg = e(b)

predict double uehat_xtreg   , ue  
predict double ehat_xtreg   , e 
predict double xbhat_xtreg  , xb
predict double xbuhat_xtreg  , xbu
predict double uhat_xtreg  , u

assert abs(ehat_xtreg-ehat)<10e-8 | (missing(ehat_xtreg) | missing(ehat))
assert abs(uehat_xtreg-uehat)<10e-8 | (missing(uehat_xtreg) | missing(uehat)) 
assert abs(xbhat_xtreg-xbhat)<10e-8 | (missing(xbhat_xtreg) | missing(xbhat))
assert abs(xbuhat_xtreg-xbuhat)<10e-8 | (missing(xbuhat_xtreg) | missing(xbuhat))
assert abs(uhat_xtreg-uhat)<10e-8 | (missing(uhat_xtreg) | missing(uhat))



********************************************************************************
*** finish                                                                   ***
********************************************************************************

cap log close
set rmsg off
