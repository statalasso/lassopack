clear all 


********************************************************************************
***  verify fweights 														 ***
********************************************************************************

$loadspam

// frequencies
set seed 123
gen freq = max(rpoisson(3),1)

// estimate with fweights
lassologit v58 v1-v57 [fw=freq]
mat A = e(betas)
mat LA = e(lambdas)

// estimate with expanded data
expand freq
lassologit v58 v1-v57
mat B = e(betas)
mat LB = e(lambdas)

assert mreldif(LA,LB)<1e-10
assert mreldif(A,B)<1e-10

********************************************************************************
***  verify fweights (without constant)										 ***
********************************************************************************

$loadspam

// frequencies
set seed 123
gen freq = max(rpoisson(3),1)

// estimate with fweights
lassologit v58 v1-v57 [fw=freq], nocons
mat A = e(betas)
mat LA = e(lambdas)

// estimate with expanded data
expand freq
lassologit v58 v1-v57, nocons
mat B = e(betas)
mat LB = e(lambdas)

assert mreldif(LA,LB)<1e-10
assert mreldif(A,B)<1e-10

********************************************************************************
***  verify aweights														 ***
********************************************************************************

$loadspam

// aweights / pweights
// check by setting lambda very small
// => lasso logit and post-logit should be similar
set seed 123
gen double wtvar=runiform()
lassologit v58 v1-v10 [aw=wtvar], l(0.000001) tolopt(1e-15) tolzero(1e-15)
mat bL=e(beta_dense)
mat bPL=e(beta_post_dense)
assert mreldif(bL,bPL)<0.002

// confirm pw and aw are equivalent
lassologit v58 v1-v5 [aw=wtvar], l(0.05)
mat bL=e(beta_dense)
mat bPL=e(beta_post_dense)
lassologit v58 v1-v5 [pw=wtvar], l(0.05)
assert mreldif(bL,e(beta_dense))<1e-10
assert mreldif(bPL,e(beta_post_dense))<1e-10

