
clear all 

cap cd "/Users/kahrens/MyProjects/lassopack"
cap cd "C:\LocalStore\ecomes\Dropbox\StataLasso\lassoLogit"
adopath + "`c(pwd)'"

global loadspam  insheet using https://archive.ics.uci.edu/ml/machine-learning-databases/spambase/spambase.data, clear comma
global loadprostate insheet using https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data, clear tab

cap log close
log using "cs_lassologit_all.smcl", replace

which lassologit
which lassologit_p
which cvlassologit
which rlassologit

do "cs_lassologit.do" // ok

do "cs_cvlassologit.do" // 

do "cs_rlassologit.do"

do "cs_lassologit_predict.do"

do "cs_lassologit_weights.do"

cap log close


