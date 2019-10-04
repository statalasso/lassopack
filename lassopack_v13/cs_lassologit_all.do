
clear all 

cap cd "C:\Users\Achim.Ahrens\Dropbox\StataLasso\lassoLogit"
cap cd "/home/achim/Dropbox/StataLasso/lassoLogit/"
cap cd "C:\Users\achim\Dropbox\StataLasso\lassoLogit"
adopath + "`c(pwd)'"


cap log close
log using "cs_lassologit_all.smcl", replace

which lassologit
which lassologit_p
which cvlassologit
which rlassologit

do "cs_lassologit.do"

do "cs_cvlassologit.do"

do "cs_rlassologit.do"

do "cs_predict.do"

do "cs_weights.do"

cap log close


