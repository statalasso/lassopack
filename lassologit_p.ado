*! lassologit_p
*! part of lassopack v1.4.2
*! last edited: 21mar2021
*! authors: aa/ms

program define lassologit_p, rclass

	syntax namelist(min=1 max=2) [if] [in] , [	XB ///
												Pr ///
												Class ///
												lopt lse ///
												POSTLogit ///
												NOIsily ///
												lic(string) /// 
												est ]
	
	* create variable here
	tokenize `namelist'
	if "`2'"=="" {					//  only new varname provided
		local varlist `1'
	}
	else {							//  datatype also provided
		local vtype `1'
		local varlist `2'
	}
	*
	
	
	// show estimation
	if "`noisily'"=="" {
		local qui qui
	}
	
	// check if lambda is defined
	if ("`e(cmd)'"=="cvlassologit") {
		if ("`lopt'"!="") {
			`qui' cvlassologit, lopt postresults
		}
		else if ("`lse'"!="") {
			`qui' cvlassologit, lse postresults
		}
		else {
			di as err "predict requires 'lse' or 'lopt' after cvlassologit."
			exit 1
		}

		// use post-logit or just logistic lasso?
		tempname betaused
		if ("`postlogit'"!="") {
			di as text "Uses post-logit for prediction."
			mat `betaused' = e(beta_post_dense)	
		}
		else {
			mat `betaused' = e(b)	
		}
		//
	}
	else if ("`e(cmd)'"=="lassologit"&"`est'"!="") {
	
		// information criteria
		if ("`lic'"=="ebic") | ("`lic'"=="aic") | ("`lic'"=="aicc") | ("`lic'"=="bic") {
			lassologit, lic(`lic') postresults
		}
		else if ("`lic'"!="ebic") & ("`lic'"!="aic") & ("`lic'"!="aicc") & ("`lic'"=="bic") & ("`lic'"!="") {
			di as err "lic(`lic') not allowed."
			exit 1
		}
		
		// check if lambda is defined
		if (`e(lcount)'>1) {
			di as err "No single lambda specified."
			di as err "Use predict with lic() option, or lassologit with postresults." 
		}

		// use post-logit or just logistic lasso?
		tempname betaused
		if ("`postlogit'"!="") {
			di as text "Uses post-logit for prediction."
			mat `betaused' = e(betas)	
		}
		else {
			mat `betaused' = e(b)	
		}
		//
		
	}
	else if ("`e(cmd)'"=="lassologit"&"`est'"=="") {

		// get ID of optimal lambda
		if ("`lic'"=="ebic") {
			local optid = e(ebicid)
		}
		else if ("`lic'"=="aic") {
			local optid = e(aicid)
		}
		else if ("`lic'"=="aicc") {
			local optid = e(aiccid)
		}
		else if ("`lic'"=="bic") {
			local optid = e(bicid)
		}
		else if ("`lic'"!="ebic") & ("`lic'"!="aic") & ("`lic'"!="aicc") & ("`lic'"=="bic") & ("`lic'"!="") {
			di as err "lic(`lic') not allowed."
			exit 1
		}

		local postlogit_est = `e(postlogit)'
		local betascols = `e(p)'+`e(cons)'
		tempname betaused

		if "`postlogit'"!=""&`e(postlogit)'==0 {
			di as err "postlogit option ignored"
		}
		if `e(postlogit)'==0 di as text "using lasso coefficients stored in e(betas)"
		if `e(postlogit)'==1 di as text "using post-lasso coefficients stored in e(betas)"
		mat `betaused'=e(betas)
		mat `betaused'=`betaused'[`optid',1..`betascols']

	}
	else if ("`e(cmd)'"!="rlassologit") {
		
		di as err "command not recognized."
		exit 198
		
	}
	

	
	if ("`noisily'"!="") {
		di as text "Beta used for prediction:"
		mat list `betaused'
	}
	//
	
	*** obtain prediction
	tempvar xbvar
	qui matrix score `vtype'  `xbvar'= `betaused'  `if'

	if ("`xb'"!="") {
	
		qui gen `vtype' `varlist' = `xbvar' `if'
		label var `varlist' "Linear prediction"
		
	}
	else if ("`pr'"!="") {
	
		qui gen `vtype' `varlist' = exp(`xbvar')/(1+exp(`xbvar')) `if' 
		label var `varlist' "Predicted probabilities"
		
	}
	else if ("`class'"!="") {
	
		qui gen `vtype' `varlist' = exp(`xbvar')/(1+exp(`xbvar')) `if' 
		qui replace `varlist' = (`varlist'>0.5) `if' 
		label var `varlist' "Predicted class"
		
	}
	*
	
end
