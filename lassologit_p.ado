*! lassologit_p
*! part of lassopack v1.4.4
*! last edited: 27mar2021
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
			`qui' cvlassologit, lopt postresults `postlogit'
		}
		else if ("`lse'"!="") {
			`qui' cvlassologit, lse postresults `postlogit'
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
			di as text "Uses lasso for prediction"
			mat `betaused' = e(beta_dense)	
		}
		//
	}
	else if ("`e(cmd)'"=="lassologit"&`e(lcount)'==1&"`est'"=="") {

		// use post-logit or just logistic lasso?
		tempname betaused
		if ("`postlogit'"!="") {
			di as text "Uses post-logit for prediction."
			mat `betaused' = e(beta_post_dense)	
		}
		else {
			di as text "Uses lasso for prediction"
			mat `betaused' = e(beta_dense)	
		}
		//

	}
	else if ("`e(cmd)'"=="lassologit"&"`est'"!="") {
		
		if ("`lic'"=="ebic") | ("`lic'"=="aic") | ("`lic'"=="aicc") | ("`lic'"=="bic") {
			di as text "re-estimating the model with lic(`lic')"
			lassologit, lic(`lic') postresults  
		}
		else if ("`lic'"!="ebic") & ("`lic'"!="aic") & ("`lic'"!="aicc") & ("`lic'"=="bic") & ("`lic'"!="") {
			di as err "lic(`lic') not allowed."
			exit 198
		}
		else if ("`lic'"=="") {
			di as err "No single lambda specified. Use lic() option." 
			exit 198
		}
		
		// use post-logit or just logistic lasso?
		tempname betaused
		if ("`postlogit'"!="") {
			di as text "Uses post-logit for prediction."
			mat `betaused' = e(beta_post_dense)	
		}
		else {
			di as text "Uses lasso for prediction"
			mat `betaused' = e(beta_dense)	
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
		else if ("`lic'"=="") {
			di as err "lic() option required"
			exit 198
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

	if ("`pr'`xb'`class'"=="") local pr pr

	if ("`pr'"!="") {
	
		di as text "storing predicted probabilities"
		qui gen `vtype' `varlist' = exp(`xbvar')/(1+exp(`xbvar')) `if' 
		label var `varlist' "Predicted probabilities"
		
	}
	else if ("`xb'"!="") {
	
		di as text "storing linear predictions"
		qui gen `vtype' `varlist' = `xbvar' `if'
		label var `varlist' "Linear prediction"
		
	}
	else if ("`class'"!="") {
	
		di as text "storing predicted class"
		qui gen `vtype' `varlist' = exp(`xbvar')/(1+exp(`xbvar')) `if' 
		qui replace `varlist' = (`varlist'>0.5) `if' 
		label var `varlist' "Predicted class"
		
	}
	*
	
end
