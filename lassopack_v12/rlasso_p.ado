*! rlasso_p 1.0.00 30/11/2017
*! lassopack package 1.2 15jan2019
*! authors aa/cbh/ms

* postestimation predict for rlasso

program define rlasso_p, rclass

	version 12.1

	syntax namelist(min=1 max=2) [if] [in], ///
											///
				[XB 						/// [default]
				Residuals 					///
				lasso						///
				ols							///
											///
				]

	* create variable here
	tokenize `namelist'
	if "`2'"=="" {					//  only new varname provided
		local varlist `1'
		qui gen `1' = .
	}
	else {							//  datatype also provided
		local vtype `1'
		local varlist `2'
		qui gen `1' `2' = .
	}
	*

	local command=e(cmd)
	if ("`command'"~="rlasso") {
		di as err "error: -rlasso_p- supports only the -rlasso- command"
		exit 198
	}
	*
	
	if "`lasso'"~="" & "`ols'"~="" {
		di as err "error: incompatible options -lasso- and -ols-"
		exit 198
	}
	
	marksample touse, novarlist
	
	*** warning messages
	if ("`xb'`residuals'"=="") {
		di as text "No xb or residuals options specified. Assume xb (fitted values)."
		local xb xb
	}
	*** fe currently not supported.
	if `e(fe)' {
		di as err "predict not currently supported after FE estimation"
		exit 198
	}
	*
	
	*** obtain prediction/residuals
	tempname b
	tempvar xbvar
	if "`lasso'`ols'"=="" {					//  default = posted e(b) matrix
		mat `b'		=e(b)
	}
	else if "`lasso'"~="" {
		mat `b'		=e(beta)
	}
	else {
		mat `b'		=e(betaOLS)
	}
	
	qui matrix score `typlist' `xbvar'= `b'  if `touse'
	if ("`xb'"!="") {
		qui replace `varlist' = `xbvar' if `touse'
		label var `varlist' "Predicted values"
	}
	else if ("`residuals'"!="") {
		qui replace `varlist' = `e(depvar)' - `xbvar' if `touse'
		label var `varlist' "Residuals"
	}
	*

end
