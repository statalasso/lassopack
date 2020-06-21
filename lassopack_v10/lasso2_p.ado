*! lasso2_p 1.0.02 05/04/2018
*! authors aa/ms
*
* post-estimation predict for both lasso2 and cvlasso.
*
* Updates (release date):
* 1.0.02  (5apr2018 - not released)
*         Code cleaning. Removed old, dysfunctional 'pe' option.


program define lasso2_p, rclass

	syntax namelist(min=1 max=2) [if] [in] [, lse lopt NOIsily POSTEst *]
	 
	if "`noisily'"=="" {
		local qui qui 
	}
	*
	
	local cmd `e(cmd)'
	
	// XB or R 
	if ("`lse'`lopt'"=="") { // xb or r
		if ("`cmd'"=="lasso2") {
			_lasso2_p `namelist' `if' `in', `qui' `options' `postest'
		}
		else if ("`cmd'"=="cvlasso") {
			di as "lse or lopt required."
			exit 198
		}
	}
	if ("`lse'`lopt'"!="") { // xb or r
		if ("`cmd'"=="cvlasso") {
			if ("`postest'"=="") {
				tempname m
				qui estimates store `m'
			}
			// return lasso2 results with lse or lopt
			// postest option ensures that lasso2 results are being posted
			cvlasso, `lse' `lopt' postest 
			// run predict command
			_lasso2_p `namelist' `if' `in', `qui' `options'
			if ("`postest'"=="") {
				qui estimates restore `m'
			}
		}
		else {
			di as err "lse or lopt not allowed after lasso2"
			exit 198
		}
	}  
end

// program for calculating xb/r 
program define _lasso2_p, rclass

	syntax namelist(min=1 max=2) [if] [in], ///
											///
				[XB 						/// [default]
				Residuals 					///
											///
				Lambda(numlist >0 max=1)	/// Lambda value
				LID(numlist integer max=1) 	/// Lambda ID [not applicable for CV]
											///
											///
				ols 						/// use post-OLS coefficients
											///
				APPRox 	 					/// use linear approximation [not applicable for CV]
				qui 						/// display estimation output
				POSTEst 					///
				] 							// used for info message below, should not be documented


				
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
	
	*** after cross-validation
	local command=e(cmd)
	*
	
	marksample touse, novarlist 
	
	*** warning messages
	if ("`xb'`residuals'"=="") {
		di as gr "No xb or residuals options specified. Assume xb (fitted values)."
		local xb xb
	}
	*** fe currently not supported.
	local fe = e(fe)
	if ((`fe'==1) & ("`xb'"!="")) {
		di as err "xb option currently not supported after fe option."
		exit 1
	}
	*
	
	*** obtain beta-hat
	local lcount = e(lcount)
	//local alpha = e(alpha) // retrieve elastic net parameter
	tempname betaused
	if (`lcount'==1) { // only one lambda
	
		*** syntax checks
		if ("`lambda'"!="") {
			di as error "Warning: lambda() option is ignored."
		}
		if ("`lid'"!="") {
			di as error "Warning: lid option is ignored."
		}
		if ("`approx'"!="") {
			di as error "Warning: approx option is ignored."
		}
		if ("`noisely'"!="") {
			di as err "Warning: noisely option is ignored."
		}
		
		*** for return
		local lambda = e(lambda)
	
		*** lasso or post-lasso?
		if ("`ols'"=="") {
			di as text "Use e(b) from previous lasso2 estimation (lambda=`lambda')."
			mat `betaused' = e(b)
		}
		else {
			di as text "Use e(betaOLS) from previous lasso2 estimation (lambda=`lambda')."
			mat `betaused' = e(betaOLS)
		}
		*
	
	}
	else { // list of lambdas
	
		*** syntax checks
		// either lid or lambda() option required.
		if ("`lambda'"=="") & ("`lid'"=="") {
			di as error "lambda() or lid() option required."
			exit 198
		}
		*
	
		if ("`lambda'"!="") & ("`approx'"!="") { // linear approximation
		
			di as text "Use linear approximation based on two closest lambda values."
			
			*** checks
			if (`e(alpha)'!=1) {
				di as error "Warning: Linear approximation only exact for Lasso."
			}
			if ("`ols'"!="") {
				di as error "Post option not supported with approx." 
			}	
			*
			
			*** check if lambda in range
			tempname lambdas betas
			mat `lambdas'=e(lambdamat)
			mat `betas' = e(betas)
			local lmax = e(lmax)
			local lmin = e(lmin)
			if (`lambda' < `lmin') | (`lambda' > `lmax') {
				di as error "Lamba is not in range. `lmin'<=Lambda<=`lmax' is required."
				error 198
			}
			*

			***  find smallest/largest matrix value larger/smaller 
			*** than the lambda specified by user
			local j=2
			local lminus=`lmax'
			while ((`lminus'>=`lambda') & (`j'<=`lcount')) {
				local lplusid = `j'-1
				local lminusid = `j'
				local lplus = `lambdas'[`lplusid',1]
				local lminus = `lambdas'[`lminusid',1]	
				local j=`j'+1
			}
			*

			*** extract corresponding beta vectors
			local xdim = colsof(`betas')
			tempname betaplus betaminus
			mat `betaplus' = `betas'[`lplusid',1..`xdim'] 
			mat `betaminus' = `betas'[`lminusid',1..`xdim'] 

			*** approximate beta
			local Lconstant = (`lplus'-`lambda')/(`lambda'-`lminus')
			tempname betaused
			mat `betaused' = (`betaplus'+`betaminus'*`Lconstant')/(1+`Lconstant')
			return scalar lplus=`lplus'
			return scalar lminus=`lminus'
			return scalar lplusid=`lplusid'
			return scalar lminusid=`lminusid'
			
		} 
		else if ("`lid'"!="") { // extract beta using lambda is
			
			*** syntax checks
			if ("`ols'"!="") {
				di as error "Warning: postols option not supported with lid." 
			}	
			if ("`approx'"!="") {
				di as error "Warning: approx option ignored."
			}
			*
			
			tempname lambdas betas betaused
			mat `betas' = e(betas)
			mat `lambdas'=e(lambdamat)
			local xdim = colsof(`betas')
			local lcount=rowsof(`lambdas')
			if (`lid'>`lcount') {
				di as error "lid out of range"
				error 198
			}
			mat `betaused' = `betas'[`lid',1..`xdim'] 
			//local estimator "Lasso"
			local lambda = `lambdas'[`lid',1]
			
			di as text "Use lambda with id=`lid'. lambda=`lambda'."
		
		}
		else if ("`lambda'"!="") & ("`approx'"=="") { // re-estimate
		
			*** this is used after cvlasso or lasso2 (if lcount>1)

			// store e() items 
			if ("`postest'"=="") {
				tempname origest
				estimates store `origest'
			}
			
			*** do estimation (using replay syntax)
			di as text "Re-estimate model with lambda=`lambda'."
			lasso2, newlambda(`lambda')
			
			if `e(s0)'==0 {
				di as err "No variables selected."
				exit 498
			}
			
			// get the beta used for prediction
			mat `betaused' = e(b) 
			
			if ("`postest'"=="") {
				qui estimates restore `origest'
			}
			
			//return matrix Ups = `Upsused'
		}
		else {
			di as err "internal error"
			exit 1	
		}	
	}
	*
	
	*** obtain prediction/residuals
	local depvar `e(depvar)'
	if "`depvar'"=="" {
		di as err "internal lasso2_p error. no depvar found."
	}
	tempvar xbvar
	qui matrix score `vtype'  `xbvar'= `betaused'  `if'
	if ("`xb'"!="") {
		qui gen `vtype' `varlist' = `xbvar' `if'
		label var `varlist' "Predicted values"
	}
	else if ("`residuals'"!="") {
		qui gen `vtype'  `varlist' = `depvar' - `xbvar' `if'
		label var `varlist' "Residuals"
	}
	*

	*** return beta
	return matrix beta = `betaused'
	return scalar lambda = `lambda'
end
