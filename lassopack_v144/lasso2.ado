*! lasso2 1.0.14 26oct2024
*! lassopack package 1.4.4
*! authors aa/ms

* additional notes
* return dof for 1 lambda

* eclass wrapper for elastic net & sqrt-lasso estimation
* all partialling, transformations, standardization, FE, tempvars handled here
* keeps lists of original and temp varlists
* plotting supported using "plotpath" program
* display of output handled by "DisplayPath" and "DisplayCoefs" programs
* marksample and markout used here so that e(sample) can be saved
* options relevant to eclass and saved results spelled out here
* all varlists and options passed to lassoshooting
* supports replay syntax
* lasso2 accommodates two cases: scalar lambda or list of lambdas (default)

* Updates (release date):
* 1.0.03  (30jan2018)
*         First public release.
*         Promoted to require version 13 or higher.
*         Added holdout option and related changes.
*         Replaced noprestd with prestd and related changes; added unitloadings option.
*         Recoding of cons and demeaning flags. Std loadings based on demeaned vars even with nocons.
*         partial and nocons no longer compatible.
* 1.0.04  (10feb2018)
*         Support for Sergio Correia's FTOOLS FE transform (if installed).
* 1.0.05  (4apr2018)
*         Support for information criteria added. DisplayPath program has been modified accordingly.
* 		  lic(string) and ic(string) options were added. Added various results that are stored in e().
*		  See lassoutils 1.0.08 for underlying technical changes.
*		  Added ic(none) and noic options; both suppress calculation of information criteria.
* 1.0.06  (4sep2018)
*         Fixed small typo in error message ("type 'lasso2, lic(ebic)' after estimation").
* 		  lasso2 can now be directly called with lic() option, e.g. "lasso2 y x*, lic(ebic)" [24/04/18]
* 		  Added plotpath(lnlambda) and adjusted parameters of plotlabel option. [27/04/18]
*		  Added "ebicgamma" option, see notes in lassoutils.ado. [04/09/2018]
*		  Fixed bug when fe used with partial(.).
* 1.0.07  (8nov2018)
*         Added saved value of objective function for estimation with a single lambda/alpha.
*         Added version option
*         Replaced "postest" option with name "postresults"; legacy support for postest.
* 1.0.08  (12jan2019)
*         Replace Ups terminology with Psi (penalty loadings)
*         Bug fix - FE + weights would fail if data were not sorted on xtset panel var.
* 1.0.09  (28jun2019)
*         Minor tweaks to display table (to get line to default 80 chars) and message.
*         Bug fix - plotvar(.) option would fail with FV interactions.
*         Fix to accommodate if inrange(.) syntax.
* 1.0.10  (4oct2019)
*		  Fixed bug that occured when partial was used.
*		  Fixed bug: post-lasso coefficients were sometimes not returned in e(b)
* 1.1.11  (29july2020)
*         Misc bug fixes: norecover option (should be ignored in DisplayCoefs - ignore if nothing
*           partialled-out); stdcoef+prestd (stdcoef implies prestd); fe option requires earlier
*           pre-check for markout; plot no longer assumes varabbrev=on.
*         Model s (#selected) now excludes constant; consistently excludes constant, #partialled, #FE.
*         Model df includes constant (if not FE), #partialled.
*         Added dofminus/sdofminus to capture lost degrees of freedom from partial/FE.
*         Added e(.) macros r2, df & untrunc lists lambdamat0, lmin0, lmax0.
*         Added support for psolver option.
* 1.1.12  (27sept2020)
*         stdcoef option implies norecover (can't recover std coefs for partialled out vars or constant).
*         stdcoef option implemented here; lassoutils returns both std and unstd coefs.
*         Combination of ploadings(.)+unitloadings and adaptive+unitloadings now allowed.
*         Added nostd option - synonym for unitloadings but clearer when used with lglmnet.
*         Added stdall option - standardized lambda, L1, ICs, as well as coefficients. stdall=>stdcoef.
*         Fixed bug in display in partialled-out vs factor variables.
*         e(objfn) replaces e(pmse) and e(prmse).
* 1.1.13  (5jan2024)
*         Misc code snippets to support sklearn option.
* 1.1.14  (26oct2024)
*         Bug fixe for FE with holdout sample when holdout has panel units not in the estimation sample.
*         Now drop these observations from the holdout sample and output an info message.
*         Added absorb(.) option (alternative to fe option).


program lasso2, eclass sortpreserve
	version 13
	syntax [anything] [if] [in] [,					///
			PLOTpath(string)						/// "norm" or "lambda" allowed
			PLOTVar(varlist min=1 fv ts)			/// optional subset of variables to plot
			partial(varlist fv ts)					///
			psolver(string)							/// optional solver for partialling out
			PLOTOpt(string)							/// options to pass to graph command
			PLOTLabel 								///
			POSTRESults								///
			POSTEst 								/// legacy option, now replaced by postresults
			NEWLambda(numlist >0 min=1 max=1)		///
			NEWAlpha(numlist >=0 min=1 max=1)		///
			wnorm									///
			NOPATH									/// suppress display
			displayall 								/// display zeros of beta, only applicable for one lambda
			NORecover								/// don't recover partialled-out coeffs
			long 									///
			lic(string)								/// for replay
			ic(string)								/// for lasso2 output
			NOIC 									///
			VERsion									///
			OLS 									///
			* 										///
			]
	
	local lversion 1.0.12
	local pversion 1.4.1

	if "`version'" != "" {							//  Report program version number, then exit.
		di in gr "lasso2 version `lversion'"
		di in gr "lassopack package version `pversion'"
		ereturn clear
		ereturn local version		`lversion'
		ereturn local pkgversion	`pversion'
		exit
	}
	
	*** legacy option postest replaced by postresults
	if "`postest'" != "" {
		local postresults postresults
		di as err "'postest' option has been renamed to 'postresults'. Please use 'postresults' instead."
	}
	*

	*** noic and ic(none)
	// omit calculation of IC if either noic or ic(none) is used
	if "`noic'"!="" {
		local ic none
	}
	if "`ic'"=="none" {
		local noic noic
	}
	*
	
	*** initialise local that saves whether est results are in hold
	local inhold=0
	
	*** first run of _lasso2; there is a 2nd _lasso2 call if lic() is specified
	// 3 cases:
	// (a) no replay syntax (with or w/o lic): standard case; also used by cvlasso
	// (b) replay syntax and lic(): re-run with lambda selected by IC
	// (c) replay syntax and newlambda(): re-run with newlambda, used by lasso2_p
	
	if (~replay()) {
		// no replay. estimate model.
		// this is the standard case of a fully specified model.
		// newlambda() and newalpha() are required for cvlasso.
		
		// noic and lic are incompatible
		if ("`noic'"!="") & ("`lic'"!="") {
			di as err "lic and noic are incompatible. noic ignored."
			local noic 
		}
		if ("`lic'"!="") {
			local notypemessage notypemessage
		}
		*
		
		// estimate
		_lasso2 `anything' `if' `in', `options'		///
						newlambda(`newlambda')		///	
						newalpha(`newalpha')		///
						partial(`partial')			///
						psolver(`psolver')			///
						`norecover'					///
						`noic'						///
						`ols' 
		ereturn local lasso2opt `options'
		
	}
	else if (replay()) & ("`newlambda'`newalpha'"!="") & ("`lic'"=="") {
		// replay syntax. 
		// re-estimate model with (new) single lambda (and alpha) value.
		// newlambda() and newalpha() options are undocumented.
		// this case is (primarily) intended for lasso2_p.
		
		// check for lasso2 results
		if ("`e(cmd)'"!="lasso2") {
			di as error "lasso2 estimation results not found"
			exit 301
		}
		
		// estimate
		local depvar `e(depvar)'
		local varXmodel `e(varXmodel)'
		local lasso2opt `e(lasso2opt)'
		local partial_vars `e(partial_var)'
		tempvar esample
		gen `esample' = e(sample) // ensure same sample is used
		_lasso2 `depvar' `varXmodel' `partial_vars' if `esample', 	///
								partial(`partial_vars')				///
								psolver(`psolver')					///
								`lasso2opt' 						///
								newlambda(`newlambda') 				///
								newalpha(`newalpha')				///
								`ols'
		ereturn local lasso2opt `lasso2opt' 
	}
	else if (replay()) & ("`newlambda'`newalpha'"=="") & ("`lic'"!="") {
		// replay syntax. 
		// re-estimate model with lic value.
		
		// check for lasso2 results
		if ("`e(cmd)'"!="lasso2") {
			di as error "lasso2 estimation results not found"
			exit 301
		}
		*

		// set newlambda to lambda selected by IC
		if ("`lic'"=="bic") {
			local newlambda = e(lbic)
		}
		else if ("`lic'"=="aic") {
			local newlambda = e(laic)
		}
		else if ("`lic'"=="aicc") {
			local newlambda = e(laicc)
		}
		else if ("`lic'"=="ebic") {
			//local ic ebic // is this required?
			local newlambda = e(lebic)	
		}
		else {
			di as err "lic(`lic') not allowed. Select aic, bic, aicc or ebic."
			exit 198		
		}
		*
		
		// estimate
		local depvar `e(depvar)'
		local varXmodel `e(varXmodel)'
		local lasso2opt `e(lasso2opt)'
		local partial_vars `e(partial_var)'
		local licstrupper=strupper("`lic'")
		di as text ""
		di as text "Use lambda=`newlambda' (selected by `licstrupper')."
		tempvar esample
		gen `esample' = e(sample) // ensure same sample is used
		if ("`postresults'"=="") {
			tempname model0
			_estimates hold `model0'
			local inhold = 1
		}
		_lasso2 `depvar' `varXmodel' `partial_vars' if `esample', 	///
								partial(`partial_vars')				///
								psolver(`psolver')					///
								`lasso2opt' 						///
								newlambda(`newlambda')				///
								`ols'
		ereturn local lasso2opt `lasso2opt' 
	}
	else {
		if ("`newlambda'`newalpha'"!="") & ("`lic'"!="") {
			di as error "internal lasso2 error. newlambda and lic specified."
			exit 301
		}
	}
	*
	
	*** show output if lambda is a list
	if (`e(lcount)'>1) & !missing(`e(lcount)') {
		// display should be the same as lic()
		if "`lic'"!="" {
			local ic `lic'
		}
		*
		if "`nopath'"=="" {
			DisplayPath, `wnorm' `long' ic(`ic') `notypemessage'
		}
		if ("`plotpath'`plotvar'`plotopt'"!="")  {
			plotpath, plotpath(`plotpath') 		///
					  plotvar(`plotvar')   		///
					  plotopt(`plotopt') 		///
					  `plotlabel'				///
					  `wnorm' 
		}
	}
	*
		
	*** second run of _lasso2
	// only applicable if lasso2 is called with lic option
	// re-estimate for single lambda	
	if (~replay()) & ("`lic'"!="") {
	
		* check that lambda was a list in previous estimation
		if (`e(lcount)'==1) {
			di as err "lic() only allowed if lambda() is a list."
			exit 198
		}
		* set newlambda to lambda selected by IC
		if ("`lic'"=="aic") {
			local newlambda = e(laic)
		}
		else if ("`lic'"=="bic") {
			local newlambda = e(lbic)
		} 
		else if ("`lic'"=="aicc") {
			local newlambda = e(laicc)
		}
		else if ("`lic'"=="ebic") {
			local newlambda = e(lebic)
		}
		else {
			di as err "lic(`lic') not allowed. Select aic, bic, aicc or ebic."
			exit 198		
		}

		local depvar `e(depvar)'
		local varXmodel `e(varXmodel)'
		local lasso2opt `e(lasso2opt)'
		local partial_vars `e(partial_var)'
		local licstrupper=strupper("`lic'")
		di as text ""
		di as text "Use lambda=`newlambda' (selected by `licstrupper')."
		tempvar esample
		gen `esample' = e(sample) // ensure same sample is used
		if ("`postresults'"=="") {
			tempname model0
			_estimates hold `model0'
			local inhold = 1
		}
		_lasso2 `depvar' `varXmodel' `partial_vars' if `esample', 	///
								partial(`partial_vars')				///
								psolver(`psolver')					///
								`lasso2opt' 						///
								newlambda(`newlambda') 				///
								`ols'
		ereturn local lasso2opt `lasso2opt' 
	}
	*
 
	*** Show ouput if lambda is a scalar
	if `e(lcount)'==1 {
		// norecover should be ignored in DisplayCoefs, depends only on what was estimated
		DisplayCoefs, `displayall'
		if ("`plotpath'`plotvar'`plotopt'`plotlabel'"!="") {
			di as error "Plotting only supported for list of lambda values."
			di as error "Plotting options ignored."
		}
	}
	*

	*** unhold estimation results
	if ("`postresults'"=="") & (`inhold'==1) {
		_estimates unhold `model0'
	}
end

program _lasso2, eclass sortpreserve
	version 13

	syntax varlist(numeric min=2 fv ts) [if] [in] [,	///
			NOTPen(string) 							/// list of variables not penalised
			PARtial(string)							/// string so that list can contain "_cons"
			psolver(string)							/// optional solver for partialling out
			fe										/// do within-transformation
			absorb(varname)							/// provide variable for FEs
			NOCONStant								///
			NORecover 								/// recover partialled out coefficients
			///
			/// debug & more info
			debug									/// used for debugging
			Verbose									/// pass to lassoshooting
			VVerbose								/// pass to lassoshooting
			displaynames_o(string)					/// dictionary with names of vars as supplied in varlist
			displaynames_d(string)					/// corresponding display names of vars
			pminus(int 0)							/// not used; just means rlasso code also works here
			///
			/// lambda
			Lambda(numlist >0 min=1 descending)		/// L1 penalty, either list or scalar
			lambda2(numlist >0 min=1 descending)	/// optional L2 penalty, either list or scalar
			LFactor(real 1) 						/// 
			LAMBDAMat(string)						/// alternative: specify L1 lambda as matrix
			lambda2mat(string)						/// alternative: specify L2 lambda as matrix
			NEWLambda(numlist >0 min=1  max=1 )		/// scalar
			NEWPloadings(string) 					///
			Ploadings(string) 						/// L1 norm loadings
			ploadings2(string) 						/// L2 norm loadings
			UNITLoadings							///
			lglmnet									/// use glmnet parameterization
			sklearn									/// use sklearn code
			///
			/// standardization
			PREStd 									///
			STDCoef 								/// 
			STDAll									///
			NOSTD									/// synonym for unitloadings; for use with lglmnet
			///
			/// choice of estimator
			ADAptive  								/// adaptive lasso
			ADATheta(real 1) 						/// gamma paramater for adapLASSO
			ADALoadings(string)						///
			ALPha(numlist >=0 ascending) 			/// elastic net parameter
			NEWAlpha(numlist >=0 min=1  max=1) 		///
			SQRT 									/// square-root lasso
			OLS										///
													///
			POSTAll									///
			holdout(varlist numeric min=1 max=1) 	///
													///
			NOFTOOLS								///
													///
			NOIC									///
			*										///
			]

	*** convenience option for lglmnet
	if "`nostd'"~="" {
		local unitloadings unitloadings
	}
	*
	
	*** flags
	local feflag		=("`fe'"~="")
	local absorbflag	=("`absorb'"~="")
	local debugflag		=("`debug'"~="")
	local lglmnetflag	=("`lglmnet'"~="")
	local sklearnflag	=("`sklearn'"~="")
	local prestdflag	=("`prestd'"~="")
	*
	
	** reset lambda, used for predict & replay
	if ("`newlambda'"!="") {
		local lambda = `newlambda'
	}
	if ("`newalpha'"!="") {
		local alpha = `newalpha'
	}
	if ("`newploadings'"!="") {
		tempname ploadings
		mat `ploadings' = `newploadings'
		// clear these macros
		local adaptive
		local prestd
	}
	// set alpha default. 
	local alphacount	: word count `alpha'
	if (`alphacount'==0) {
		local alpha = 1	
	}
	else if (`alphacount'>1) {
		di as err "alpha() must be a scalar."
		exit 198
	}
	// adapative - any adaptive options/variants implies adaptive
	if ("`adaloadings'"~="") | (`adatheta'~=1) {
		local adaptive adaptive
	}
	*
	
	****** syntax checks *******************************************************
	if (`alpha'>1) | (`alpha'<0) {
		di as err "alpha is out of range."
		exit 198
	}
	if ("`sqrt'"!="") & (`alpha'!=1) {
		di as error "sqrt-lasso only allowed with alpha=1."
		exit 198
	}
	local notpenpar : list notpen & partial
	if ("`notpenpar'"!="") {
		di as error "`notpenpar' listed in both notpen(.) and partial(.)"
		exit 198
	}
	local checkflag 	= ("`ploadings'"!="")+("`adaptive'"!="")
	if `checkflag'>1 {
		di as error "error: cannot combine options ploadings(.) and adaptive"
		exit 198
	}
	if `feflag' & `absorbflag' {
		di as error "incompatible options: fe and absorb(.)"
		exit 198
	}
	*
	****************************************************************************
	
	*** Record which observations have non-missing values
	marksample touse
	// need to check panel var up here
	cap xtset
	local ivar	`r(panelvar)'	// empty if not xtset
	local tvar	`r(timevar)'	// empty if not xtset
	markout `touse' `varlist' `ivar' `holdout'
	tempvar toest
	qui gen `toest' = `touse'
	if ("`holdout'"!="") {
		tempvar tohold
		qui gen `tohold' = `holdout'
		assert `tohold' == 1 | `tohold'==0 if `touse'
		qui replace `toest' = 0 if `tohold'
	}
	*

	*** absorb(.)
	// if already xtset with absorb variable, then set fe option and feflag
	// if not xtset, then xtset with absorb variable but clear xtset when exiting
	// if xtset with another variable, save xtset info and reset xt when exiting
	local xtclear=0
	local xtrestore=0
	if `absorbflag' {
		if "`ivar'"~="" {
			// data already xtset
			// update fe option/flag; this is enough if panel var = absorb var
			local fe fe
			local feflag=1
			if "`ivar'"~="`absorb'" {
				// data already xtset to another setting, so save and update xtset
				local xtivar `ivar'
				local xttvar `tvar'
				xtset `absorb'
				local ivar `absorb'
				local tvar
				local xtrestore=1
			}
		}
		else {
			// data not xtset; xtset and set flag for clearing xtset when exiting
			qui xtset `absorb'
			local ivar `absorb'
			local fe fe
			local feflag=1
			local xtclear=1
		}
	}
	*
	
	*** FEs.
	if `feflag' {
		if "`ivar'"=="" {
			di as err "Error: fe option requires data to be xtset"
			exit 459
		}
		// save current sort variables
		local sortvar	: sortedby
		// if panel id appears in holdout but not estimation sample, drop from holdout
		if "`holdout'"~="" {
			tempvar hcount hmiss
			sort `ivar' `tohold'
			qui by `ivar': gen `hcount'=sum(`tohold')
			qui by `ivar': egen `hmiss'=min(`hcount'), missing
			// if hmiss=1, means no obs for that panel in estimation sample
			qui count if `hmiss'==1
			local hobs=r(N)
			qui tab `ivar' if `hmiss'==1
			local hpanels=r(r)
			if r(N) > 0 {
				di as text "note: panels missing in the estimation sample are dropped from the holdout sample"
				di as text "      `hpanels' panels with `hobs' holdout observations dropped"
				qui replace `tohold'=0 if `hmiss'==1
				qui replace `touse'=0 if `hmiss'==1
			}
			// restore sort
			sort `sortvar'
		}
		// fe transformation may expect data to be sorted on ivar
		local sortvar_1	: word 1 of `sortvar'				// in case sorted on multiple variables
		if "`ivar'"~="`sortvar_1'" {
			sort `ivar'
		}
	}
	*

	*** sample size	
	sum `touse' if `touse', meanonly		//  will sum weight var when weights are used
	local N		= r(N)
	*
	
	*** sklearn
	// sklearn requires lglmnet
	if `sklearnflag' & ~`lglmnetflag' {
		di as res "note - sklearn option implies lglmnet option/parameterization"
		local lglmnet lglmnet
		local lglmnetflag=1
	}
	*
	
	*** lglmnet
	// glmnet treats penalty loadings and standardization separately
	if `lglmnetflag' {
		// requires either prestandardization or unit loadings
		if "`unitloadings'"=="" {
			local prestd	prestd
		}
	}
	*

	*** constant, partial, etc.
	// conmodel: constant in original model
	// consflag: constant in transformed equation to estimate
	local consmodel		=("`noconstant'"=="") & ~`feflag'	// if fe, then cons=0 & partialcons=""
	local partial		: subinstr local partial "_cons" "", all word count(local pconscount)
	local notpen		: subinstr local notpen "_cons" "", all word count(local notpenconscount)
	if (`notpenconscount'>1) {
		di as err "Warning: notpen(_cons) not supported. Constant is always partialled out."
	}
	local partialflag	= ("`partial'"~="")   	// =1 if regressor other than constant is partialled out
	local notpenflag	= ("`notpen'"~="")  	// =1 if regressor other than constant is not penalised
	local stdallflag	= ("`stdall'"~="")		// return everything in std units
	if `stdallflag' {
		local stdcoef stdcoef					// stdall => stdcoef
	}
	local stdcoefflag	= ("`stdcoef'"~="")		// return coef estimates in std units
	local prestdflag 	= ("`prestd'"~="")		// =1 if data to be pre-standardized
	// default is to use standardization loadings; overridden by ploadings, unitloadings, pre-standardization, adaptive
	local stdloadflag	= ("`ploadings'`unitloadings'`prestd'`adaptive'"=="")
	local sqrtflag 		= ("`sqrt'"!="")
	// ignore norecover if no partialled-out variables
	local parrflag		= ("`norecover'"=="") & (`partialflag' | `prestdflag')
	*
	
	// if partial list has factor vars, will need to be replaced with tempvars
	cap _fv_check_depvar `partial'
	local partialfvflag	=(_rc==198)
	// Tell estimation code if cons has been partialled out or there isn't one in the first place
	if `feflag' | `partialflag' | `prestdflag' | (~`consmodel') {
		local consflag	0
	}
	else {
		local consflag	1
	}
	*
	
	*** create main varlist and tempvars
	// remove duplicates from varlist
	// _o list is vars with original names
	fvexpand `varlist' if `touse'  
	local varlist_o	`r(varlist)'
	// check for duplicates has to follow expand
	local dups			: list dups varlist_o
	if "`dups'"~="" {
		di as text "Dropping duplicates: `dups'"
	}
	local varlist_o		: list uniq varlist_o
	*

	*** Create separate _o varlists: Y, X, notpen, partial
	// Y, X
	local varY_o		: word 1 of `varlist_o'
	local varX_o		: list varlist_o - varY_o				//  incl notpen/partial
	// notpen
	fvexpand `notpen' if `touse'
	local notpen_o		`r(varlist)'
	local dups			: list dups notpen_o
	if "`dups'"~="" {
		di as text "Dropping duplicates: `dups'"
	}
	local notpen_o		: list uniq notpen_o
	// partial
	fvexpand `partial' if `touse'
	local partial_o		`r(varlist)'
	local dups			: list dups partial_o
	if "`dups'"~="" {
		di as text "Dropping duplicates: `dups'"
	}
	local partial_o		: list uniq partial_o
	// "model" = vars without partialled-out
	local varXmodel_o	: list varX_o - partial_o
	*
	
	*** syntax checks
	// check that notpen vars are in full list
	local checklist	: list notpen_o - varX_o
	local checknum	: word count `checklist'
	if `checknum' {
		di as err "syntax error - `checklist' in notpen(.) but not in list of regressors"
		exit 198
	}
	// check that partial vars are in full list
	local checklist	: list partial_o - varX_o
	local checknum	: word count `checklist'
	if `checknum' {
		di as err "syntax error - `checklist' in partial(.) but not in list of regressors"
		exit 198
	}
	// check that ivar (FE) is not a used variable
	if `feflag' {
		fvrevar `varY_o' `varX_o', list					//  list option means we get only base vars
		local vlist `r(varlist)'
		local checklist	: list ivar - vlist
		local checknum	: word count `checklist'
		if `checknum'==0 {
			di as err "syntax error - `ivar' is xtset variable and cannot be used in model"
			exit 198
		}
	}
	// other checks
	if `pconscount' & `feflag' {
		di as err "error: incompatible options, partial(_cons) and fe"
		exit 198
	}
	if "`partial'"~="" & "`noconstant'"~="" {
		di as err "error: incompatible options, partial and nocons"
		exit 198
	}
	if `feflag' & "`noconstant'"~="" {
		di as err "error: incompatible options, fe and nocons"
		exit 198
	}
	*
	
	*** Create _t varlists: Y, X, notpen, partial
	// _o list is vars with original names
	// _t list is temp vars if transform needed, original vars if not
	if `feflag' {												//  everything needs to be transformed including partial
		local temp_ct : word count `varlist_o'
		mata: s_maketemps(`temp_ct')
		local varlist_t `r(varlist)'
	}
	else if `partialflag' | `prestdflag' {						//  everything except partial_o needs to be transformed
		local varYXmodel_o `varY_o' `varXmodel_o'
		local temp_ct : word count `varYXmodel_o'
		mata: s_maketemps(`temp_ct')
		local varYXmodel_t `r(varlist)'
		matchnames "`varlist_o'" "`varYXmodel_o'" "`varYXmodel_t'"
		local varlist_t		`r(names)'
	}
	else {														//  no transformation needed but still need temps
		fvrevar `varlist_o' if `touse'							//  fvrevar creates temps only when needed
		local varlist_t		`r(varlist)'
	}
	// dictionary is now varlist_o / varlist_t
	// now create separate _o and _t varlists using dictionary
	foreach vlist in varY varX varXmodel notpen partial {
		matchnames "``vlist'_o'" "`varlist_o'" "`varlist_t'"
		local `vlist'_t		`r(names)'						//  corresponding tempnames; always need this because of possible fvs
	}
	*

	******************* Display names ***********************************************************
	//  may be called by another program with tempvars and display names for them
	//  if display names option not used, use _o names as provided in rlasso command
	//  if display names option used, use display names matched with _o names
	//  if display names macros are empty, has no effect
	matchnames "`varY_o'" "`displaynames_o'" "`displaynames_d'"
	local varY_d		`r(names)'
	matchnames "`varXmodel_o'" "`displaynames_o'" "`displaynames_d'"
	local varXmodel_d	`r(names)'
	matchnames "`varX_o'" "`displaynames_o'" "`displaynames_d'"
	local varX_d		`r(names)'
	matchnames "`notpen_o'" "`displaynames_o'" "`displaynames_d'"
	local notpen_d		`r(names)'
	matchnames "`partial_o'" "`displaynames_o'" "`displaynames_d'"
	local partial_d		`r(names)'
	*

	*** summary varlists and flags:
	* cons			= 1 if constant, 0 if not
	* varY_o		= dep var
	* varY_t		= dep var, temp var
	* varX_o		= full, expanded set of RHS, original names, includes partial
	* varX_t		= as above but with temp names for all variables
	* varXmodel_o	= full, expanded set of RHS, original names, excludes partial
	* varXmodel_t	= as above but with temp names for all variables
	* notpen_o		= full, expanded set of not-penalized
	* notpen_t		= as above but with temp names for all variables

	//  p is calculated in lassoutils as number of model vars excluding constant
	//  here we calculate which of the model vars are omitted/base vars
	//  to provide as `pminus' to lassoutils
	//  use _o names / display names since they have info on whether var is omitted/base/etc.
	if ~`pminus' {
		foreach vn of local varXmodel_d {								//  display names
			_ms_parse_parts `vn'
			// increment pminus if model variable is MISSING
			if r(omit) {
				local ++pminus
			}
		}
	}
	//  p0 here is total number of variables provided to model EXCLUDING constant
	local p0	: word count `varXmodel_o'
	local p		=`p0'-`pminus'
	*

	******************* FE, partialling out, standardization ************************************
	//  If FE:    partial-out FEs from temp variables, then preserve,
	//            then partial-out low-dim ctrls from temp variables
	//            restore will restore all temp vars with only FEs partialled-out
	//  If no FE: leave original variables unchanged.
	//            partial-out low-dim ctrls from temp variables.
	//            if no FE/low-dim ctrls, no transform needed
	// dofminus/sdofminus captures lost degrees of freedom from FE/partialling

	local dofminus	=0										//  overwritten by FE count
	local sdofminus	=`consmodel'							//  initial "small" df count is cons/nocons
	local dmflag	=0										//  initialize demeaned flag
	if `feflag' {											//  FE-transform all variables
		fvrevar `varY_o' `varX_o' if `touse'				//  in case any FV or TS vars in _o list
		local vlist `r(varlist)'
		lassoutils `vlist',									/// call on _o list
						touse(`touse') 						///
						toest(`toest') 						///
						tvarlist(`varY_t' `varX_t')			/// overwrite/initialize these
						`noftools'							///
						fe(`ivar')							//  triggers branching to FE utility
		local dofminus	=r(N_g)								//  overwrite dofminus
		local N_g		=r(N_g)								//  N_g will be empty if no FEs
		local noftools `r(noftools)'						//  either not installed or user option
		local dmflag=1										//  data are now demeaned
		if `partialflag' {									//  And then partial out any additional vars	
			preserve										//  preserve the original values of tempvars before partialling out
			lassoutils `varY_t' `varXmodel_t',				/// _t vars have been created and filled so use here
							touse(`touse')					/// don't need tvarlist because vars already created
							toest(`toest')					/// don't need tvarlist because vars already created
							partial(`partial_t')			/// _t vars have been created and filled so use here
							partialflag(`partialflag')		/// triggers branching to partial utility
							psolver(`psolver')				/// optional choice of solver
							dmflag(1)						//  FE => mean zero
			local sdofminus	=r(rank)						//  small dof is #partialled
		}
		if `prestdflag' {
			tempname prestdY prestdX
			lassoutils `varY_t',							/// _t vars have been created and filled so use here
							touse(`touse')					/// don't need tvarlist because vars already created
							toest(`toest')					///
							std								///
							dmflag(1)						//  FE => data already mean zero
			mat `prestdY'=r(stdvec)
			lassoutils `varXmodel_t',						/// 
							touse(`touse')					/// 
							toest(`toest')					///
							std								///
							dmflag(1)						//  FE => data already mean zero
			mat `prestdX'=r(stdvec)
		}
	}
	else if `partialflag' {									//  Just partial out
		fvrevar `varY_o' `varXmodel_o' if `touse'			//  in case any FV or TS vars in _o list
		local vlist `r(varlist)'
		fvrevar `partial_o' if `touse'						//  in case any FV or TS vars in _o list
		local pvlist `r(varlist)'
		lassoutils `vlist',									/// call on _o list
						touse(`touse') 						///
						toest(`toest') 						///
						partial(`pvlist')					///
						tvarlist(`varY_t' `varXmodel_t')	/// overwrite/initialize these
						partialflag(`partialflag')			/// triggers branching to partial utility
						psolver(`psolver')					/// optional choice of solver
						dmflag(0)							//  data are not yet demeaned
		local sdofminus	=r(rank)							//  overwrite sdofminus with #partialled
		local dmflag	=1									//  data are now demeaned
		if `prestdflag' {
			tempname prestdY prestdX
			lassoutils `varY_t',							/// _t vars have been created and filled so use here
							touse(`touse')					/// don't need tvarlist because vars already created
							toest(`toest')					///
							std								///
							dmflag(1)						//  partial => already mean zero
			mat `prestdY'=r(stdvec)
			lassoutils `varXmodel_t',						/// 
							touse(`touse')					///
							toest(`toest')					///
							std								///
							dmflag(1)						//  partial => already mean zero
			mat `prestdX'=r(stdvec)
		}
	}
	else if `prestdflag' {
		tempname prestdY prestdX
		lassoutils `varY_o',								/// call on _o list
						touse(`touse')						///
						toest(`toest')						///
						std									///
						tvarlist(`varY_t')					/// overwrite/initialize these
						consmodel(`consmodel')				/// =1 => data should be demeaned
						dmflag(0)							//  data not (yet) mean zero
		mat `prestdY'=r(stdvec)
		fvrevar `varXmodel_o' if `touse'					//  in case any FV or TS vars in _o list
		local vlist `r(varlist)'
		lassoutils `vlist',									/// call on _o list
						touse(`touse')						///
						toest(`toest')						///
						std									///
						tvarlist(`varXmodel_t')				/// overwrite/initialize these
						consmodel(`consmodel')				/// =1 => data should be demeaned
						dmflag(0)							//  data not yet mean zero
		mat `prestdX'=r(stdvec)
		if `consmodel' {
			local dmflag	=1								//  if cons in model, data are now demeaned
		}
	}

 	*************** lambda to matrix **************************************************
	// lambda provided either as scalar(s) or matrix; convert to matrix
	// macro lambdamat0 will be empty if not provided
	// optional adjustment using undocumented lfactor option used for CV
	if "`lambda'`lambdamat'"!="" {
		tempname lambdamat0
		getlambdamat, lscalar(`lambda') lmatrix(`lambdamat') lfactor(`lfactor')
		mat `lambdamat0'	= r(lambdamat)
	}
	// optional L2 norm lambda
	if "`lambda2'`lambda2mat'"!="" {
		tempname lambda2mat0
		getlambdamat, lscalar(`lambda2') lmatrix(`lambda2mat') lfactor(`lfactor')
		mat `lambda2mat0'	= r(lambdamat)
	}
	*

	************* Partialling/standardization END ***********************************************
	
	*** Lasso estimation with transformed/partialled-out vars
	if "`verbose'`vverbose'"=="" {
		local quietly "quietly"							//  don't show lassoutils output
	}	

	*** Lasso estimation
	`quietly' lassoutils `varY_t',							///
					path									/// branches to _lassopath
					toest(`toest')							///
					xnames_o(`varXmodel_d')					/// display name
					xnames_t(`varXmodel_t')					///
					consflag(`consflag')					/// =0 if cons already partialled out or if no cons
					stdallflag(`stdallflag')				/// =1 if lambdas etc. are provided in the standardized metric
					dmflag(`dmflag')						/// =1 if data have been demeaned
					dofminus(`dofminus')					/// dofs lost from FEs
					sdofminus(`sdofminus')					/// dofs lost from partialling
					pminus(`pminus')						///
					notpen_o(`notpen_d') 					/// not penalised (display name)
					notpen_t(`notpen_t')					///
					lambda(`lambdamat0')					/// 
					lambda2(`lambda2mat0')					/// 
					`adaptive'								///
					adatheta(`adatheta')					///
					adaloadings(`adaloadings')				///
					`sqrt'									///
					`ols'									///
					alpha(`alpha')							///
					stdy(`prestdY')							///
					stdx(`prestdX')							///
					stdl(`stdloadflag')						/// use standardization loadings
					ploadings(`ploadings') 					/// L1 norm loadings
					ploadings2(`ploadings2') 				/// L2 norm loadings
					`verbose' `vverbose'					///
					holdout(`tohold')						///
					`noic' 									///
					`lglmnet'								/// use glmnet parameterization
					`sklearn'								/// use sklearn code
					`options'

	************* Finish up ********************************************************

	*** Create macros etc.
	local lcount	=r(lcount)
	if (`lcount'==1) { //------- scalar lambda -----------------------------------------------//	

		// message relevant for single lambda only
		if `stdcoefflag' {
			di as text "note: option stdcoef implies norecover; no constant reported" 
			// set partial-recovery flag to 0
			local parrflag 0
		}
		
		*** e-return lasso estimation results
		tempname b beta betaOLS sbeta sbetaOLS Psi stdvec
		tempname betaAll betaAllOLS sbetaAll sbetaAllOLS
		tempname lambda slambda lambda0 rmse rmseOLS objfn srmse srmseOLS sobjfn r2 df
		if "`cluster'" ~= "" {
			local N_clust		=r(N_clust)
		}
		mat `beta'			=r(beta)		//  may be empty!
		mat `betaOLS'		=r(betaOLS)		//  may be empty!
		mat `betaAll'		=r(betaAll)
		mat `sbeta'			=r(sbeta)
		mat `sbetaOLS'		=r(sbetaOLS)
		mat `betaAllOLS'	=r(betaAllOLS)
		mat `sbetaAll'		=r(sbetaAll)
		mat `sbetaAllOLS'	=r(sbetaAllOLS)
		mat `Psi'			=r(Psi)
		//*//mat `sPsi'			=r(sPsi)
		mat `stdvec'		=r(stdvec)
		scalar `lambda'		=r(lambda)
		scalar `slambda'	=r(slambda)
		scalar `lambda0'	=r(lambda0)
		scalar `rmse'		=r(rmse)		//  Lasso RMSE
		scalar `rmseOLS'	=r(rmseOLS)		//  post-Lasso RMSE
		scalar `srmse'		=r(srmse)		//  Standardized Lasso RMSE
		scalar `srmseOLS'	=r(srmseOLS)	//  Standardized post-Lasso RMSE
		scalar `r2'			=r(r2)
		scalar `df'			=r(df)
		scalar `objfn'		=r(objfn)
		scalar `sobjfn'		=r(sobjfn)
		local selected		`r(selected)'	//  EXCL NOTPEN/CONS
		local selected0		`r(selected0)'	//  INCL NOTPEN, EXCL CONS
		local s				=r(s)			//  EXCL NOTPEN/CONS; number of elements in selected
		local s0			=r(s0)			//  INCL NOTPEN, EXCL CONS; number of elements in selected0
		local k				=r(k)			//  number of all variables in beta INCL NOTPEN/CONS (if present)
		local p0			=r(p0)			//  number of all variables in betaAll INCL NOTPEN/CONS (if present)
		//*//local clustvar		`r(clustvar)'
		//*//local robust		`r(robust)'
		//*//local center		=r(center)
		local sqrtflag 		=r(sqrt)
		local alpha			=r(alpha)
		local olsflag 		= r(olsflag)
		local method		`r(method)'		//  lasso or sqrt-lasso
		local niter			=r(niter)
		local maxiter		=r(maxiter)
		//*//local nupsiter		=r(nupsiter)
		//*//local maxupsiter	=r(maxupsiter)
		// issue warning if lasso max iteration limit hit
		if `niter'==`maxiter' {
			di as text "Warning: reached max shooting iterations w/o achieving convergence."
		}
		// fix depvar (rownames) of beta vectors to use _o (or _d if display names provided) not _t
		mat rownames `beta'			= `varY_d'
		mat rownames `betaOLS'		= `varY_d'
		mat rownames `betaAll'		= `varY_d'
		mat rownames `betaAllOLS'	= `varY_d'
		mat rownames `sbeta'		= `varY_d'
		mat rownames `sbetaOLS'		= `varY_d'
		mat rownames `sbetaAll'		= `varY_d'
		mat rownames `sbetaAllOLS'	= `varY_d'
		// used below
		if `k'>0 {										// cnames will be empty if k=0
			local cnames_o	: colnames `beta'
			fvstrip `cnames_o'							//  colnames may insert b/n/o operators - remove
			local cnames_o	`r(varlist)'
			matchnames "`cnames_o'" "`varlist_o'" "`varlist_t'"
			local cnames_t	`r(names)'
		}
		if `debugflag' {
			di as text "selected: `selected'"
			di as text "Returned results from lassoutils:"
			return list
			di as text "beta and betaOLS:"
			mat list `beta'
			mat list `betaOLS'
		}
		*

		*********** Get coeff estimates for partialled-out vars. ********************
		if `feflag' & `partialflag' {					//  FE case and there are partialled-out notpen vars
			restore										//  Restores dataset with tempvars after FE transform but before notpen partialled out
		}
		if (`partialflag' | (`prestdflag' & `consmodel')) & (`parrflag') {	//  standardization removes constant so must enter for that
			if `feflag' {
				local depvar `varY_t'					//  use FE-transformed depvar and X vars
				local scorevars `cnames_t'
			}
			else {
				local depvar `varY_o'					//  use original depvar and X vars
				local scorevars `cnames_o'
			}
			lassoutils `depvar',						///
				unpartial								///
				touse(`toest')							///
				beta(`beta')							///
				scorevars(`scorevars')					///
				partial(`partial_t')					///
				names_o(`varlist_o')					/// dictionary
				names_t(`varlist_t')					///	dictionary
				consmodel(`consmodel')
			mat `beta'			= r(b)
			mat `betaAll'		= `betaAll', r(bpartial)
			lassoutils `depvar',						///
				unpartial								///
				touse(`toest')							///
				beta(`betaOLS')							///
				scorevars(`scorevars')					///
				partial(`partial_t')					///
				names_o(`varlist_o')					/// dictionary
				names_t(`varlist_t')					///	dictionary
				consmodel(`consmodel')
			mat `betaOLS'		= r(b)
			mat `betaAllOLS'	= `betaAllOLS', r(bpartial)
			// finish by adding partialled-out to k
			local k				=colsof(`beta')
		}
		*	

		*** Post results
		if `stdcoefflag' {
			// post standardized coeffs
			if "`ols'"=="" & "`postall'"=="" {
				mat `b' = `sbeta'		//  selected post-lasso coeffs by default
			}
			else if "`ols'"!="" & "`postall'"=="" {
				mat `b' = `sbetaOLS'
			}
			else if "`ols'"=="" & "`postall'"!="" {
				mat `b' = `sbetaAll'
			}
			else {
				mat `b' = `sbetaAllOLS'
			}
		}
		else {
			// post unstandardized coeffs
			if "`ols'"=="" & "`postall'"=="" {
				mat `b' = `beta'		//  selected post-lasso coeffs by default
			}
			else if "`ols'"!="" & "`postall'"=="" {
				mat `b' = `betaOLS'
			}
			else if "`ols'"=="" & "`postall'"!="" {
				mat `b' = `betaAll'
			}
			else {
				mat `b' = `betaAllOLS'
			}
		}

		if `k'==0 {				//  no vars selected
			ereturn post    , obs(`N') depname(`varY_d') esample(`toest')		// display name
		}
		else {
			ereturn post `b', obs(`N') depname(`varY_d') esample(`toest')		// display name
		}
		
		// additional returned results
		ereturn local noftools		`noftools'
		ereturn local postall		`postall'
		ereturn scalar niter		=`niter'
		ereturn scalar maxiter		=`maxiter'
		//*//ereturn scalar nupsiter		=`nupsiter'
		//*//ereturn scalar maxupsiter	=`maxupsiter'
		//*//ereturn local robust		`robust'
		ereturn local ivar			`ivar'
		ereturn local selected		`selected'			//  selected only
		ereturn local varXmodel		`varXmodel_d'		//  display name
		ereturn local varX			`varX_d'			//  display name
		ereturn local method		`method'
		ereturn local predict		lasso2_p
		ereturn local cmd			lasso2
		ereturn scalar pminus		=`pminus'
		//*//ereturn scalar center		=`center'
		ereturn scalar stdcoef		=`stdcoefflag'
		ereturn scalar cons			=`consmodel'
		ereturn scalar slambda		=`slambda'
		//*//ereturn scalar lambda0		=`lambda0'
		ereturn scalar lambda		=`lambda'

		if "`N_clust'" ~= "" {
			ereturn local clustvar	`clustvar'
			ereturn scalar N_clust	=`N_clust'
		}
		if "`N_g'" ~= "" {
			ereturn scalar N_g		=`N_g'
		}
		ereturn scalar fe			=`feflag'
		ereturn local  absorb		`absorb'
		ereturn scalar rmse			=`rmse'
		ereturn scalar rmseOLS		=`rmseOLS'
		ereturn scalar srmse		=`srmse'
		ereturn scalar srmseOLS		=`srmseOLS'
		ereturn scalar r2			=`r2'
		ereturn scalar df			=`df'
		ereturn scalar objfn		=`objfn'
		ereturn scalar sobjfn		=`sobjfn'
		ereturn scalar p			=`p'
		ereturn scalar k			=`k'				//  number of all estimated coefs INCLUDING PARTIALLED-OUT AND CONSTANT
		ereturn scalar s			=`s'				//  number of selected

		ereturn matrix stdvec		=`stdvec'
		//*//ereturn matrix sPsi 		=`sPsi'
		ereturn matrix Psi 			=`Psi'
		ereturn matrix sbetaAllOLS	=`sbetaAllOLS'
		ereturn matrix sbetaAll		=`sbetaAll'
		ereturn matrix sbetaOLS		=`sbetaOLS'
		ereturn matrix sbeta		=`sbeta'
		ereturn matrix betaAllOLS	=`betaAllOLS'
		ereturn matrix betaAll		=`betaAll'
		ereturn matrix betaOLS		=`betaOLS'
		ereturn matrix beta			=`beta'

		// lasso2-specific:
		// constant is always considered partialled out
		if (`consmodel') {
			local selected0		`selected0' _cons
			local partial_d		`partial_d' _cons		//  display name
		}
		local notpen_ct		: word count `notpen_d'		//  number of notpen INCLUDING CONSTANT (if not partialled-out)
		local partial_ct	: word count `partial_d'	//  number of partialled-out INCLUDING CONSTANT
		ereturn local selected0		`selected0'			//  selected or notpen, INCL CONS
		local s0			: word count `selected0'	//  overwrite s0 to include partialled-out/constant count
		ereturn scalar s0			=`s0'				//  number of selected or notpen, INCL CONS
		ereturn local partial		`partial_d'			//  display name
		ereturn local partial_var 	`partial'
		ereturn local notpen		`notpen_d'			//  display name
		ereturn scalar notpen_ct	=`notpen_ct'		//  number of notpen INCLUDING CONS (unless partialled-out)
		ereturn scalar partial_ct	=`partial_ct'		//  number of partialled-out INCLUDING CONS
		*		
		
		*** more lasso2 ereturns
		ereturn scalar alpha		=`alpha'
		ereturn scalar fe 			=`feflag'
		ereturn local  absorb		`absorb'
		ereturn scalar sqrt  		= `sqrtflag'
		ereturn scalar prestd		= `prestdflag'
		ereturn scalar ols 			= `olsflag'
		ereturn scalar adaptive		= "`adaptive'"!=""
		ereturn scalar lcount		=`lcount'

		ereturn scalar dofminus		=`dofminus'
		ereturn scalar sdofminus	=`sdofminus'

	}
	else if (`lcount'>1) { //------- list of lambdas -------------------------------------------------//

		*** Create macros etc.
		local nobs		=r(N)
		local lcount	=r(lcount)
		local method	`r(method)'	//  "lasso", "sqrt-lasso", "elastic net"
		local alpha 	= r(alpha)
		local sqrt		= r(sqrt)
		tempname Psi sups stdvec  
		mat `Psi' = r(Psi)
		mat `stdvec'	= r(stdvec)
		local olsflag 	= r(olsflag)
		local xvars		`varXmodel_o'
		local depvar	`r(depvar)'
		local lmin		=r(lmin)
		local lmin0		=r(lmin0)
		local lmax		=r(lmax)
		local lmax0		=r(lmax0)
		tempname betas sbetas lambdas lambdas0 slambdas slambdas0
		tempname l1norm sl1norm wl1norm swl1norm Rsquared dof shat shat0
		mat `betas'		=r(betas)
		mat `sbetas'	=r(sbetas)
		mat `lambdas'	=r(lambdalist)
		mat `lambdas0'	=r(lambdalist0)
		mat `slambdas'	=r(slambdalist)
		mat `slambdas0'	=r(slambdalist0)
		mat `l1norm'	=r(l1norm)
		mat `sl1norm'	=r(sl1norm)
		mat `wl1norm'	=r(wl1norm)
		mat `swl1norm'	=r(swl1norm)
		mat `Rsquared'	=r(Rsquared)
		mat `dof'		=r(dof)
		mat `shat'		=r(shat)
		mat `shat0'		=r(shat0)
		local s			=r(s)
		if ("`holdout'"!="") {
			tempname mspe0
			mat `mspe0' = r(mspe)
			//mat list `mspe0'
		}	
		else if "`noic'"=="" {
			tempname rss ess tss rsq
			tempname aicmat bicmat aiccmat ebicmat saicmat sbicmat saiccmat sebicmat IC sIC
			mat `rss'			= r(rss)
			mat `ess'			= r(ess)
			mat `tss'			= r(tss)
			mat `rsq'			= r(rsq)	
			// aic
			local laicid		= r(laicid)
			mat `aicmat'		= r(aic)
			local aicmin		= r(aicmin)
			local laic			= `lambdas'[`laicid',1]
			mat `saicmat'		= r(saic)
			local saicmin		= r(saicmin)
			local slaic			= `slambdas'[`laicid',1]
			// aicc
			local laiccid		= r(laiccid)
			mat `aiccmat'		= r(aicc)
			local aiccmin		= r(aiccmin)
			local laicc			= `lambdas'[`laiccid',1]
			mat `saiccmat'		= r(saicc)
			local saiccmin		= r(saiccmin)
			local slaicc		= `slambdas'[`laiccid',1]
			// bic
			local lbicid		= r(lbicid)
			mat `bicmat'		= r(bic)
			local bicmin		= r(bicmin)
			local lbic			= `lambdas'[`lbicid',1]
			mat `sbicmat'		= r(sbic)
			local sbicmin		= r(sbicmin)
			local slbic			= `slambdas'[`lbicid',1]
			// ebic
			local ebicgamma		= r(ebicgamma)
			local lebicid		= r(lebicid)
			mat `ebicmat'		= r(ebic)
			local ebicmin		= r(ebicmin)
			local lebic			= `lambdas'[`lebicid',1]
			mat `sebicmat'		= r(sebic)
			local sebicmin		= r(sebicmin)
			local slebic		= `slambdas'[`lebicid',1]
			
			mat `IC'			= `aicmat', `aiccmat', `bicmat', `ebicmat'
			mat colnames `IC'	= "AIC" "AICC" "BIC" "EBIC"
			mat `sIC'			= `saicmat', `saiccmat', `sbicmat', `sebicmat'
			mat colnames `sIC'	= "sAIC" "sAICC" "sBIC" "sEBIC"
		}

		// standardization removed constant if present so k is number of elements in vectors
		local cnames_o	: colnames `betas'
		local pmodel	: word count `cnames_o' 

		*********** Get coeff estimates for constant and/or unpenalized vars. ********************
		if (!`feflag') & ((`partialflag') | (`consmodel' & `prestdflag')) & (`parrflag') {	
		// recovery of constant 
		// only supported if there are no other partialled out regressors and no FE

			fvstrip `cnames_o'					//  colnames may insert b/n/o operators - remove
			local cnames_o	`r(varlist)'
			matchnames "`cnames_o'" "`varlist_o'" "`varlist_t'"
			local cnames_t	`r(names)'

			tempname bi binew betasnew
			forvalues i= 1/`lcount' {
				mat `bi' = `betas'[`i',1..`pmodel']

				lassoutils `varY_o',						///
					unpartial								///
					touse(`toest')							///
					beta(`bi')								///
					scorevars(`cnames_o')					///
					partial(`partial_t')					///
					names_o(`varlist_o')					/// dictionary
					names_t(`varlist_t')					///	dictionary
					consmodel(`consmodel')

				mat `binew' = r(b)
				if `i'==1 {
					mat `betasnew' = `binew'
				}
				else {
					mat `betasnew' = (`betasnew' \ `binew')
				}
			}
			mat `betas' = `betasnew'

		}		
		*
	
		*** ereturns
		ereturn post    , obs(`N') depname(`varY_d') esample(`toest')	//  display name
		ereturn scalar stdcoef		=`stdcoefflag'
		ereturn scalar stdall		=`stdallflag'
		ereturn scalar cons 		=`consmodel'
		ereturn scalar fe 			=`feflag'
		ereturn local  absorb		`absorb'
		ereturn scalar alpha		=`alpha'
		ereturn scalar sqrt  		=`sqrtflag'
		ereturn scalar ols	 		=`olsflag' 
		ereturn scalar adaptive		="`adaptive'"!=""
		ereturn scalar p 			=`p'
		local notpen_ct		: word count `notpen_d'		//  number of notpen INCLUDING CONSTANT (if not partialled-out)
		local partial_ct	: word count `partial_d'	//  number of partialled-out INCLUDING CONSTANT
		ereturn scalar notpen_ct	=`notpen_ct'		//  number of notpen INCLUDING CONS (unless partialled-out)
		ereturn scalar partial_ct	=`partial_ct'		//  number of partialled-out INCLUDING CONS
		ereturn scalar prestd		=`prestdflag'
		ereturn scalar lcount 		=`lcount'
		ereturn scalar lmax			=`lmax'
		ereturn scalar lmax0		=`lmax0'
		ereturn scalar lmin			=`lmin'
		ereturn scalar lmin0		=`lmin0'
		ereturn local noftools		`noftools'
		ereturn local method		`method'
		ereturn local predict		lasso2_p
		ereturn local cmd			lasso2
		ereturn local varXmodel		`varXmodel_d'		//  display name
		ereturn local varX			`varX_d'			//  display name
		ereturn local partial		`partial_d'			//  display name
		ereturn local partial_var 	`partial'
		ereturn local notpen		`notpen_d'			//  display name
		ereturn matrix l1norm		=`l1norm'
		ereturn matrix sl1norm		=`sl1norm'
		ereturn matrix wl1norm		=`wl1norm'
		ereturn matrix swl1norm		=`swl1norm'
		ereturn matrix Psi			=`Psi'
		ereturn matrix betas		=`betas' 	 
		ereturn matrix sbetas		=`sbetas' 	 
		ereturn matrix dof			=`dof'
		ereturn matrix s			=`shat'
		ereturn matrix s0			=`shat0'
		ereturn matrix lambdamat	=`lambdas'
		ereturn matrix lambdamat0	=`lambdas0'
		ereturn matrix slambdamat	=`slambdas'
		ereturn matrix slambdamat0	=`slambdas0'
		ereturn scalar dofminus		=`dofminus'
		ereturn scalar sdofminus	=`sdofminus'
		ereturn matrix stdvec		=`stdvec'
		if ("`holdout'"!="") {
			ereturn matrix mspe = `mspe0'
		}
		else if "`noic'"=="" {
			//ereturn scalar laicid = `laicid'
			//ereturn scalar laiccid = `laiccid'
			//ereturn scalar lbicid = `lbicid'
			//ereturn scalar lebicid = `lebicid'
			ereturn scalar saicmin		= `saicmin'
			ereturn scalar saiccmin		= `saiccmin'
			ereturn scalar sbicmin		= `sbicmin'
			ereturn scalar sebicmin 	= `sebicmin'
			ereturn scalar aicmin		= `aicmin'
			ereturn scalar aiccmin		= `aiccmin'
			ereturn scalar bicmin		= `bicmin'
			ereturn scalar ebicmin		= `ebicmin'
			ereturn scalar ebicgamma	= `ebicgamma'
			ereturn matrix rss			= `rss'
			ereturn matrix ess			= `ess'
			ereturn matrix tss			= `tss'
			ereturn matrix rsq			= `rsq'
			ereturn scalar slaic		= `slaic'
			ereturn scalar slaicc		= `slaicc'
			ereturn scalar slbic		= `slbic'
			ereturn scalar slebic	 	= `slebic'
			ereturn scalar laic			= `laic'
			ereturn scalar laicc		= `laicc'
			ereturn scalar lbic			= `lbic'
			ereturn scalar lebic	 	= `lebic'
			// ereturn matrix aic		= `aicmat'
			// ereturn matrix bic		= `bicmat'
			// ereturn matrix aicc		= `aiccmat'
			// ereturn matrix ebic		= `ebicmat'
			// ereturn matrix saic		= `saicmat'
			// ereturn matrix saicc		= `saiccmat'
			// ereturn matrix sbic		= `sbicmat'
			// ereturn matrix sebic		= `sebicmat'
			
			ereturn matrix IC			= `IC'
			ereturn matrix sIC			= `sIC'
		}

	}
	
	// finish up
	// if we xtset the data to support absorb, remove this
	if `xtclear' {
		qui xtset, clear
	}
	// if we changed xtset to support absorb, restore to original settings
	if `xtrestore' {
		qui xtset `xtivar' `xttvar'
	}
	
end


program define plotpath

	syntax [anything] [, 	plotvar(string)		///
							plotpath(string)	///
							plotopt(string)		///
							plotlabel			///
							wnorm				///
							]
	
	version 12
	
	if ("`plotpath'"!="") {
	
		if (("`plotpath'"!="lambda") & ("`plotpath'"!="norm") & ("`plotpath'"!="lnlambda")) {
			di as err "Plotpath() allows 'lambda', `lnlambda' or 'norm'."
			error 198
		}
	
		// Contents of b matrix and lambda vector made into Stata variables for plotting.
		// Varnames taken from colnames of b matrix.
		// Strip out constant (if it's there) since creating a variable called _cons not alllowed.
		// If `plotvar' macro is empty, graph all regressors.
		tempname lambdas l1norm 
		tempvar touse
		gen `touse'=e(sample)
		
		// standardized or not
		if e(stdall) {
			local s s
		}

		mat b=e(`s'betas)
		mat `lambdas' = e(`s'lambdamat)

		if "`wnorm'"=="" {
			mat `l1norm' = e(`s'l1norm)
		}
		else {
			mat `l1norm' = e(`s'wl1norm)
		}
		local lcount = e(lcount)
		local cons = e(cons)
		if `cons' {
			local rb1 = colsof(b) - 1	//  constant is in the last column
			mat b = b[1...,1..`rb1']
		}
		local bnames : colnames b
		fvstrip `bnames'				//  annoying - Stata inserts b/n etc. in first factor variable etc.
		local bnames `r(varlist)'
		// process pv names
		if "`plotvar'"=="" {
			local pvnames `bnames'		//  plot all
		}
		else {							//  plot user-provided
			fvstrip `plotvar' if `touse', expand dropomit
			local pvnames	`r(varlist)'
		}
		foreach pvn in `pvnames' {		//  unab one-by-one to standardise, get i prefix etc.
			fvunab pvn_unab	: `pvn'
			local pvnames_unab `pvnames_unab' `pvn_unab'
		}
		// process b names
		foreach bn in `bnames' {		//  unab one-by-one to standardise, get i prefix etc.
			fvunab bn_unab	: `bn'
			local bnames_unab `bnames_unab' `bn_unab'
		}
		// now that unabbreviated varlists are prepared, check that plotvars are in regressors
		if "`plotvar'"~="" {
			local nplotvarcheck	 : list pvnames_unab - bnames_unab
			if ("`nplotvarcheck'"!="") {								
				di as error "Variable(s) `nplotvarcheck' of plotvar() not listed as regressor(s)." 
				exit 198
			}
		}
		// in case there are any . or # operators included, change to "_" or "__"
		local bnames	: subinstr local bnames_unab "." "_", all count(local numsubs)
		local pvnames	: subinstr local pvnames_unab "." "_", all count(local numsubs)
		local bnames	: subinstr local bnames "#" "__", all count(local numsubs)
		local pvnames	: subinstr local pvnames "#" "__", all count(local numsubs)
		// check for max number of variables to plot
		local npv : word count `pvnames'
		if `npv' >= 100 {
			di as err "Error: lassopath can graph at most 99 regressors"
			di as err "       use plotvar(.) option to specify subset of regressors"
			exit 103
		}

		// create graphing data and then plot
 		preserve						//  do this here so that above vars exist
		clear
		qui svmat b
		foreach var of varlist b* {
			tokenize `bnames'
			rename `var' `1'
			mac shift
			local bnames `*'
		}
		if "`plotpath'"=="lnlambda" {
			qui svmat `lambdas', names("lambda")
			replace lambda1 = ln(lambda1)
			if ("`plotlabel'"!="") {
				local txt
				local xcoord = lambda1[1]-abs(lambda1[_N]-lambda1[1])*1.03
				local xcoordminus = lambda1[1]-abs(lambda1[_N]-lambda1[1])*1.1
				foreach var of varlist `pvnames' {	
					local ycoord = `var'[_N] 
					local vn = abbrev("`var'",8)
					local txt `txt' text(`ycoord' `xcoord' `"`vn'"', place(w) just(left) size(small))	
				}
				local yscalealt yscale(alt)
				local xscale xscale(range(`xcoordminus'))		//  extend plot area on left to allow room for varnames
			}  
			twoway line `pvnames' lambda1, `plotopt' `txt' `yscalealt' xtit("ln(Lambda)") `graphr' `xscale'
		}
		else if "`plotpath'"=="lambda" {
			qui svmat `lambdas', names("lambda")
			if ("`plotlabel'"!="") {
				local txt
				local xcoord = -abs(lambda1[1])*0.03
				local xcoordminus = -abs(lambda1[1])*0.15
				foreach var of varlist `pvnames' {	
					local ycoord = `var'[_N] 
					local vn = abbrev("`var'",8)
					local txt `txt' text(`ycoord' `xcoord' `"`vn'"', place(w) just(left) size(small))	
				}
				local yscalealt yscale(alt)
				local xscale xscale(range(`xcoordminus'))		//  extend plot area on left to allow room for varnames
			}  
			twoway line `pvnames' lambda1, `plotopt' `txt' `yscalealt' xtit("Lambda") `graphr' `xscale'
		}
		else {
			qui svmat `l1norm', names("l1norm")
 			sort l1norm1
			if ("`plotlabel'"!="") {
				local txt
				local xcoord = l1norm1[_N]*1.02		//  extend plot area on right to allow room for varnames
				local xcoordplus = l1norm1[_N]*1.1
				foreach var of varlist `pvnames' {
					local ycoord = `var'[_N]
					local vn = abbrev("`var'",8)
					local txt `txt' text(`ycoord' `xcoord' `"`vn'"', place(e) just(left) size(small))
				}
				local xscale xscale(range(`xcoordplus'))
			}
			if "`wnorm'"=="" {
				local xtitle L1 Norm
			}
			else {
				local xtitle Weighted L1 Norm
			}

			line `pvnames' l1norm, `plotopt' `txt' xtit("`xtitle'") `graphr' `xscale'
		}

 		restore
	}
	*
	
end


// Display table of path with knots, lambda, vars added/removed etc.
program define DisplayPath
	//syntax [anything] [, stdcoef(int 0)]
	syntax [anything] [, wnorm long ic(string) NOTYpemessage ]

	version 12
	tempname betas r1 r2 vnames d addedM removedM lambdas dof l1norm vnames0 allselected allsec
	tempname all_ic icmat rsq
	
	***** information criteria *************************************************
	if ("`ic'"=="") {
		local ic ebic
	}
	if ("`ic'"!="aic") & ("`ic'"!="bic") & ("`ic'"!="aicc") & ("`ic'"!="ebic") & ("`ic'"!="none") {
		di as err "Option ic(`ic') not allowed. Using the default ic(ebic)."
		local ic ebic
	}
	if e(stdall) {
		// prefix for standardized IC
		local s s
	}
	
	// has all IC or standardized IC vectors
	mat `all_ic'		= e(`s'IC)
	
	if ("`ic'"=="ebic") {
		mat `icmat' 	=`all_ic'[.,"`s'EBIC"]
		local icmin 	=e(`s'ebicmin)
		local ic EBIC
	}
	else if ("`ic'"=="aic") {
		mat `icmat' 	=`all_ic'[.,"`s'AIC"]
		local icmin 	=e(`s'aicmin)
		local ic AIC
	}
	else if ("`ic'"=="bic") {
		mat `icmat' 	=`all_ic'[.,"`s'BIC"]
		local icmin 	=e(`s'bicmin)
		local ic BIC
	}
	else if ("`ic'"=="aicc") {
		mat `icmat' 	=`all_ic'[.,"`s'AICC"]
		local icmin 	=e(`s'aiccmin)	
		local ic AICc
	}
	else {
		mat `icmat' 	=.
		local icmin 	=.
		local ic IC
	}
	****************************************************************************
	
	mat `lambdas'	=e(`s'lambdamat)

	mat `dof'		=e(s)
	mat `rsq' 		=e(rsq)
	mata: `vnames'	=st_matrixcolstripe("e(betas)")		// starts as k x 2
	mata: `vnames'	=(`vnames'[.,2])'					// take 2nd col and transpose into row vector
	mata: `betas'	=st_matrix("e(betas)")
	mata: `r1'		=(`betas'[1,.] :!= 0)
	mata: `addedM'	=select(`vnames', (`r1' :== 1))
	mata: st_local("added",invtokens(`addedM'))
	local knot		=0
	mata: `r1'		=J(1,cols(`betas'),0)
	di
	if "`wnorm'"=="" {
		mat `l1norm' = e(`s'l1norm)
		di as text "  Knot{c |}  ID     Lambda    s      L1-Norm        `ic'" _c
		di as text _col(58) "R-sq   {c |} Action"
	}
	else {
		mat `l1norm' = e(`s'wl1norm)
		di as text "  Knot{c |}  ID     Lambda    s     wL1-Norm        `ic'" _c
		di as text _col(58) "R-sq   {c |} Action"  
	}
	// line is always 80 chars
	di as text "{hline 6}{c +}{hline 57}{c +}{hline 15}"
	forvalues i=1/`e(lcount)' {
		mata: `r2'			=(`betas'[`i',.] :!= 0)
		mata: `d'			=`r2'-`r1'
		mata: `addedM'		=select(`vnames',(`d' :== 1))
		mata: `allselected' = (sum(`r2':==0))==0 // = 1 if all selected
		mata: st_numscalar("`allsec'",`allselected')
		mata: st_local("added",invtokens(`addedM'))
		mata: `removedM'	=select(`vnames',(`d' :== -1))
		mata: st_local("removed",invtokens(`removedM'))
		if ("`added'`removed'" ~= "") | ("`long'"!="") { 
			if ("`added'`removed'" ~= "") {
				local ++knot
				di as res %6.0f `knot' _c
			}
			di as text _col(7) "{c |}" _c
			di as res %4.0f `i' _c
			di as res _col(13) %10.5f el(`lambdas',`i',1) _c
			di as res _col(25) %4.0f el(`dof',`i',1) _c
			di as res _col(31) %10.5f el(`l1norm',`i',1) _c
			di as res _col(43) %11.5f el(`icmat',`i',1) _c			//  can be negative so add a space
			// changed from 10^-5 to 10^-12
			if ("`long'"!="") & (reldif(`icmin',el(`icmat',`i',1))<10^-12) & ("`icmin'"!=".") {
				di as text "*" _c
			}
			else {
				di as text " " _c
			}
			di as res _col(56) %7.4f el(`rsq',`i',1) _c
			di as text _col(65) "{c |}" _c
			// clear macro
			macro drop _dtext
			if (`i'==1) & (`allsec') {
				local dtext All selected.
			}
			else {
				if "`added'" ~= "" {
					local dtext Added `added'.
				}
				if "`removed'" ~= "" {
					local dtext `dtext' Removed `removed'.
				}
			}
			DispVars `dtext', _lc1(7) _lc2(65) _col(67)
		}
		mata: `r1'		=`r2'
	}
	local iclower = strlower("`ic'")
	if ("`long'"=="") {
		di as text "Use {bf:long} option for full output."
	}
	else if ("`ic'"!="IC") {
		di as text "{helpb lasso2##aicbic:*}indicates minimum `ic'."
	}
	if ("`ic'"!="IC") & ("`notypemessage'"=="") {
		di as text											///
			"{p 0 6 2}"										///
			"Type e.g. {stata lasso2, lic(`iclower')}"		///
			" to run the model selected by `ic'."			///
			"{p_end}"
	}
	mata: mata drop `betas' `r1' `r2' `vnames' `d' `addedM' `removedM'
	
end

// Display varlist with specified indentation
program define DispVars
	version 11.2
	syntax [anything] [, _col(integer 15) _lc1(integer 0) _lc2(integer 0) ]
	local maxlen = c(linesize)-`_col'
	local len = 0
	local first = 1
	foreach vn in `anything' {
		local vnlen		: length local vn
		if `len'+`vnlen' > `maxlen' {
			di
			local first = 1
			local len = `vnlen'
			if `_lc1' {
				di as text _col(`_lc1') "{c |}" _c
			}
			if `_lc2' {
				di as text _col(`_lc2') "{c |}" _c
			}
		}
		else {
			local len = `len'+`vnlen'+1
		}
		if `first' {
			local first = 0
			di as res _col(`_col') "`vn'" _c
			}
		else {
			di as res " `vn'" _c
		}
	}
* Finish with a newline
	di
end

// version 2020-09-11
prog DisplayCoefs

	syntax	,								///
		[									///
		displayall							///  full coef vector in display (default=selected only)
		varwidth(int 17)					///
		]
	
	local cons			=e(cons)

	local alist			: colnames e(betaAll)
	fvstrip `alist'
	local alist			`r(varlist)'
	local plist			`e(partial)'
	fvstrip `plist'
	local plist			`r(varlist)'
	// vplist is partialled-out vars that don't appear as regressors
	// if this is empty, then they appear as regressors and should be displayed
	// if this not empty, then there is no separate partialled-out section to display
	// stdcoef overrides this - no partialled-out including constant
	local vplist		: list plist - alist
	if e(stdcoef) {
		local partial
		local partial_ct	=0
	}
	else if `: word count `vplist'' == 0 {
		// betaAll includes partialled-out coefs to display
		local partial		`e(partial)'
		local partial_ct	=e(partial_ct)
	}
	else {
		local partial
		local partial_ct	=0
	}
	
	// varlists
	local selected		`e(selected)'
	fvstrip `selected'
	local selected		`r(varlist)'
	local notpen		`e(notpen)'
	fvstrip `notpen'
	local notpen		`r(varlist)'
	local selected0		`e(selected0)'
	fvstrip `selected0'
	local selected0		`r(varlist)'
	// coef vectors
	tempname beta betaOLS
	if "`displayall'"~="" {						//  there must be some vars specified even if nothing selected
		if e(stdcoef) {
			mat `beta'		=e(ssbetaAll)
			mat `betaOLS'	=e(ssbetaAllOLS)
		}
		else {
			mat `beta'		=e(betaAll)
			mat `betaOLS'	=e(betaAllOLS)
		}
		local col_ct	=colsof(`beta')
		local vlist		: colnames `beta'
		local vlistOLS	: colnames `betaOLS'
		local baselevels baselevels
	}
	else if e(k)>0 {							//  display only selected, but only if there are any
		if e(stdcoef) {
			mat `beta'		=e(sbeta)
			mat `betaOLS'	=e(sbetaOLS)
		}
		else {
			mat `beta'		=e(beta)
			mat `betaOLS'	=e(betaOLS)
		}
		local col_ct	=colsof(`beta')
		local vlist		: colnames `beta'
		local vlistOLS	: colnames `betaOLS'
	}
	else {										//  nothing selected, zero columns in beta
		local col_ct	=0
	}
	if e(k)>0 {
		_ms_build_info `beta' if e(sample)
		_ms_build_info `betaOLS' if e(sample)
	}

	*** (Re-)display coefficients including constant/partial
	local varwidth1		=`varwidth'+1
	local varwidth3		=`varwidth'+3
	local varwidth4		=`varwidth'+4
	local varwidthm7	=`varwidth'-7
	local varwidthm13	=`varwidth'-13
	di
	di as text "{hline `varwidth1'}{c TT}{hline 32}"
	if "`e(method)'"=="sqrt-lasso" {
		di as text _col(`varwidthm7') "Selected {c |}      Sqrt-lasso   Post-est OLS"
	}
	else if "`e(method)'"=="ridge" {
		di as text _col(`varwidthm7') "Selected {c |}           Ridge   Post-est OLS"
	}
	else if "`e(method)'"=="elastic net" {
		di as text _col(`varwidthm7') "Selected {c |}     Elastic net   Post-est OLS"
		di as text _col(`varwidthm7') "         {c |}" _c
		di as text "   (alpha=" _c
		di as text %4.3f `e(alpha)' _c
		di as text ")"
	}
	else if "`e(method)'"=="lasso" {
		di as text _col(`varwidthm7') "Selected {c |}           Lasso   Post-est OLS"
	}
	else {
		di as err "internal DisplayCoefs error. unknown method."
		exit 1
	}
	di as text "{hline `varwidth1'}{c +}{hline 32}"
	local anynotpen = 0
	local i 1
	local lastcol = `col_ct' - `partial_ct'
	tokenize `vlist'								//  put elements of coef vector into macros 1, 2, ...
	while `i' <= `lastcol' {
		local vn ``i''
		fvstrip `vn'								// get rid of o/b/n prefix for display purposes
		local vn		`r(varlist)'
		_ms_display, element(`i') matrix(`beta') width(`varwidth') `baselevels'
		// in selected or notpen list?
		local isselnotpen	: list posof "`vn'" in selected0
		local isnotpen		: list posof "`vn'" in notpen
		local anynotpen		= `anynotpen' + `isnotpen'
		// note attached? base, empty, omitted
		qui _ms_display, element(`i') matrix(`beta')
		local note `r(note)'
		qui _ms_display, element(`i') matrix(`betaOLS')
		local noteOLS `r(note)'
		// if notpen, add footnote
		if `isnotpen' & "`note'"=="" {
			di as text "{helpb rlasso##notpen:*}" _c
		}
		if `isselnotpen' {
			// lasso coef
			if "`note'"=="" {
				di _col(`varwidth4') as res %15.7f el(`beta',1,`i') _c
			}
			else {
				di _col(`varwidth4') as text %15s "`note'" _c
			}
			// post-lasso coef - can be omitted if collinear
			if "`noteOLS'"=="" {
				di as res %15.7f el(`betaOLS',1,`i')
			}
			else {
				di as text %15s "`noteOLS'"
			}
		}
		else if "`note'"=="(omitted)" {
			// not selected
			di _col(`varwidth4') as text %15s "(not selected)" _c
			di                   as text %15s "(not selected)"
		}
		else {
			// other eg base var
			di as text %15s "`note'" _c
			di as text %15s "`noteOLS'"
		}
		local ++i
	}
	if `partial_ct' {
		di as text "{hline `varwidth1'}{c +}{hline 32}"
		di as text _col(`varwidthm13') "Partialled-out{help lasso2##examples_partialling:*}{c |}"
		di as text "{hline `varwidth1'}{c +}{hline 32}"
		local i = `lastcol'+1
		while `i' <= `col_ct' {
			local vn ``i''
			fvstrip `vn'								// get rid of o/b/n prefix for display purposes
			local vn		`r(varlist)'
			_ms_display, element(`i') matrix(`beta') width(`varwidth') `baselevels'
			// note attached? base, empty, omitted
			qui _ms_display, element(`i') matrix(`beta')
			local note `r(note)'
			qui _ms_display, element(`i') matrix(`betaOLS')
			local noteOLS `r(note)'
			// lasso coef
			if "`note'"=="" {
				di _col(`varwidth4') as res %15.7f el(`beta',1,`i') _c
			}
			else {
				di _col(`varwidth4') as text %15s "`note'" _c
			}
			// post-lasso coef - can be omitted if collinear
			if "`noteOLS'"=="" {
				di as res %15.7f el(`betaOLS',1,`i')
			}
			else {
				di as text %15s "`noteOLS'"
			}
			local ++i
		}
	}
	di as text "{hline `varwidth1'}{c BT}{hline 32}"
	
	if `anynotpen' {
		di "{help lasso2##examples_partialling:*Not penalized}"
	}
	
end

// subroutine to process user-supplied lambda(s)
prog define getlambdamat, rclass
	version 13
	syntax ,					///
	[							///
	lscalar(string)				///
	lmatrix(string)				///
	lfactor(real 1)				///
	]

	tempname lambdamat
	if "`lscalar'"~="" {
		// provided as scalars
		foreach lambda_i of local lscalar {
			mat `lambdamat'	= nullmat(`lambdamat'), `lambda_i'
		}
	}
	else if "`lmatrix'"~="" {
		// provided as matrix
		mat `lambdamat'		= `lmatrix'
	}
	else {
		di as err "internal lasso2 error - parsing lambda option"
		exit 198
	}
	// optional adjustment using undocumented lfactor option
	// used for CV
	mat `lambdamat'			= `lambdamat' * `lfactor'
	return matrix lambdamat	= `lambdamat'
	
end


// internal version of fvstrip 1.01 ms 24march2015
// takes varlist with possible FVs and strips out b/n/o notation
// returns results in r(varnames)
// optionally also omits omittable FVs
// expand calls fvexpand either on full varlist
// or (with onebyone option) on elements of varlist

program define fvstrip, rclass
	version 11.2
	syntax [anything] [if] , [ dropomit expand onebyone NOIsily ]
	if "`expand'"~="" {												//  force call to fvexpand
		if "`onebyone'"=="" {
			fvexpand `anything' `if'								//  single call to fvexpand
			local anything `r(varlist)'
		}
		else {
			foreach vn of local anything {
				fvexpand `vn' `if'									//  call fvexpand on items one-by-one
				local newlist	`newlist' `r(varlist)'
			}
			local anything	: list clean newlist
		}
	}
	foreach vn of local anything {									//  loop through varnames
		if "`dropomit'"~="" {										//  check & include only if
			_ms_parse_parts `vn'									//  not omitted (b. or o.)
			if ~`r(omit)' {
				local unstripped	`unstripped' `vn'				//  add to list only if not omitted
			}
		}
		else {														//  add varname to list even if
			local unstripped		`unstripped' `vn'				//  could be omitted (b. or o.)
		}
	}
// Now create list with b/n/o stripped out
	foreach vn of local unstripped {
		local svn ""											//  initialize
		_ms_parse_parts `vn'
		if "`r(type)'"=="variable" & "`r(op)'"=="" {			//  simplest case - no change
			local svn	`vn'
		}
		else if "`r(type)'"=="variable" & "`r(op)'"=="o" {		//  next simplest case - o.varname => varname
			local svn	`r(name)'
		}
		else if "`r(type)'"=="variable" {						//  has other operators so strip o but leave .
			local op	`r(op)'
			local op	: subinstr local op "o" "", all
			local svn	`op'.`r(name)'
		}
		else if "`r(type)'"=="factor" {							//  simple factor variable
			local op	`r(op)'
			local op	: subinstr local op "b" "", all
			local op	: subinstr local op "n" "", all
			local op	: subinstr local op "o" "", all
			local svn	`op'.`r(name)'							//  operator + . + varname
		}
		else if"`r(type)'"=="interaction" {						//  multiple variables
			forvalues i=1/`r(k_names)' {
				local op	`r(op`i')'
				local op	: subinstr local op "b" "", all
				local op	: subinstr local op "n" "", all
				local op	: subinstr local op "o" "", all
				local opv	`op'.`r(name`i')'					//  operator + . + varname
				if `i'==1 {
					local svn	`opv'
				}
				else {
					local svn	`svn'#`opv'
				}
			}
		}
		else if "`r(type)'"=="product" {
			di as err "fvstrip error - type=product for `vn'"
			exit 198
		}
		else if "`r(type)'"=="error" {
			di as err "fvstrip error - type=error for `vn'"
			exit 198
		}
		else {
			di as err "fvstrip error - unknown type for `vn'"
			exit 198
		}
		local stripped `stripped' `svn'
	}
	local stripped	: list retokenize stripped						//  clean any extra spaces
	
	if "`noisily'"~="" {											//  for debugging etc.
di as result "`stripped'"
	}

	return local varlist	`stripped'								//  return results in r(varlist)
end


// Internal version of matchnames
// Sample syntax:
// matchnames "`varlist'" "`list1'" "`list2'"
// takes list in `varlist', looks up in `list1', returns entries in `list2', called r(names)
program define matchnames, rclass
	version 11.2
	args	varnames namelist1 namelist2

	local k1 : word count `namelist1'
	local k2 : word count `namelist2'

	if `k1' ~= `k2' {
		di as err "namelist error"
		exit 198
	}
	foreach vn in `varnames' {
		local i : list posof `"`vn'"' in namelist1
		if `i' > 0 {
			local newname : word `i' of `namelist2'
		}
		else {
* Keep old name if not found in list
			local newname "`vn'"
		}
		local names "`names' `newname'"
	}
	local names	: list clean names
	return local names "`names'"
end


version 13
mata:
void s_maketemps(real scalar p)
{
	(void) st_addvar("double", names=st_tempname(p), 1)
	st_global("r(varlist)",invtokens(names))
}


// END MATA SECTION
end
