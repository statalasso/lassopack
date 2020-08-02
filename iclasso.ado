*! iclasso 1.0.02 31july2020
*! based on lasso2 1.0.10 14oct2019
*! lassopack package 1.3.1
*! authors aa/ms

* Notes:
* Add display of selected lambda to rlasso.
* Make display consistent with rlasso (cons).
* Sort out lambda vs newlambda vs lambdamat etc.
* fix unused options in iclasso wrapper.
* bug in lasso2 - stdcoef w/o prestd is ignored rather than triggering stdcoef
* allow stdcoef w/o prestd (use std on-the-fly for speed)
* now stdcoef => prestd automatically
* add stdcoef to rlasso
* ebicgamma
* active set convergence issue (p. 7 FHT JSS 2010)
* just warn about lambda edges, or more active response?
* in lasso2, emphasise that adaloadings should be for standardized coeffs
* state coefs are standarized in output?
* to fix - prestd + stdcoef + cons and any partialled out wrong - cons and partialled-out vars need rescaling
* in fact ... should just perhaps not report them. stdcoef = norecover.
* norecover + stdcoef applies to lasso2, cvlasso etc. as well.
* but has implications for partialling-out...

program iclasso, eclass sortpreserve
	version 13
	syntax [anything] [if] [in] [,					///
			displayall 								/// display zeros of beta, only applicable for one lambda
			NOGrid									/// no display of grid
			long 									///
			VERsion									///
			IC(string)								///
			* 										///
			]
	
	local lversion 1.0.00
	local pversion 1.3.1

	if "`version'" != "" {							//  Report program version number, then exit.
		di in gr "iclasso version `lversion'"
		di in gr "lassopack package version `pversion'"
		ereturn clear
		ereturn local version		`lversion'
		ereturn local pkgversion	`pversion'
		exit
	}
	
	if ~replay() {
		// if not replay, estimate / grid search
		_iclasso `anything' `if' `in', ic(`ic') `options'
	}
	else if "`ic'"~="" & "`ic'"~="`e(ic)'" {
		// replay with ic(.) option
		// if ic is different, need to reestimate with lasso3
		// new estimation results will be in e(.) macros
		_reestimate, ic(`ic')
	}

	// display results
	if "`nogrid'"=="" {
		DispGrid
	}
	DispHeader
	DisplayCoefs, `displayall'

	*
 
end

program _iclasso, eclass sortpreserve
	version 13

	syntax varlist(numeric min=2 fv ts) [if] [in] [,	///
			NOTPen(string) 							/// list of variables not penalised
			PARtial(string)							/// string so that list can contain "_cons"
			fe										/// do within-transformation
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
			/// info criterion
			IC(string)								///
			///
			/// lambda
			Lambda(numlist >0 min=1 descending)		/// expect list (syntax check captured below)
			LCount0(integer 100)					/// 0 to preserve original
			LMAX(real 0) 							///
			LMINRatio(real 1e-4)					/// ratio of maximum to minimum lambda
			LFactor(real 1) 						/// 
			LAMBDAMat0(string)						/// alternative: specify lambda as matrix; 0 to preserve original
			Ploadings(string) 						///
			UNITLoadings							///
			///
			/// standardization
			PREStd 									///
			STDCoef 								/// 
			/// choice of estimator
			SQRT 									/// square-root lasso
			RIDGE									///
			LASSO									///
			ENET									///
			OLS										///
			///
			/// adaptive lasso
			ADAptive  								/// adaptive lasso
			ADATheta(numlist >=0 ascending) 		/// theta parameter for adaLASSO
			ADALoadings(string)						///
			///
			///
			ALPha(numlist >=0 ascending) 			/// elastic net parameter
			AMIN(real 0)							///
			AMAX(real 1)							///
			ACount0(integer 5)						/// 0 to preserve original
													///
			POSTAll									///
			holdout(varlist numeric min=1 max=1) 	///
													///
			NOFTOOLS								///
			psolver(string)							/// optional choice of solver
													///
			*										///
			]


	** flags
	local ridgeflag		="`ridge'"~=""
	local sqrtflag		="`sqrt'"~=""
	local enetflag		="`enet'"~=""
	local lassoflag		="`lasso'"~="" | (`ridgeflag' + `sqrtflag' + `enetflag')==0
	local feflag=("`fe'"~="")
	*

	** set IC default
	if "`ic'"=="" {
		local ic bic
	}
	else {
		local ic = strlower("`ic'")
	}
	*
	
	** alpha list / elastic net
	// elastic net
	// alpha list implies elastic net
	// alphalist0 is original list in original order
	// alphalist is reordered for estimation purposes
	local acount			: word count `alpha'
	if `acount' > 0 {
		// user provided a list of alphas
		local enetflag		= 1
		local lassoflag		= 0
		// overwrite default acount0
		local acount0		= `acount'
		local alphalist0	`alpha'
	}
	else if `enetflag' & `acount'==0 {
		// user indicated elastic net but with no list of alphas
		// calculate default alpha grid using default acount0
		local adelta		= (`amax'-`amin')/(`acount0'-1)
		numlist "`amin'(`adelta')`amax'"
		local alphalist0	`r(numlist)'
		// update acount
		local acount		: word count `alphalist0'
	}
	else {
		// overwrite default acount0
		local acount0		= `acount'
	}
	// process alphas, update alpha count, etc.
	// alpha defined for all estimators so acount must be updated
	if (`acount'==0) {
		if (`lassoflag' | `sqrtflag') {
			local alpha		= 1
			local alphalist	1
		}
		else if `ridgeflag' {
			local alpha		= 0
			local alphalist	0
		}
		// and update acount
		local acount	= 1
	}
	else {
		// the lambda grid created by lassoutils depends on alpha
		// we want the grid created in the first iteration to have the smallest lambda endpoint
		// ensured by making the first alpha the biggest
		forvalues i=1/`acount' {
			local a	: word `i' of `alphalist0'
			local alphalist `a' `alphalist'
		}
		// overwrite alpha
		// first alpha to iterate over
		local alpha 	: word 1 of `alphalist'
	}
	*

	** theta list / adalasso
	local tcount0			: word count `adatheta'		// count of adaptive lasso parameter; 0 for consistency with lambda
	local adaflag			="`adaptive'"~="" | `tcount0'>0 | "`adaloadings'"~=""
	if `adaflag' {
		local adaptive adaptive
	}
	// maintain thetalist0/tcount0 and thetalist0/thetalist; 0 fixed, no zero may change
	if (`tcount0'==0) {
		// default list is [ 1 ]
		local thetalist0	1
		local thetalist		`thetalist0'
		local tcount		= 1
		// update
		local tcount0		= 1
		local tcount		= `tcount0'
		// first theta to iterate over
		local theta			= 1
	}
	else {
		local tcount		= `tcount0'
		local thetalist0	`adatheta'
		local thetalist		`thetalist0'
		local theta	: word 1 of `thetalist0'
	}
	// any theta not equal to 1 implied prestd
	local one				1
	local thetanotone		: list thetalist - one
	local thetanotone		: word count `thetanotone'
	if `thetanotone' {
		di as text "note: theta <1 or >1 implies prestd; data will be pre-standardized" 
		local prestd		prestd
	}
	*

	****** syntax checks *******************************************************
	if (`ridgeflag'+`lassoflag'+`enetflag'+`sqrtflag')>1 {
		di as err "incompatible options: `ridge' `lasso' `enet' `sqrt'"
		exit 198
	}
	local iclist aic aicc ebic bic
	local checklist	: list ic - iclist
	local checknum	: word count `checklist'
	if `checknum' {
		di as err "syntax error - `ic' invalid IC option"
		exit 198
	}
	if `acount'==1 {
		if (`alpha'>1) | (`alpha'<0) {
			di as err "error - alpha is out of range."
			exit 198
		}
		if (`lassoflag' | `sqrtflag') & (`alpha'~=1)  {
			di as err "error - `lasso'`sqrt' requires alpha=1"
			exit 198
		}
		if (`ridgeflag') & (`alpha'~=0)  {
			di as err "error - ridge requires alpha=0"
			exit 198
		}
		if ("`sqrt'"!="") & (`alpha'!=1) {
			di as error "error - sqrt-lasso allowed only with alpha=1."
			exit 198
		}
	}
	else if `enetflag'==0 {
		di as err "error - alpha list allowed only with elastic net estimator"
		exit 198
	}
	else if (`amin'<0) | (`amax'>1) {
		di as err "error - alpha is out of range."
		exit 198
	}
	local notpenpar : list notpen & partial
	if ("`notpenpar'"!="") {
		di as error "`notpenpar' listed in both notpen(.) and partial(.)"
		exit 198
	}
// ms bug fix: 
*	if ("`stdcoef'"!="") & ("`prestd'"!="") {
	if ("`stdcoef'"!="") & ("`prestd'"=="") {
		di as text "note: option stdcoef implies prestd; data will be pre-standardized" 
		local prestd prestd
	}
	local checkflag 	= ("`ploadings'"!="")+("`unitloadings'"!="")+("`adaptive'"!="")
	if `checkflag'>1 {
		di as error "error: cannot combine options ploadings(.), unitloadings and/or adaptive"
		exit 198
	}
	*
	****************************************************************************
	
	*** debug mode; create flag
	local debugflag	=("`debug'"~="")
	*
	
	*** Record which observations have non-missing values
	marksample touse
	// need to check panel var as well
	if `feflag' {
		cap xtset
		local ivar	`r(panelvar)'
	}
	markout `touse' `varlist' `ivar' `holdout'
	sum `touse' if `touse', meanonly		//  will sum weight var when weights are used
	local N		= r(N)
	tempvar toest
	qui gen `toest' = `touse'
	if ("`holdout'"!="") {
		assert `holdout' == 1 | `holdout'==0 if `touse'
		qui replace `toest' = 0 if `holdout'
	}
	*

	*** FEs. Create 1/0 flag.
	if `feflag' {
		if "`ivar'"=="" {
			di as err "Error: fe option requires data to be xtset"
			exit 459
		}
		// fe transformation may expect data to be sorted on ivar
		local sortvar	: sortedby
		local sortvar	: word 1 of `sortvar'				// in case sorted on multiple variables
		if "`ivar'"~="`sortvar'" {
			di as text "(sorting by xtset panelvar `ivar')"
			sort `ivar'
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
	local stdcoefflag	= ("`stdcoef'"~="")		// return coef estimates in std units
	local prestdflag 	= ("`prestd'"~="")		// =1 if data to be pre-standardized
	// default is to use standardization loadings; overridden by ploadings, unitloadings, pre-standardization, adaptive
	local stdloadflag	= ("`ploadings'`unitloadings'`prestd'`adaptive'"=="")
	local sqrtflag 		= ("`sqrt'"!="")
	local parrflag		= ("`norecover'"=="")	
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

	//  p is number of penalized vars in the model; follows convention in BCH papers
	//  p is calculated in lassoutils/_rlasso as number of model vars excluding constant
	//  here we calculate which of the model vars are unpenalized or omitted/base vars
	//  to provide as `pminus' to lassoutils/_rlasso (unless provided by user)
	//  do here so that code above is compatible with iclasso
	//  use _o names / display names since they have info on whether var is omitted/base/etc.
	if ~`pminus' {
		foreach vn of local varXmodel_d {								//  display names
			_ms_parse_parts `vn'
			// increment pminus if model variable is MISSING
			if r(omit) {
				local ++pminus
			}
		}
		foreach vn of local notpen_d {									//  display names
			_ms_parse_parts `vn'
			// increment pminus if notpen variable is NOT MISSING
			if ~r(omit) {
				local ++pminus
			}
		}
	}
	//  p0 here is total number of variables provided to model EXCLUDING constant
	local p0	: word count `varXmodel_o'
	local p		=`p0'-`pminus'
	// warn
	if `p'<=0 {
		di as text "warning: no penalized regressors; results are OLS"
	}
	//  now for error-checking below, p0 should INCLUDE constant unless partialled-out etc.
	local p0	=`p0'+`consflag'
	*

	******************* FE, partialling out, standardization ************************************
	//  If FE:    partial-out FEs from temp variables, then preserve,
	//            then partial-out low-dim ctrls from temp variables
	//            restore will restore all temp vars with only FEs partialled-out
	//  If no FE: leave original variables unchanged.
	//            partial-out low-dim ctrls from temp variables.
	//            if no FE/low-dim ctrls, no transform needed

	local dofminus	=0										//  overwritten by FE count
	local sdofminus	=`consmodel'							//  initial "small" df count is cons/nocons
	local dmflag	=0										//  initialize demeaned flag
	if `feflag' {											//  FE-transform all variables
		fvrevar `varY_o' `varX_o' if `touse'				//  in case any FV or TS vars in _o list
		local vlist `r(varlist)'
		lassoutils `vlist',								/// call on _o list
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
						psolver(`psolver')					/// optional choice of solver
						tvarlist(`varY_t' `varXmodel_t')	/// overwrite/initialize these
						partialflag(`partialflag')			/// triggers branching to partial utility
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
// ms note - now the original lcount is lcount0 and original matrix is lambdamat0
// macro lcount gets overwritten
	if ("`lambda'"!="") {
		// lambda list provided in macro form
		// overwrite default lcount
		local lcount0		: word count `lambda'
		local lcount		= `lcount0'
		if `lcount0'==1 {
			di as err "error - lambda list must have at least 2 entries"
			exit 198
		}
		// overwrite default lambda0
		tempname lambdamat0 lambdamat
		mat `lambdamat0'	= J(1,`lcount0',.)
		local j = 1
		foreach lami of local lambda {
			mat `lambdamat0'[1,`j'] = `lami'  
			local j=`j'+1
		}
		// optional adjustment using undocumented lfactor option
		mat `lambdamat0'	=`lambdamat0'*`lfactor'
		mat `lambdamat'		= `lambdamat0'
	}
	else if ("`lambdamat0'"!="") {
		// lambdamat0 provided
		// overwrite default lcount
		local lcount0	= colsof(`lambdamat0')
		local lcount	= `lcount0'
		if `lcount0'==1 {
			di as err "error - lambda list must have at least 2 entries"
			exit 198
		}
		// optional adjustment using undocumented lfactor option
		mat `lambdamat0'=`lambdamat0'*`lfactor'
		tempname lambdamat
		mat `lambdamat'=`lambdamat0'
	}
	else {
		// macro `lambdamat0' is empty so will be constructed by lassoutils
		local lcount	= `lcount0'
	}
	*

 	*************** alpha to matrix **************************************************
	if (`acount0'>1) {
		tempname alphamat0
		mat `alphamat0'	= J(1,`acount0',.)
		local j = 1
		foreach ai of local alphalist0 {
			mat `alphamat0'[1,`j'] = `ai'  
			local j=`j'+1
		}
	}
	*

 	*************** theta to matrix **************************************************
	if (`tcount0'>1) {
		tempname thetamat0
		mat `thetamat0'	= J(1,`tcount0',.)
		local j = 1
		foreach ti of local thetalist0 {
			mat `thetamat0'[1,`j'] = `ti'  
			local j=`j'+1
		}
	}
	*
	
	************* Partialling/standardization END ***********************************************

	************** Lasso estimation with transformed/partialled-out vars
	if "`verbose'`vverbose'"=="" {
		local quietly "quietly"							//  don't show lassoutils output
	}	

	************* Lasso/Ridge/Enet estimation - first iteration *********************************
	// will be the only iteration if looping through lambda only
	// will be initial iteration if then looping through alpha/theta/lambda
	
	if `adaflag' {
		`quietly' di as text "Initial iteration: theta=" %5.3f `theta' " ..."
	}
	`quietly' di as text "Initial iteration: alpha=" %5.3f `alpha' " ..."
	`quietly' lassoutils `varY_t',							///
					path									/// branches to _lassopath
					toest(`toest')							///
					xnames_o(`varXmodel_d')					/// display name
					xnames_t(`varXmodel_t')					///
					consflag(`consflag')					/// =0 if cons already partialled out or if no cons
					dmflag(`dmflag')						/// =1 if data have been demeaned
					dofminus(`dofminus')					/// dofs lost from FEs
					sdofminus(`sdofminus')					/// dofs lost from partialling
					notpen_o(`notpen_d') 					/// not penalised (display name)
					notpen_t(`notpen_t')					///
					lambda(`lambdamat0')					/// original (untruncated) lambdamat
					`adaptive'								///
					adatheta(`theta')						///
					adaloadings(`adaloadings')				///
					`sqrt'									///
					`ols'									///
					alpha(`alpha')							/// initial alpha
					stdy(`prestdY')							///
					stdx(`prestdX')							///
					stdl(`stdloadflag')						/// use standardisation loadings
					stdcoef(`stdcoefflag')					/// return in standard units
					ploadings(`ploadings') 					///
					`verbose' `vverbose'					///
					holdout(`holdout')						///
					lcount(`lcount0')						/// first iteration only
					lmax(`lmax')							/// first iteration only
					lminratio(`lminratio')					/// first iteration only
					`options'

	tempname lambdas lambdas0
	// truncated and untruncated lambda list
	mat `lambdas'	=r(lambdalist)
	mat `lambdas0'	=r(lambdalist0)
	// lambda grid may be shortened if endpoint encountered
	local lcount	=r(lcount)
	// min and max of untruncated grid
	local lmin0		=r(lmin0)
	local lmax0		=r(lmax0)

	tempname fullmat
	local cnames	lambda alpha theta aic aicc bic ebic s df rss ess r2
	mat `fullmat' = r(lambdalist), J(`lcount',1,`alpha'), J(`lcount',1,`theta')
	mat `fullmat'	= `fullmat', r(aic), r(aicc), r(bic), r(ebic), r(shat), r(dof), r(rss), r(ess), r(rsq)
	mat colnames `fullmat' = `cnames'
	// initialize rcount
	local rcount 1
	forvalues i=1/`lcount' {
		local rnames	`rnames' r`rcount'
		local ++rcount
	}
	// need this only once
	local ebicgamma = r(ebicgamma)
	// initialize mat lists and minimizing values
	tempname aicmat aiccmat bicmat ebicmat
	foreach m in aic aicc bic ebic {
		mat ``m'mat'			= J(1,7,.)
		mat colnames ``m'mat'	= icmin lambda alpha theta min max endpoint
		mat ``m'mat'[1,1]		= r(`m'min)
		mat ``m'mat'[1,2]		= `lambdas'[r(l`m'id),1]
		mat ``m'mat'[1,3]		= `alpha'
		mat ``m'mat'[1,4]		= `theta'
		mat ``m'mat'[1,5]		= `lmin0'
		mat ``m'mat'[1,6]		= `lmax0'
		// endpoint refers to truncated list
		// -1 <=> lower endpoint, +1 <=> upper endpoint, 0 <=> interior of lambda grid
		if ``m'mat'[1,2] == `lambdas'[1,1] {
			mat ``m'mat'[1,7]	= +1				//  upper endpoint (max)
		}
		else if ``m'mat'[1,2] == `lambdas'[`lcount',1] {
			mat ``m'mat'[1,7]	= -1				//  lower endpoint (min)
		}
		else {
			mat ``m'mat'[1,7]	= 0					//  interior of grid
		}
		
	}

	************* Elastic net / Adaptive LASSO estimation - 2nd+ iterations *****************

	// outer loop ii over theta
	// inner loop jj over alpha
	// ii=jj=1 initialization has already taken place so do not execute block if ii=jj=1
	// note we can now use the existing (untruncated) list of lambdas

	forvalues ii=1/`tcount' {
		local theta	: word `ii' of `thetalist'
		
		forvalues jj=1/`acount' {
			local alpha	: word `jj' of `alphalist'
		
			if ~((`ii'==1) & (`jj'==1)) {		// do not enter if ii=jj=1
						
				if `adaflag' {
					`quietly' di as text "Iterating: theta=" %5.3f `theta' " ..."
				}
				if `enetflag' {
					`quietly' di as text "Iterating: alpha=" %5.3f `alpha' " ..."
				}
			
				`quietly' lassoutils `varY_t',							///
								path									/// branches to _lassopath
								toest(`toest')							///
								xnames_o(`varXmodel_d')					/// display name
								xnames_t(`varXmodel_t')					///
								consflag(`consflag')					/// =0 if cons already partialled out or if no cons
								dmflag(`dmflag')						/// =1 if data have been demeaned
								dofminus(`dofminus')					/// dofs lost from FEs
								sdofminus(`sdofminus')					/// dofs lost from partialling
								notpen_o(`notpen_d') 					/// not penalised (display name)
								notpen_t(`notpen_t')					///
								lambda(`lambdamat0')					/// original (untruncated) lambdamat
								`adaptive'								///
								adatheta(`theta')						///
								adaloadings(`adaloadings')				///
								`sqrt'									///
								`ols'									///
								alpha(`alpha')							///
								stdy(`prestdY')							///
								stdx(`prestdX')							///
								stdl(`stdloadflag')						/// use standardisation loadings
								stdcoef(`stdcoefflag')					/// return in standard units
								ploadings(`ploadings') 					///
								`verbose' `vverbose'					///
								holdout(`holdout')						///
								lcount(`lcount0')						/// 
								lmax(`lmax')							/// 
								lminratio(`lminratio')					/// 
								`options'
		
				mat `lambdas'	=r(lambdalist)
				mat `lambdas0'	=r(lambdalist0)
				// lambda grid may be shortened if endpoint encountered
				local lcount	=r(lcount)
				// min and max of untruncated grid
				local lmin0		=r(lmin0)
				local lmax0		=r(lmax0)
				// append results
				tempname thismat
				mat `thismat' = r(lambdalist)
				mat `thismat'	= `thismat', J(`lcount',1,`alpha'), J(`lcount',1,`theta')
				mat `thismat'	= `thismat', r(aic), r(aicc), r(bic), r(ebic), r(shat), r(dof), r(rss), r(ess), r(rsq)
				mat `fullmat'	= `fullmat' \ `thismat'
				forvalues i=1/`lcount' {
					local rnames	`rnames' r`rcount'
					local ++rcount
				}
				
				foreach m in aic aicc bic ebic {
					if r(`m'min) < ``m'mat'[1,1] {
						mat ``m'mat'[1,1]	= r(`m'min)
						mat ``m'mat'[1,2]	= `lambdas'[r(l`m'id),1]
						mat ``m'mat'[1,3]	= `alpha'
						mat ``m'mat'[1,4]	= `theta'
						mat ``m'mat'[1,5]	= `lmin0'
						mat ``m'mat'[1,6]	= `lmax0'
						if ``m'mat'[1,2] == `lambdas'[1,1] {
							mat ``m'mat'[1,7]	= +1				//  upper endpoint (max)
						}
						else if ``m'mat'[1,2] == `lambdas'[`lcount',1] {
							mat ``m'mat'[1,7]	= -1				//  lower endpoint (min)
						}
						else {
							mat ``m'mat'[1,7]	= 0					//  interior of grid
						}
					}
				}
			}	// end block of code to execute in loop
	
		}	// end alpha loop

	}	// end theta loop

	mat rownames `fullmat' = `rnames'
	
	************* Finish up ********************************************************

	// set depending on which IC used
	local icmin		= ``ic'mat'[1,1]
	local llambda	= ``ic'mat'[1,2]
	local aalpha	= ``ic'mat'[1,3]
	local tmin		= ``ic'mat'[1,4]

	tempname lambdamat1
	mat `lambdamat1' = `llambda'

	*** Lasso re-estimation with selected hyperparameters
	`quietly' lassoutils `varY_t',							///
					path									/// branches to _lassopath
					toest(`toest')							///
					xnames_o(`varXmodel_d')					/// display name
					xnames_t(`varXmodel_t')					///
					consflag(`consflag')					/// =0 if cons already partialled out or if no cons
					dmflag(`dmflag')						/// =1 if data have been demeaned
					dofminus(`dofminus')					/// dofs lost from FEs
					sdofminus(`sdofminus')					/// dofs lost from partialling
					notpen_o(`notpen_d') 					/// not penalised (display name)
					notpen_t(`notpen_t')					///
					lambda(`lambdamat1')					/// chosen lambda
					`adaptive'								///
					adatheta(`tmin')						/// chosen theta
					adaloadings(`adaloadings')				///
					`sqrt'									///
					`ols'									///
					alpha(`aalpha')							/// chosen alpha
					stdy(`prestdY')							///
					stdx(`prestdX')							///
					stdl(`stdloadflag')						/// use standardisation loadings
					stdcoef(`stdcoefflag')					/// return in standard units
					ploadings(`ploadings') 					///
					`verbose' `vverbose'					///
					holdout(`holdout')						///
					`options'

		*** e-return lasso estimation results
		tempname b beta betaOLS Psi stdvec
		tempname betaAll betaAllOLS
		tempname rmse rmseOLS objfn r2 df
		if "`cluster'" ~= "" {
			local N_clust		=r(N_clust)
		}
		mat `beta'			=r(beta)		//  may be empty!
		mat `betaOLS'		=r(betaOLS)		//  may be empty!
		mat `betaAll'		=r(betaAll)
		mat `betaAllOLS'	=r(betaAllOLS)
		mat `Psi'			=r(Psi)
		//*//mat `sPsi'			=r(sPsi)
		mat `stdvec'		=r(stdvec)
		//*//scalar `lambda'		=r(lambda)
		//*//scalar `slambda'	=r(slambda)
		//*//scalar `lambda0'	=r(lambda0)
		scalar `rmse'		=r(rmse)		//  Lasso RMSE
		scalar `rmseOLS'	=r(rmseOLS)		//  post-Lasso RMSE
		scalar `r2'			=r(r2)
		scalar `df'			=r(df)
		scalar `objfn'		=r(objfn)
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
		ereturn local predict		iclasso_p
		ereturn local cmd			iclasso
 		//*//ereturn scalar pminus		=`pminus'
		//*//ereturn scalar center		=`center'
		ereturn scalar cons			=`consmodel'
		//*//ereturn scalar slambda		=`slambda'
		//*//ereturn scalar lambda0		=`lambda0'

		if "`N_clust'" ~= "" {
			ereturn local clustvar	`clustvar'
			ereturn scalar N_clust	=`N_clust'
		}
		if "`N_g'" ~= "" {
			ereturn scalar N_g		=`N_g'
		}
		ereturn scalar fe			=`feflag'
		ereturn scalar rmse			=`rmse'
		ereturn scalar rmseOLS		=`rmseOLS'
		ereturn scalar r2			=`r2'
		ereturn scalar df			=`df'
		if "`method'"=="sqrt-lasso" {
			ereturn scalar prmse	=`objfn'
		}
		else {
			ereturn scalar pmse		=`objfn'
		}
		ereturn scalar p			=`p'
		ereturn scalar k			=`k'				//  number of all estimated coefs INCLUDING PARTIALLED-OUT AND CONSTANT
		ereturn scalar s			=`s'				//  number of selected

		ereturn matrix stdvec		=`stdvec'
		//*//ereturn matrix sPsi 		=`sPsi'
		ereturn matrix Psi 			=`Psi'
		ereturn matrix betaAllOLS	=`betaAllOLS'
		ereturn matrix betaAll		=`betaAll'
		ereturn matrix betaOLS		=`betaOLS'
		ereturn matrix beta			=`beta'

		// iclasso-specific:
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
		
		*** more iclasso ereturns
		ereturn scalar lambda		=`lambdamat1'[1,1]	//  IC-maximizing lambda
		ereturn scalar alpha		=`aalpha'
		ereturn scalar theta		=`tmin'
		ereturn scalar fe 			=`feflag'
		ereturn scalar sqrt  		=`sqrtflag'
		ereturn scalar prestd		=`prestdflag'
		ereturn scalar ols 			=`olsflag'
		ereturn scalar adaptive		= "`adaptive'"!=""
		if `adaflag' {
			ereturn scalar tmin		=`tmin'
		}
		if "`adaloadings'"~="" {
			ereturn matrix adaloadings	=`adaloadings', copy	//  copy so that original stays in memory
		}
		ereturn local ic			`ic'
		ereturn scalar lasso		=`lassoflag'
		ereturn scalar ridge		=`ridgeflag'
		ereturn scalar enet			=`enetflag'
		ereturn scalar dofminus		=`dofminus'
		ereturn scalar sdofminus	=`sdofminus'

		*** iclasso grid ereturns
// refers to untruncated grid
		ereturn scalar lcount 		=`lcount0'
// ms note - lambdas0 is a col vector, lambdamat0 is a row vector, alphamat0 is a row vector ... we save col vectors
//		ereturn matrix lambdamat	=`lambdas0'
		if `acount0'>1 {
			mat `alphamat0'			= `alphamat0''
			ereturn matrix alphamat	= `alphamat0'
		}
		if `tcount0'>1 {
			mat `thetamat0'			= `thetamat0''
			ereturn matrix thetamat	= `thetamat0'
		}
//		ereturn scalar lmax			=`lmax'
//		ereturn scalar lmin			=`lmin'
		ereturn matrix aicmat		= `aicmat'
		ereturn matrix aiccmat		= `aiccmat'
		ereturn matrix bicmat		= `bicmat'
		ereturn matrix ebicmat		= `ebicmat'
		ereturn scalar ebicgamma	= `ebicgamma'
		ereturn matrix gridmat		= `fullmat'


end

program define _reestimate, eclass
	version 13
	syntax [anything] , IC(string)

		tempname icresults icmat
		tempname b beta betaOLS betaAll betaAllOLS
		tempname rmse prmse rmseOLS pmse r2 df
		tempvar esample
		local y						`e(depvar)'
		local X						`e(varX)'
		// icmat has maximizing hyperparameter values for specified IC criterion
		mat `icmat'					= e(`ic'mat)
		local lambda				= `icmat'[1,2]
		local alpha					= `icmat'[1,3]
		local theta					= `icmat'[1,4]
		if e(sqrt) {
			local l3options `l3options' sqrt
		}
		if e(adaptive) {
			local l3options `l3options' adaptive adatheta(`theta')
			tempname aload
			mat `aload'				= e(adaloadings)
			if ~matmissing(`aload') {
				local l3options `l3options' adaloadings(`aload')
			}
		}
		local l3options				`l3options' notpen(`e(notpen)')
		local consname				_cons
		local partial				`e(partial)'
		local partial				: list partial - consname
		local l3options				`l3options' partial(`partial')
		qui gen `esample'			= e(sample)
		_estimates hold `icresults'
		
		qui lasso3 `y' `X' if `esample', lambda(`lambda') alpha(`alpha') `l3options'
		
		mat `b'						= e(b)
		mat `beta'					= e(beta)
		mat `betaOLS'				= e(betaOLS)
		mat `betaAll'				= e(betaAll)
		mat `betaAllOLS'			= e(betaAllOLS)
		local k						= e(k)
		local selected				`e(selected)'
		local selected0				`e(selected0)'
		local s						= e(s)
		local s0					= e(s0)
		local niter					= e(niter)
		local method				`e(method)'
		scalar `rmse'				= e(rmse)
		scalar `rmseOLS'			= e(rmseOLS)
		scalar `pmse'				= e(pmse)
		scalar `prmse'				= e(prmse)
		scalar `r2'					= e(r2)
		scalar `df'					= e(df)
		_estimates unhold `icresults'
		ereturn repost b = `b', resize
		ereturn matrix betaAllOLS	=`betaAllOLS'
		ereturn matrix betaAll		=`betaAll'
		ereturn matrix betaOLS		=`betaOLS'
		ereturn matrix beta			=`beta'
		ereturn scalar k			=`k'
		ereturn scalar s			=`s'
		ereturn scalar s0			=`s0'
		ereturn scalar rmse			=`rmse'
		ereturn scalar rmseOLS		=`rmseOLS'
		ereturn scalar r2			=`r2'
		ereturn scalar df			=`df'
		if e(sqrt) {
			ereturn scalar prmse		=`prmse'
		}
		else {
			ereturn scalar pmse			=`pmse'
		}
		ereturn scalar niter		=`niter'
		ereturn scalar lambda		=`lambda'
		ereturn scalar alpha		=`alpha'
		if e(adaptive) {
			ereturn scalar theta	=`theta'
		}
		ereturn local selected		`selected'
		ereturn local selected0		`selected0'
		ereturn local ic			`ic'
		ereturn local method		`method'

end

program define DispGrid
	version 11.2

	tempname gridmat
	mat `gridmat'		= e(gridmat)
	local gridrows		= rowsof(`gridmat')
	
	tempname aicmat aiccmat bicmat ebicmat
	mat `aicmat'		= e(aicmat)
	mat `aiccmat'		= e(aiccmat)
	mat `bicmat'		= e(bicmat)
	mat `ebicmat'		= e(ebicmat)
	
	local alphachar = uchar(945)
	local thetachar = uchar(952)
	local bline {hline 9}
	local header {space 3}Lambda{space 1}
	if e(enet) {
		local bline `bline'{hline 5}
		local header `header'{space 2}`alphachar'{space 2}
	}
	if e(adaptive) {
		local bline `bline'{hline 5}
		local header `header'{space 2}`thetachar'{space 2}
	}
	local bline `bline'{hline 1}{c +}{hline 40}{c +}{hline 18}
	local header `header'{c |}{space 4}AIC{space 7}AICC{space 6}BIC{space 7}EBIC{space 2}{c |}{space 3}s{space 4}df{space 4}R-sq

	di
//	di as text "{hline 80}"
	di as text "`header'"
	di as text "`bline'"
	forvalues i=1/`gridrows' {
	
		if `i'==1 {
			// initialize
			if e(enet) {
				local this_alpha	=`gridmat'[`i',"alpha"]
			}
			if e(adaptive) {
				local this_theta	=`gridmat'[`i',"theta"]
			}
		}
		else {
			// print line if alpha or theta change (but not twice)
			local blineflag	= 0
			if e(enet) {
				if `this_alpha'~=`gridmat'[`i',"alpha"] {
					// display line and update alpha
					local this_alpha	=`gridmat'[`i',"alpha"]
					di as text "`bline'"
					local blineflag		= 1
				}
			}
			if e(adaptive) {
				if `this_theta'~=`gridmat'[`i',"theta"] {
					// display line and update theta
					local this_theta	=`gridmat'[`i',"theta"]
					if `blineflag'==0 {
						di as text "`bline'"
					}
					local blineflag		= 1
				}
			}
		}
	
		local c	=0
		di as res %9.4f `gridmat'[`i',"lambda"] _c
		local c	=`c'+11
		di _col(`c') _c
		if e(enet) {
			local c	=`c'+1
			di _col(`c') _c
			di as res %3.1f `gridmat'[`i',"alpha"] _c
			local c	=`c'+4
			di _col(`c') _c
		}
		if e(adaptive) {
			local c	=`c'+1
			di _col(`c') _c
			di as res %3.1f `gridmat'[`i',"theta"] _c
			local c	=`c'+4
			di _col(`c') _c
		}
		di as text _col(`c') "{c |}" _c
		local c	=`c'+2
		di _col(`c') _c
		
		di as res %8.2f `gridmat'[`i',"aic"] _c
		if `gridmat'[`i',"aic"]==`aicmat'[1,"icmin"] {
			di as text "*" _c
		}
		local c	=`c'+10
		di _col(`c') _c
		
		di as res %8.2f `gridmat'[`i',"aicc"] _c
		if `gridmat'[`i',"aicc"]==`aiccmat'[1,"icmin"] {
			di as text "*" _c
		}
		local c	=`c'+10
		di _col(`c') _c

		di as res %8.2f `gridmat'[`i',"bic"] _c
		if `gridmat'[`i',"bic"]==`bicmat'[1,"icmin"] {
			di as text "*" _c
		}
		local c	=`c'+10
		di _col(`c') _c

		di as res %8.2f `gridmat'[`i',"ebic"] _c
		if `gridmat'[`i',"ebic"]==`ebicmat'[1,"icmin"] {
			di as text "*" _c
		}
		local c	=`c'+9
		di _col(`c') _c

		di as text _col(`c') "{c |}" _c
		local c	=`c'+2
		di _col(`c') _c

		di as res %3.0f `gridmat'[`i',"s"] _c
		local c	=`c'+3
		di _col(`c') _c

		di as res %7.1f `gridmat'[`i',"df"] _c
		local c	=`c'+7
		di _col(`c') _c
		
		di as res %7.3f `gridmat'[`i',"r2"] _c
		di
	}
	di as text "{helpb iclasso##minic:*}indicates minimum IC."
end

program define DispHeader
	version 11.2

	tempname icmat
	local ic 			`e(ic)'
	mat `icmat'			= e(`ic'mat)
	local ebicgamma		= e(ebicgamma)
	local icmin			= `icmat'[1,1]
	local endpoint		= `icmat'[1,7]
	
	tempname lambdamat alphamat thetamat
	local lmin			= `icmat'[1,5]
	local lmax			= `icmat'[1,6]
	mat `alphamat'		= e(alphamat)
	mat `thetamat'		= e(thetamat)
	local acount		= rowsof(`alphamat')
	if `acount'>1 {
		mata: st_local("alphalist",invtokens(strofreal(st_matrix("e(alphamat)"))'))
		// insert commas
		local alphalist	: subinstr local alphalist " " ", ", all
	}
	local tcount		= rowsof(`thetamat')
	if `tcount'>1 {
		mata: st_local("thetalist",invtokens(strofreal(st_matrix("e(thetamat)"))'))
		// insert commas
		local thetalist	: subinstr local thetalist " " ", ", all
	}

	local icstrupper=strupper("`ic'")
	di as text ""
	di as text "Minimized `icstrupper':" _col(26) as res %9.2f `icmin'
	if "`ic'"=="ebic" {
		if `ebicgamma'==0 {
			local extratext " (equivalent to BIC)"
		}
		di as text "EBIC gamma:" _col(27) as res %5.3f `ebicgamma' as text "`extratext'"
	}
	if `e(adaptive)' {
		di as text "Adaptive theta:" _col(26) as res %5.1f `e(tmin)'
	}
	if `e(enet)' {
		if `e(alpha)'==0 {
			local extratext " (Ridge)"
		}
		else if `e(alpha)'==1 {
			local extratext " (Lasso)"
		}
		di as text "Elastic net with lambda = " as res %8.4f `e(lambda)' as text " and alpha = " as res %5.3f `e(alpha)' as text "`extratext'"
	}
	else {
		di as text "Lambda:" _col(29) as res %8.4f `e(lambda)'
	}

	di as text "Lambda grid:"  as res %4.0f `e(lcount)' as text " points in" _col(30) "[" as res %5.3f `lmin' as text ", " as res %7.3f `lmax' as text "]"
	
	if `endpoint' ~= 0 {
		if `endpoint' == -1 {
			local whichend "lower"
		}
		else {
			local whichend "upper"
		}
		di as text "Warning: lambda is at the `whichend' boundary of the evaluated gridpoints."
	}

	if `acount' > 1 {
		di "{txt}Alpha grid:"					 _col(30) "[{res}`alphalist'{txt}]"
	}
	if `tcount' > 1 {
		di "{txt}Theta grid:"					 _col(30) "[{res}`thetalist'{txt}]"
	}
	di "{txt}N:"								_col(25) as res %7.0f e(N)
	di "{txt}p (ex. partialled, cons):"			_col(28) as res %4.0f e(p)
	di "{txt}s (ex. partialled, cons):"			_col(28) as res %4.0f e(s)
	di "{txt}#partialled+cons:"					_col(28) as res %4.0f e(sdofminus)
	if e(alpha)<1 {
		di "{txt}Effective df (incl. p, c):"	_col(28) as res %9.4f e(df)
	}
	else {
		di "{txt}Effective df (incl. p, c):"	_col(28) as res %4.0f e(df)
	}
	di "{txt}In-sample R-sq: "					_col(28) as res %9.4f e(r2)
	
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

// Used in rlasso and lasso2.
// version  2020-02-21
prog DisplayCoefs

	syntax	,								///
		[									///
		displayall							///  full coef vector in display (default=selected only)
		varwidth(int 17)					///
		]
	
	local cons			=e(cons)

	if colsof(e(betaAll)) > e(p) {
		// betaAll has partialled-out coefs to display
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
		mat `beta'		=e(betaAll)
		mat `betaOLS'	=e(betaAllOLS)
		local col_ct	=colsof(`beta')
		local vlist		: colnames `beta'
		local vlistOLS	: colnames `betaOLS'
		local baselevels baselevels
	}
	else if e(k)>0 {							//  display only selected, but only if there are any
		mat `beta'		=e(beta)
		mat `betaOLS'	=e(betaOLS)
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
		di as text _col(`varwidthm13') "Partialled-out{help iclasso##examples_partialling:*}{c |}"
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
		di "{help iclasso##examples_partialling:*Not penalized}"
	}
	
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
