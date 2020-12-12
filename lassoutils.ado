*! lassoutils 1.2.04 8dec2020
*! lassopack package 1.4.1
*! authors aa/cbh/ms

* Adapted/expanded from lassoShooting version 12, cbh 25/09/2012
* mm_quantile from moremata package by Ben Jann:
*   version 1.0.8  20dec2007  Ben Jann

* Notes:
* Partialling out, temp vars, factor vars, FE transform all handled
* by calling programs.
* names_o and names_t (via options) are original and (possibly) temp names of Xs
* lassoutils branches to internal _rlasso, _lassopath, _fe, _unpartial, _partial, _std
* current cluster code is memory-intensive since it works with a full n x p matrix instead of cluster-by-cluster

* Updates (release date):
* 1.0.05  (30jan2018)
*         First public release.
*         Added seed(.) option to rlasso to control rnd # seed for xdep & sup-score.
*         Fixed up return code for lassoutils (method, alpha).
*         Promoted to required version 13 or higher.
*         Introduced centerpartial(.) Mata program for use with rlasso; returns X centered and with Xnp partialled-out.
*         Separate fields in datastruct: sdvec, sdvecpnp, sdvecp, sdvecnp. Latter two conformable with Xp and Xnp.
*         Added dots option for simulations (supscore, xdep).
*         Recoding relating to different treatment of cross-validation.
*         Changes to _std, _fe, _partial relating to holdout vs full sample.
*         Recoding of cons flag; now also dmflag to indcate zero-mean data.
* 1.0.06  (10feb2018)
*         Misc Mata speed tweaks following advice at http://scorreia.com/blog/2016/10/06/mata-tips.html,
*         e.g., evaluating for limits before loop; referring to vector elements with a single subscript.
*         Rewrite of FE transform to use Sergio Correia's FTOOLS package (if installed).
* 1.0.07  (17feb2018)
*         Bug fix related to ftools - leaves matastrict set to on after compilation, causing rest of lassoutils
*         to fail to load. Fix is to reset matastrict to what it was before calling ftools,compile.
* 1.0.08  (5apr2018)
*	      Changed the default maximum lambda for the elastic net: it was 2*max(abs(X'y)),
*		  and is now 2*max(abs(X'y))/max(0.001,alpha); see Friedman et al (2010, J of Stats Software). 
*		  Added Mata programs getInfoCriteria() and getMinIC(). getInfoCriteria() calculates information criteria
*		  along with RSS and R-squared. getMinIC() obtains minimum IC values and corresponding lambdas.
*		  Degrees of freedom calculation was added to DoLassoPath() and DoSqrtLassoPath().
* 		  Misc changes to the initial beta (= Ridge estimate), which in some cases didn't account for 
*		  penalty loadings.
* 		  Added lglmnet option to facilitate comparison with glmnet (was dysfunctional 'ladjust' option).
*         Restructured calc of lambda with cluster to match JBES paper; prev used sqrt(nclust) rather than sqrt(N),
*         now always 2c*sqrt(N)*invnormal(1-(gamma/log(N))/(2*p))). Fixed bug in calc of lambda for sqrt-lasso with cluster.
*         Definition of xdep lambda + cluster also changed slightly (by factor of sqrt(nclust/(nclust-1)).
*         Undocumented nclust1 option mimics behavior of CBH lassoCluster.ado; uses (nclust-1)*T rather than nclust*T.
* 1.0.14  (04sep2018)
*         Bug with initial (ridge) beta if XX+lambda2/2*diag(Ups2) singular, caused by using lusolve.
*         Fixed with subroutine luqrsolve(.). Uses lusolve by default; if singular, then uses qrsolve.
*         Dropped unpenalized variables now has more informative error message.
*         Modified code to accommodate cluster variables that are strings.
*         Added weighting subroutine (data assumed to be preweighted).
*         Consistent use of d.dmflag and d.cons:
*           dmflag=1 => treat data as zero mean
*           cons=1 => constant in model to be recovered by estimation code; also used in variable counts
*           dmflag=1 => cons=0 but cons=0 /=> dmflag=1 because model may not have a constant
*		  Fixed wrong formula for R-squared in getInfoCriteria(). RSQ =1:-RSS:/TSS (was ESS:/TSS).
*         Fixed initial resids routine for rlasso for model with no constant.
*         Fixed bug in sup score test + prestd. Fixed bug that didn't allow sup score test + pnotpen(.).
*         Added separate storage of prestandardization SDs on data struct as prestdx and prestdy.
*		  Added "ebicgamma" option. Default choice is now ebicgamma=1-1/(2*kappa) and p=n^kappa.
*			See Chen & Chen (2008, Sec. 5, p. 768).
* 1.1.01  (08nov2018)
*         Rewrite of supscore code. Now reports sqrt(n)*[] rather than n*[] version as per CCK.
*           Bug fix in cluster version (was using nclust*[] rather than n*[]).
*           Now handles weighted/demeaned data, center option correctly.
*         Rewrite of lambdacalc code. Shares subroutines with supscore code.
*         Rewrite of xdep and lasso/sqrt-lasso weights code. Shares subroutines with supscore code.
*           sqrt-lasso xdep simulation now normalizes by empirical SD of simulated N(0,1).
*           Bug fix in xdep gamma; distribution was based on max(Xe) instead of max(abs(Xe));
*           effect was appx equivalent to using a gamma 1/2 times the value specified; now fixed.
*         Rewrite of cluster weights code to use panelsum(.) and avoid large temp matrices.
*         Added saved value of minimized objective function (standardized).
*         Initial weights for sqrt-lasso + inid data as per BCW 2014 p. 769 = colmax(abs(X))
*           now available as undocumented maxabsx option.  For cluster, use colmax(panelmean(abs(X)).
*           maxabsx => default maxupsiter=3 in line with BCW 2014 recommendation to iterate.
*         Bug fix in reporting of sqrt-lasso with vverbose (wasn't reporting correct value of obj fn).
*         Bug fix in sqrt-lasso with nocons (was demeaning despite no constant present).
*         Bug fix in rlasso verbose reporting - wasn't reporting selected variables.
*         Removed factor of 0.5 from initial penalty (loadings) for lasso + post-lasso OLS resids.
*         Added undocumented c0(.) option to replicate previous factor of 0.5 for initial penalty loadings.
*         Added undocumented sigma(.) option for std lasso to override estimation of initial sigma.
*         Added undocumented ssiid option for supscore test to override use of mult bootstrap in iid case.
* 1.1.02  (2jan2018)
*         Saved value of minimized objective function is now in y units (not standardized).
* 		  Adaptive lasso: added check for zeros in e(b). Use univariate OLS instead.
* 1.1.03 (13jan2019)
*         Updated options and terminology from Ups to Psi (penalty loadings).
*         gamma option now controls overall level, = 0.1/log(N) or 0.1/log(N_clust) by default
*         gammad option now undocumented
* 1.1.04 (14oct2019)
*         Bug fix - lambda0(.) option was being ignored for standard lasso.
*         Algorithm update for rlasso - now maxupsiter=1 will exit with the initial penalty loadings,
*         initial lambda, and lasso estimates based on this.
*         NB: Also specifying lambda0(.) will control the initial lambda.
* 1.2.01 (29july2020)
*         Substantial rewrite and simplification of code, misc bug fixes, new features.
*         HAC/AC rlasso added. Imports code from livreg2 with modifications.
*         2-way cluster-robust added. rlasso now also returns r2.
*         DoSqrtLasso(.) now incorporated into DoLasso(.); similarly for DoSqrtLassoPath(.).
*         Shooting algo convergence now in standardized metric
*         DoLassoPath(.) lambda iteration report; small coefs not set to zero if rige; compares to std coeffs;
*            saves untruncated lambda list; dof and model dimension fixes (counts all partialled out or none);
*            exits path early if dev criteria met; override and obtain old behavior with nodevcrit option.
*         Bug fix in lassopath with stdcof+supplied lamba.
*         lassopath now exits if lminratio is illegal.
*         lassopath saves shat0 (shat+partialled out) and lambdamat0 (untruncated lambdamat).
*         Bug fixes in info criteria code (check for constant, handling of sqrt-lasso, demeaning of y).
*         Effective degrees of freedom bug fix (accounts for constant and partialled-out).
*         Adalasso fixes: initialization of zerosinada; scaling of adaptive weights; varnames checks;
*            initialization of Psi (default 1e-9 rather than exact 0); adaptive now matches adaptive+prestd;
*            abs(.) precedes raising to theta power; verbose behavior
*         Unpartial bug fix (could crash if no selected regressors).
*         Solver option allowed for partialling out; recoded to be string rather than numeric.
*         Used for both partialling out by calling program and for rlasso partialling.
*         Default solver for partialling is now qrxx = QR+quadcross (was SVD+QR).
*         Warning issued if collinearities encountered when partialling out.
*         Bug fix in case all ICs were missing.
*         nodevcrit option (uses full lambda list irrespective of size of deviance in path)
* 1.2.02 (4aug2020)
*         Bug fix in adaptive lasso relating to omitted vars (would get a missing penalty loading).
*         Reverted to adaptive lasso behavior where any omitted vars trigger univariate OLS.
*         Bug fix in lmax+prestd option (lmax should be rescaled).
*         Reporting bug fix with vverbose option and elastic net - wrong obj fn. (Reporting only.)
* 1.2.03 (27sept2020)
*         Added lglmnet option for _lassopath.
*         Recoded algorithms to drop use of lambda2 and psi2.
*         Removed * for options so illegals caught.
*         Bug fix in value returned for ridge objective function (was just rmse, excluded penalty).
*         Bug fix in starting ridge coefs for rlasso.
*         Added check for negative penalty loadings.
*         Bug fix for user-supplied penalty loadings and adaptive lasso with alpha<1.
*         Calculation of intercepts moved to return results subroutines.
*         Adaptive lasso recoded in Mata.
*         Adaptive elastic net now correctly uses adaptive weights for L1 norm only.
*         Standardized and unstandardized coefficients now returned for lassopath; stdcoef option unneeded.
*         Standardized lambdas, norms, ICs, etc. also returned.
*         lassopath code takes stdall option - indicates that provided lambdas are in the standardized metric
*         Fixed bug in reporting of value of maximized obj fn for elastic net/ridge (missing 1/2 on L2 norm).
*         EBIC now excludes omitted/base variables when calculating model size p.
* 1.2.04  Bug fix to SD calculation for special case of lglmnet with unit loadings (=not prestandardized)


program lassoutils, rclass sortpreserve
	version 13
	syntax [anything] ,							/// anything is actually a varlist but this is more flexible - see below
		[										/// 
		rlasso									/// branches to _rlasso
		path									/// branches to _lassopath
		unpartial								/// branches to _unpartial
		fe(string)								/// branches to _fe
		std										/// branches to _std
		wvar(string)							/// branches to _wt
		partialflag(int 0)						/// branches to _partial if =1
		partial(string)							/// used by _partial and _unpartial
		tvarlist(string)						/// used by _partial, _fe, _std; optional list of temp vars to fill
		tpvarlist(string)						/// used by _partial; optional list of temp partial vars to file
												///
		ALPha(real 1) 							/// elastic net parameter
		SQRT									/// square-root-lasso
		ols 									///
		adaptive								///
												///
												/// verbose options will be transformed to 0/1/2 before passing to subs
		VERbose									/// shows additional detail
		VVERbose								/// even more detail
		*										///
		]

	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename.
	local varlist `anything'
	*

	** verbose
	if ("`vverbose'"!="") {
		local verbose=2							//  much detail
	} 
	else if ("`verbose'"!="") {
		local verbose=1							//  some additional detail
	}
	else {
		local verbose=0							//  no additional output
	}
	*

	************* SUBROUTINES: shoot, path, CV **************************************

	if "`rlasso'" ~= "" {
		_rlasso `varlist',							 	/// branch to _rlasso
			verbose(`verbose')							///
			`sqrt'										///
			`options'									//  no penloads etc.
	}
	else if "`path'" != "" {
		_lassopath `varlist',							///  branch to _ lassopath
			verbose(`verbose')							///
			`sqrt'										///
			alpha(`alpha')								///
 			`adaptive'									///
			`ols'										///
			`options'
	}
	else if "`unpartial'" ~= "" {
		_unpartial `varlist',							/// branch to _unpartial
			partial(`partial')							///
			wvar(`wvar')								///
			`options'
	}
	else if `partialflag' {
		_partial `varlist',							 	/// branch to _partial
			partial(`partial')							///
			tvarlist(`tvarlist')						///
			wvar(`wvar')								///
			`options'
	}
	else if "`fe'" ~= "" {
		_fe `varlist',								 	/// branch to _fe
			tvarlist(`tvarlist')						///
			fe(`fe')									///
			wvar(`wvar')								///
			`options'
	}
	else if "`std'" ~= "" {
		_std `varlist',								 	/// branch to _std
			tvarlist(`tvarlist')						///
			`options'
	}
	else if "`wvar'" ~= "" {
		_wt `varlist',								 	/// branch to _wt
			tvarlist(`tvarlist')						///
			wvar(`wvar')								///
			`options'
	}
	else {
		di as err "internal lassoutils error"
		exit 499
	}

	*** return
	// so subroutine r(.) detail passed on
	// can add additional items to r(.) here
	return add

	// if estimation, return method etc.
	if "`rlasso'`path'`cv'" ~= "" {
		return scalar alpha			= `alpha'
		return scalar sqrt			= ("`sqrt'"!="")
		
		if ("`sqrt'" ~= "") {
			return local method		"sqrt-lasso"
		}
		else if `alpha'==1 {
			return local method		"lasso"
		}
		else if `alpha'==0 {
			return local method		"ridge"
		}
		else if (`alpha'<1) & (`alpha'>0)  {
			return local method		"elastic net"
		}
		else {
			di as err "internal lassoutils error - alpha outside allowed range"
			exit 499
		}
	}

end

// Subroutine for BCH _rlasso
program define _rlasso, rclass sortpreserve
	version 13
	syntax anything ,							/// anything is actually a varlist but this is more flexible - see below
		touse(varlist numeric max=1)			/// required `touse' variable (so marksample unneeded)
		[										///
		verbose(int 0)							///
												///
		ROBust									/// penalty level/loadings allow for heterosk. [homoskedasticity is the default]
		CLuster(varlist max=2)					/// penalty level/loadings allow for within-panel dependence & heterosk.
		bw(int 0)								/// HAC bandwidth
		kernel(string)							///
		tindex(varlist max=1)					///
		bsize(int 5)							/// default block size for the HAC/AC multiplier bootstrap
		maq										/// specific to HAC/AC with truncated kernel
		XDEPendent								/// penalty level is estimated depending on X (default: No)
		numsim(integer 5000)					/// number of simulations with xdep (default=5,000)
												/// 
		TOLOpt(real 1e-10)						/// was originally 1e-5 but 1e-10 gives good equiv between xstd and stdloadings
		TOLPsi(real 1e-4)						///
		TOLZero(real 1e-4)						///
		MAXIter(integer 10000)					/// max number of iterations in estimating lasso
		MAXPSIIter(int 2)						/// max number of lasso-based iterations in est penalty loadings; default=2
		CORRNumber(int 5) 						/// number of regressors used in InitialResiduals(); default=5
		maxabsx									/// initial loadings for sqrt-lasso as per BCW 2014 p. 769 = colmax(abs(X))
												///
		LASSOPSI								/// use lasso residuals to estimate penalty loadings (post-lasso is default)
												///
		c(real 1.1)								/// "c" in lambda function; default = 1.1
		c0(real 0)								/// when iterating, initial value of c in lambda function; default set below
		gamma(real 0)							/// "gamma" in lambda function (default 0.1/log(N) or 0.1/log(N_clust))
		gammad(real 0)							/// undocumented "gamma" denominator option (default log(N) or log(N_clust)
		lambda0(real 0)							/// optional user-supplied lambda0
		LALTernative 							/// use alternative, less sharp lambda0 formula
		PMINus(int 0)							/// dim(X) adjustment in lambda0 formula
		nclust0									/// no longer in use as of 1.0.08; replaced by nclust1
		nclust1									/// use (nclust-1)*T instead of nclust*T in cluster-lasso
		center									/// center x_i*e_i or x_ij*e_ij
		sigma(real 0)							/// user-specified sigma for rlasso
												///
		supscore								///
		ssgamma(real 0.05)						/// default gamma for sup-score test
		ssnumsim(integer 500)					///
		testonly								///
		ssiid									/// undocumented option - override use of mult bootstrap in iid case
												///
		seed(real -1)							/// set random # seed; relevant for xdep and supscore
		dots									/// display dots for xdep and supscore repetitions
												///
		xnames_o(string)						/// original names in varXmodel if tempvars used (string so varlist isn't processed)
		xnames_t(string)						/// temp names in varXmodel
		notpen_o(string)						/// used in checking
		notpen_t(string)						///
		consflag(int 0)							/// =0 if no cons or already partialled out
		dofminus(int 0)							/// lost degrees of freedom from FEs
		sdofminus(int 0)						/// lost degrees of freedom from partialling out ("s"=small; includes the constant)
		dmflag(int 0)							/// data have been demeaned or should be treated as demeaned
												///
		stdy(string)							///
		stdx(string)							///
												///
		SQRT									/// square-root lasso
												///
												/// LEGACY OPTIONS
		TOLUps(real 1e-4)						///
		MAXUPSIter(int 2)						/// max number of lasso-based iterations in est penalty loadings; default=2
		LASSOUPS								/// use lasso residuals to estimate penalty loadings (post-lasso is default)
		psolver(string)							///
		]

	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename.
	local varlist `anything'
	*

	*** count number of obs
	_nobs `touse'
	local nobs = r(N)
	*

	*** cluster
	// create numeric clustvar in case cluster variable is string
	// sort needed if cluster option provided; will need to re-sort afterwards
	local clustercount	: word count `cluster'
	if `clustercount' {
		local sortedby	: sortedby
	}
	if `clustercount'==1 {
		tempvar clustvar
		qui egen `clustvar'		= group(`cluster') if `touse'
		sort `clustvar'
		local allclust			`clustvar'
	}
	else if `clustercount'==2 {
		tokenize `cluster'
		tempvar clustvar clustvar2 clustvar3
		qui egen `clustvar'		= group(`1') if `touse'
		qui egen `clustvar2'	= group(`2') if `touse'
		qui egen `clustvar3'	= group(`1' `2') if `touse'
		sort `clustvar3' `clustvar'
		local allclust			`clustvar' `clustvar2' `clustvar3'
	}
	*
	
	*** HAC/AC
	if `bw'~=0 {
		capture tsset
		if "`r(timevar)'" == "" {
			di as err "must tsset data and specify timevar"
			exit 5
		}
		local tdelta	= `r(tdelta)'
		// check if valid bandwidth
		sum `r(timevar)' if `touse', meanonly
		local T = r(max)-r(min) + 1
		if `bw' > (`T'-1)/`tdelta' {
			di as err "invalid bandwidth in option bw() - cannot exceed timespan of data"
			exit 198
		}
		// check if valid kernel, standardize kernel name, also set spectral flag
		mata: s_vkernel("`kernel'", "`bw'")
		local kernel	`r(kernel)'
		local spectral	= `r(spectral)'
		local maqflag	= ("`maq'"~="" & "`kernel'"=="Truncated")
	}
	else {
		// bw=0 means ignore kernel option and set spectral to default of 0
		local tdelta	= 0
		local kernel
		local spectral	= 0
		local maqflag	= 0
	}
	*
	
	// incompatible options
	if `clustercount' & `bw' {
		di as err "incompatible options - cluster(.) and bw(.)"
		exit 198
	}

	*** define various parameters/locals/varlists
	local varY			`varlist'							// now is only element in varlist
	local varXmodel		`xnames_t'							// including notpen if any; may include tempnames
	local pen_t			: list xnames_t - notpen_t			// list of penalized regressors
	local p0			: word count `varXmodel' 			// #regressors incl notpen but excl constant
	local p				= `p0' - `pminus'					// p defined by BCH
	local p0			= `p0' + `consflag'					// #regressors but now incl constant
	*
	
	*** define various flags
	local clustflag		=("`cluster'"!="")					// =1 if cluster
	local hetero		=("`robust'"!="") & ~`clustflag'	// =1 if het-robust but NOT cluster
	local sqrtflag		=("`sqrt'"!="")
	local xdep			=("`xdependent'"!="")
	local lassopsiflag	=("`lassopsi'`lassoups'"!="")		// lassoups is the equivalent legacy option
	local lalternative	=("`lalternative'"!="")				// =1 if alternative, less sharp lambda0 should be used
	local nclust1flag	=("`nclust1'"!="")
	local center		=("`center'"!="")
	local supscoreflag	=("`supscore'`testonly'"!="")
	local testonlyflag	=("`testonly'"!="")
	local dotsflag		=("`dots'"!="")
	local prestdflag	=("`stdy'"!="")
	local ssiidflag		=("`ssiid'"!="")
	local maxabsxflag	=("`maxabsx'"!="")
	*
	
	*** defaults
	if `maxabsxflag' {
		// implement BCW 2014 p. 769 recommendation for sqrt lasso initial weights = colmax(abs(X))
		local corrnumber -1
		local maxpsiiter 3		// since BCW say to iterate
	}
	if `c0'==0 {
		// allow initial value of c in iterated penalties to vary; matches treatment in CBH and hdm code
		local c0 `c'
	}
	*

	*** legacy options
	if `tolups' ~= 1e-4 {
		local tolpsi `tolups'
	}
	if `maxupsiter' ~= 2 {
		local maxpsiiter `maxupsiter'
	}
	*

	tempname lambda slambda rmse rmseOLS objfn r2	//  note lambda0 already defined as a macro
	tempname b bOLS sb sbOLS Psi sPsi ePsi stdvec
	tempname bAll bAllOLS
	tempname supscoremat CCK_ss CCK_p CCK_cv CCK_gamma
	// initialize so that returns don't crash in case values not set
	local k				=0						//  used as flag to indicate betas are non-missing
	local N				=.
	local N_clust		=.
	local s				=.
	local s0			=.
	local niter			=.
	local npsiiter		=.
	scalar `lambda'		=.
	scalar `slambda'	=.
	scalar `rmse'		=.
	scalar `rmseOLS'	=.
	scalar `r2'			=.
	scalar `objfn'		=.
	mat `b'				=.
	mat `bAll'			=.
	mat `bOLS'			=.
	mat `bAllOLS'		=.
	mat `sb'			=.
	mat `sbOLS'			=.
	mat `Psi'			=.
	mat `sPsi'			=.
	mat `ePsi'			=.
	mat `stdvec'		=.
	mat `supscoremat'	=.
	scalar `CCK_ss'		=.
	scalar `CCK_p'		=.
	scalar `CCK_cv'		=.
	scalar `CCK_gamma'	=.

	if `p' & ~`testonlyflag' {					//  there are penalized model variables; estimate
	
		*** BCH rlasso
		mata:	EstimateRLasso(					///
							"`varY'",			///
							"`varXmodel'",		///
							"`xnames_o'",		///
							"`pen_t'",			///
							"`notpen_t'",		/// #5
							"`touse'",			///
							"`allclust'",		///
							"`stdy'",			///
							"`stdx'",			///
							`sqrtflag',			/// #10
							`hetero',			///
							`bw',				///
							"`kernel'",			///
							`maqflag',			///
							`spectral',			/// #15
							"`tindex'",			///
							`tdelta',			///
							`bsize',			///
							`xdep',				///
							`numsim',			/// #20
							`lassopsiflag',		///
							`tolopt',			///
							`maxiter',			///
							`tolzero',			///
							`maxpsiiter',		/// #25
							`tolpsi',			///
							`verbose',			///
							`c',				///
							`c0',				///
							`gamma',			/// #30
							`gammad',			///
							`lambda0',			///
							`lalternative',		///
							`corrnumber',		///
							`pminus',			/// #35
							`nclust1flag',		///
							`center',			///
							`sigma',			///
							`supscoreflag',		///
							`ssnumsim',			/// #40
							`ssgamma',			///
							`seed',				///
							`dotsflag',			///
							`consflag',			///
							`dmflag',			/// #45
							`dofminus',			///
							`sdofminus',		///
							`prestdflag',		///
							"`psolver'")

		*** Via ReturnResults(.)
		// coefs are returned as column vectors
		// convert to row vectors (Stata standard)
		mat `b'				=r(b)'					//  in original units
		mat `bOLS'			=r(bOLS)'
		mat `sb'			=r(sb)'					//  in standardized units
		mat `sbOLS'			=r(sbOLS)'
		mat `bAll'			=r(bAll)'
		mat `bAllOLS'		=r(bAllOLS)'
		mat `Psi'			=r(Psi)
		mat `sPsi'			=r(sPsi)
		mat `ePsi'			=r(ePsi)
		mat `stdvec'		=r(stdvecp)				//  stdvec for penalized vars only
		local selected0		`r(sel)'				//  selected variables INCL NOTPEN, EXCL CONSTANT
		local s0			=r(s)					//	number of selected vars INCL NOTPEN, EXCL CONSTANT; may be =0
		local k				=r(k)					//  number of all vars in estimated parameter vector INC CONSTANT; may be =0
		local niter			=r(niter)
		local npsiiter		=r(npsiiter)
		local N				=r(N)
		local N_clust		=r(N_clust)
		local N_clust1		=r(N_clust1)
		local N_clust2		=r(N_clust2)
		scalar `lambda'		=r(lambda)				//  relates to depvar in original units
		scalar `slambda'	=r(slambda)				//  relates to standardized depvar
		scalar `rmse'		=r(rmse)				//  lasso rmse
		scalar `rmseOLS'	=r(rmsePL)				//  post-lasso rmse
		scalar `r2'			=r(r2)
		scalar `objfn'		=r(objfn)				//  minimized objective function
		if `supscoreflag' {
			mat `supscoremat'	=r(supscore)			//  sup-score row vector of results
			scalar `CCK_ss'		=`supscoremat'[1,colnumb(`supscoremat',"CCK_ss")]
			scalar `CCK_p'		=`supscoremat'[1,colnumb(`supscoremat',"CCK_p")]
			scalar `CCK_cv'		=`supscoremat'[1,colnumb(`supscoremat',"CCK_cv")]
			scalar `CCK_gamma'	=`supscoremat'[1,colnumb(`supscoremat',"CCK_gamma")]
		}
		// these may overwrite existing locals
		tempname lambda0 c gamma gammad
		scalar `lambda0'	=r(lambda0)				//  BCH definition of lambda; overwrites existing/default macro
		scalar `c'			=r(c)
		scalar `gamma'		=r(gamma)
		scalar `gammad'		=r(gammad)
		// for HAC or 2-way cluster lasso (neg loadings possible)
		local psinegs		=r(psinegs)
		local psinegvars	`r(psinegvars)'
		*
	}
	else if ~`testonlyflag' {						//  there are no penalized model vars; just do OLS

		if `p0' {									//  there are unpenalized vars and/or constant
			di as text "warning: no penalized variables; estimates are OLS"
		}
		if ~`consflag' {
			local noconstant noconstant
		}
		else {
			local consname _cons
		}
		qui regress `varY' `varXmodel' if `touse', `noconstant'

		mat `b'					=e(b)
		local k					=colsof(`b')
		foreach m in `b' `bOLS' `bAll' `bAllOLS' {
			mat `m'				=e(b)
			mat colnames `m'	=`xnames_o' `consname'
		}
		scalar `rmse'			=e(rmse)
		scalar `rmseOLS'		=e(rmse)
		scalar `r2'				=e(r2)
		local selected0			`xnames_o'
		local N					=e(N)
		local N_clust			=.
		local N_clust1			=.
		local N_clust2			=.
	}
	else if `testonlyflag' {					//  supscore test only

		mata:	EstimateSupScore(				///
							"`varY'",			///
							"`varXmodel'",		///
							"`xnames_o'",		///
							"`pen_t'",			///
							"`notpen_t'",		/// #5
							"`touse'",			///
							"`allclust'",		///
							"`stdy'",			///
							"`stdx'",			///
							`sqrtflag',			/// #10
							`hetero',			///
							`bw',				///
							"`kernel'",			///
							`maqflag',			///
							`spectral',			/// #15
							"`tindex'",			///
							`tdelta',			///
							`bsize',			///
							`verbose',			///
							`ssnumsim',			/// #20
							`ssiidflag',		///
							`c',				///
							`nclust1flag',		///
							`center',			///
							`ssgamma',			/// #25
							`pminus',			///
							`seed',				///
							`dotsflag',			///
							`consflag',			///
							`dmflag',			/// #30
							`dofminus',			///
							`sdofminus',		///
							`prestdflag',		///
							"`psolver'")

		mat  `supscoremat'	=r(supscore)
		scalar `CCK_ss'		=`supscoremat'[1,colnumb(`supscoremat',"CCK_ss")]
		scalar `CCK_p'		=`supscoremat'[1,colnumb(`supscoremat',"CCK_p")]
		scalar `CCK_cv'		=`supscoremat'[1,colnumb(`supscoremat',"CCK_cv")]
		scalar `CCK_gamma'	=`supscoremat'[1,colnumb(`supscoremat',"CCK_gamma")]
		local N				=r(N)
		local N_clust		=r(N_clust)
		local N_clust1		=r(N_clust1)
		local N_clust2		=r(N_clust2)

		// for HAC or 2-way cluster lasso (neg loadings possible)
		local psinegs		=r(psinegs)
		local psinegvars	`r(psinegvars)'
	}
	else {										//  shouldn't reach here
		di as err "internal _rlasso error"
		exit 499
	}

	**************** Estimation complete ********************
	
	// cluster required sorting, so need to restore original sort if was originally sorted
	if `clustercount' & "`sortedby'"~="" {
		sort `sortedby'
	}
	
	*** check notpen is in selected0
	if ~`testonlyflag' {
		fvstrip `selected0'							//  fvstrip removes b/n/o prefixes (shouldn't need dropomit)
		local list1			`r(varlist)'
		fvstrip `notpen_o', dropomit				//  use dropomit option here
		local list2			`r(varlist)'
		local missing		: list list2 - list1
		local missing_ct	: word count `missing'
		if `missing_ct' {
			di as err "internal _rlasso error - unpenalized `missing' missing from selected vars"
			di as err "may be caused by collinearities in unpenalized variables or by tolerances"
			di as err "set tolzero(.) or other tolerances smaller or use partial(.) option"
			exit 499
		}
	}
	*

	*** conventions
	// k			= number of selected/notpen INCLUDING constant; k>0 means estimated coef vector is non-empty
	// s0			= number of selected/notpen EXCLUDING constant
	// s			= number of selected (penalized)
	// p0			= number of all variables in model INCLUDING constant
	// selected0	= varlist of selected and unpenalized EXCLUDING constant; has s0 elements
	// selected		= varlist of selected; has s elements
	// cons			= 1 if constant, 0 if no constant
	// notpen_ct	= number of unpenalized variables in the model EXCLUDING CONSTANT
	// coef vectors b, bOLS, sb, sbOLS have k elements
	// coef vectors bAll, bAllOLS have p0 elements
	// Psi, sPsi have p0-cons elements
	// stdvec has p0 - notpen_ct - cons elements (=#penalized)
	// Note that full set of regressors can including omitteds, factor base categories, etc.
	*

	*** fix colnames of beta vectors to include omitted "o." notation
	// includes post-lasso vector in case of collinearity => a selected variable was omitted
	// trick is to use _ms_findomitted utility but give it
	// diag(bAll) as vcv matrix where it looks for zeros
	// also build in fv info
	tempname tempvmat
	if `k' & ~`testonlyflag' {							//  if any vars in est coef vector (k can be zero)
		mat `tempvmat'	= diag(`bOLS')
		_ms_findomitted	`bOLS' `tempvmat'
		_ms_build_info	`bOLS' if `touse'
	}
	if `p0' & ~`testonlyflag' {							//  if any variables in model (p0 can be zero)
		mat `tempvmat'	= diag(`bAll')
		_ms_findomitted	`bAll' `tempvmat'
		_ms_build_info	`bAll' if `touse'
		mat `tempvmat'	= diag(`bAllOLS')
		_ms_findomitted	`bAllOLS' `tempvmat'
		_ms_build_info	`bAllOLS' if `touse'
	}
	*

	*** manipulate variable counts and lists
	// selected0 and s0 are selected+notpen (excluding constant)
	// selected and s are selected only
	local notpen_ct		: word count `notpen_o'			//  number of notpen EXCL CONSTANT
	local selected		: list selected0 - notpen_o		//  selected now has only selected/penalized variables
	local s				: word count `selected'			//  number of selected/penalized vars EXCL NOTPEN/CONSTANT
	*
	*** error checks
	// check col vector of b (selected) vs. k
	local col_ct		=colsof(`b')
	if `col_ct'==1 & el(`b',1,1)==. & ~`testonlyflag' {
		// coef vector is empty so k should be zero
		if `k'~=0 {
			di as err "internal _rlasso error - r(k)=`k' does not match number of selected vars/coefs=`col_ct'"
			exit 499
		}
	}
	else if `k'~=`col_ct' & ~`testonlyflag' {
		// coef vector not empty so k should match col count
		di as err "internal _rlasso error - r(k)=`k' does not match number of selected vars/coefs=`col_ct'"
		exit 499
	}
	// check col vector of bAll vs. p0
	local col_ct		=colsof(`bAll')
	if `p0'~=`col_ct' & ~`testonlyflag' {
		// full coef vector not empty so p0 should match col count
		di as err "internal _rlasso error - p0=`p0' does not match number of model vars/coefs=`col_ct'"
		exit 499
	}
	*

	*** return results
	// coef vectors are row vectors (Stata standard)
	return matrix beta		=`b'
	return matrix betaOLS	=`bOLS'
	return matrix sbeta		=`sb'
	return matrix sbetaOLS	=`sbOLS'
	return matrix betaAll	=`bAll'
	return matrix betaAllOLS=`bAllOLS'
	return matrix Psi		=`Psi'
	return matrix sPsi		=`sPsi'
	return matrix ePsi		=`ePsi'
	return matrix stdvec	=`stdvec'					//  penalized vars only
	return scalar lambda	=`lambda'					//  relates to depvar in original units
	return scalar slambda	=`slambda'					//  relates to standardized depvar
	return scalar lambda0	=`lambda0'					//  BCH definition of lambda
	return scalar rmse		=`rmse'						//  lasso rmse
	return scalar r2		=`r2'
	return scalar objfn		=`objfn'
	return scalar c			=`c'
	return scalar gamma		=`gamma'
	return scalar gammad	=`gammad'
	return scalar rmseOLS	=`rmseOLS'					//  post-lasso rmse
	return local  selected0	`selected0'					//  all selected/notpen vars INCLUDING NOTPEN (but excl constant)
	return local  selected	`selected' 					//  all selected (penalized) vars EXCLUDING NOTPEN & CONS
	return scalar k			=`k'						//  number of all vars in sel/notpen parameter vector INCLUDING CONSTANT
	return scalar s0		=`s0'						//  number of vars selected INCLUDING NOTPEN (but excl constant)
	return scalar s			=`s'						//  number of vars selected EXCLUDING NOTPEN & CONS
	return scalar p0		=`p0'						//  number of all vars in original model including constant
	return scalar p			=`p'						//  p defined by BCH; excludes notpen & cons
	return scalar N_clust	=`N_clust'
	return scalar N_clust1	=`N_clust1'
	return scalar N_clust2	=`N_clust2'
	return scalar N			=`N'
	return scalar center	=`center'
	
	// for HAC lasso
	return scalar bw		=`bw'
	return local kernel		`kernel'
	return local psinegs	=`psinegs'
	return local psinegvars	`psinegvars'

	return local  clustvar	`cluster'
	return local  robust	`robust'
	return scalar niter		=`niter'
	return scalar maxiter	=`maxiter'
	return scalar npsiiter	=`npsiiter'
	return scalar maxpsiiter=`maxpsiiter'
	
	return scalar ssnumsim		=`ssnumsim'
	return scalar supscore		=`CCK_ss'
	return scalar supscore_p	=`CCK_p'
	return scalar supscore_cv	=`CCK_cv'
	return scalar supscore_gamma=`CCK_gamma'
	*

end		// end _rlasso subroutine

// Subroutine for lassopath
program define _lassopath, rclass sortpreserve
	version 13
	syntax [anything] ,							/// anything is actually a varlist but this is more flexible - see below
		toest(varlist numeric max=1)			/// required `touse' variable (so marksample unneeded)
		[										/// 
		verbose(int 0)							///
		consflag(int 1)							/// default is constant specified
		stdallflag(int 0)						/// =1 if lambdas etc. are provided in the standardized metric
		dofminus(int 0)							/// lost degrees of freedom from FEs
		sdofminus(int 0)						/// lost degrees of freedom from partialling out ("s"=small; includes the constant)
		pminus(int 0)							/// omitted/base variables to be excluded from model count
		dmflag(int 0)							/// data have been demeaned or should be treated as demeaned
		notpen_t(string) 						///
		notpen_o(string) 						///
												///
												/// lambda settings
		ALPha(real 1) 							/// elastic net parameter
		Lambda(string)							/// overwrite default lambda
		lambda2(string)							/// overwrite default L2 lambda
		LCount(integer 100)						///
		LMAX(real 0) 							///
		LMINRatio(real 1e-4)					/// ratio of maximum to minimum lambda
												///
		TOLOpt(real 1e-10)						/// was originally 1e-5 but 1e-10 gives good equiv between xstd and stdloadings
		TOLZero(real 1e-4)						///
		MAXIter(integer 10000)					/// max number of iterations in estimating lasso
		fdev(real 1e-5)							/// minimum fractional change in deviance for stopping path a la glmnet
		devmax(real 0.999)						/// maximum fraction of explained deviance for stopping path a la glmnet
												///
		PLoadings(string)						/// set penalty loadings as Stata row vector
		ploadings2(string)						/// set optional L2 penalty loadings as Stata row vector
												///
		LGLMnet									/// use glmnet scaling for lambda/alpha
												///
		xnames_o(string)						/// original names in varXmodel if tempvars used (string so varlist isn't processed)
		xnames_t(string)						/// temp names in varXmodel
												///
		sqrt 									/// square-root lasso
		ols										///
												///
  		STDLflag(int 0) 						/// use standardisation loadings?
		stdy(string)							///
		stdx(string)							///
												///
		ADAPTive  								/// adaptive lasso
		ADATheta(real 1) 						/// gamma paramater for adapLASSO
		ADALoadings(string)						///
		ADAFirst(string)						/// specify first-step ada estimator
												///
		holdout(varlist) 						///
												///
		NOIC									///
		EBICGamma(real -99)						/// -99 leads to the default option
												///
		NODEVCRIT								/// do not use deviance as criteria to exit path
		]
		
	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename to `varlist' (usual Stata convention).
	local varlist `anything'

	*** quietly unless vverbose specified
	if `verbose'<2 {
		local quietly quietly
	}
	*
	
	*** count number of obs
	_nobs `toest'
	local nobs = r(N)
	*

	*** define various parameters/locals/varlists
	local varY			`varlist'						// now is only element in varlist
	local varXmodel		`xnames_t'						// including notpen if any; may include tempnames
	local pen_t			: list xnames_t - notpen_t		// list of penalized regressors
	local pmodel		: word count `varXmodel' 		// # regressors (excl. constant) minus partial
	local p0			: word count `varXmodel' 		// # all regressors incl notpen but excl constant
	local p0			= `p0' + `consflag'				// ditto but now incl constant
	*

	*** syntax checks
	if (`lminratio'<0) | (`lminratio'>=1) {
		di as err "lminratio should be greater than 0 and smaller than 1."
		exit 198
	}
	// txt rather than err in order to stop message unless called by lasso2
	if ("`lambda'"!="") & ((`lcount'!=100) | (`lmax'>0) | (`lminratio'!=1e-4)) {
		di as txt "lcount | lmax | lminratio option ignored since lambda() is specified."
	}
	if (`ebicgamma'!=-99) & (`ebicgamma'<0) & (`ebicgamma'>1) {
		di as err "ebicgamma must be in the interval [0,1]. Default option is used."
	}
	*
	
	*** useful flags and counts
	local prestdflag	=("`stdy'"!="")
	local sqrtflag		=("`sqrt'"!="")
	local olsflag		=("`ols'"!="")
	local adaflag 		=("`adaptive'"!="")
	local lglmnetflag	=("`lglmnet'"!="")
	local noicflag		=("`noic'"!="")
	if (~`consflag') {
		local noconstant noconstant
	}
	local nodevcritflag	=("`nodevcrit'"!="")
	*
	
	*** syntax checks
	if `sqrtflag' & `lglmnetflag' {
		di as err "lglmnet option not supported for square-root lasso"
		exit 198
	}
	
	*** special case - lglmnet with unit loadings (=not prestandardized)
	if `lglmnetflag' & "`stdy'"==""  {
		tempvar tY
		tempname stdy stdx
		qui sum `varY' if `toest'
		local sdY			= r(sd) * sqrt((r(N)-1)/r(N))
		local mY			= r(mean)
		if `consflag' {
			gen double `tY'	= `mY' + (`varY' -`mY') * 1/`sdY'  if `toest'
		}
		else {
			gen double `tY'	= `varY' * 1/`sdY'  if `toest'
		}
		local varY			`tY'
		mat `stdy'			= `sdY'
		mat `stdx'			= J(1,`p0'-`consflag',1)
		local prestd		prestd
		local prestdflag	= 1
		local stdlflag		= 0
	}
	*

	*********************************************************************************
	*** penalty loadings and adaptive lasso 									  ***
	*********************************************************************************
	
	*** syntax checks for penalty loadings
	if (`adaflag') & ("`ploadings'"!="") {
		di as error "ploadings() option not allowed with adaptive."
		exit 198
	}
	if (`adaflag') & (`stdlflag') {
		di as error "stdloadings option not allowed with adaptive."
		exit 198
	}
	if (~`adaflag') & ("`adaloadings'"!="") {
		di as error "adaloadings(`adaloadings') ignored. specify adaptive option."
	}
	if (~`adaflag') & (`adatheta'!=1) {
		di as err "adatheta(`adatheta') ignored. specify adaptive option."
	}
	*
	
	*** check dimension of penalty loading vector(s) and rename
	if ("`ploadings'"!="") {
		tempname Psi
		getPsi, ploadings(`ploadings') p(`pmodel')
		mat `Psi'		= r(Psi)
		if ("`ploadings2'")~="" {
			tempname Psi2
			getPsi, ploadings(`ploadings2') p(`pmodel')
			mat `Psi2'	= r(Psi)
		}
	}
	else if ("`adaloadings'"!="") {
		if (colsof(`adaloadings')!=(`pmodel')) | (rowsof(`adaloadings')!=1) {
			di as err "`adaloadings' should be a row vector of length `p'=dim(X) (excl constant)."
			exit 503
		}
	}
	*

	*********************************************************************************
	*** penalty loadings and adaptive lasso END									  ***
	*********************************************************************************
	
	//egen test = rowmiss(`varXmodel') if `toest'
	//sum test

	mata:	EstimateLassoPath(			///
					"`varY'",			///
					"`varXmodel'",		///
					"`xnames_o'",		///
					"`notpen_o'",		///
					"`notpen_t'",		///
					"`toest'",			///
					"`holdout'", 		/// holdout for MSPE
					`consflag',			///
					`stdallflag',		///
					`dmflag',			///
					`dofminus',			/// lost dofs from FEs
					`sdofminus',		/// lost dofs from partialling
					`pminus',			///
					`prestdflag',		///
					"`lambda'",			/// lambda matrix for L1 norm or missing (=> construct default list)
					"`lambda2'",		/// lambda matrix for L2 norm (also optional)
					`lmax',				/// maximum lambda (optional; only used if lambda not specified)
					`lcount',			/// number of lambdas (optional; only used if lambda not specified)
					`lminratio',		/// lmin/lmax ratio (optional; only used if lambda not specified)
					`lglmnetflag',		/// 
					"`Psi'",			/// Optional L1 norm loadings
					"`Psi2'",			/// Optional L2 norm loadings
					"`stdy'",			/// Stata matrix with SD of dep var
					"`stdx'",			/// Stata matrix with SDs of Xs
					`stdlflag',			/// use standardisation loadings  
					`sqrtflag',			/// sqrt lasso  
					`alpha',			/// elastic net parameter
					`adaflag',			/// 
					"`adaloadings'",	/// 
					`adatheta',			/// 
					"`adafirst'",		/// 
					`olsflag',			/// post-OLS estimation
					`tolopt',			///
					`maxiter',			///
					`tolzero',			///
					`fdev',				///
					`devmax',			///
					`verbose',			///
					`noicflag',			///
					`ebicgamma',		///
					`nodevcritflag'		///
					)

	if (`r(lcount)'>1) { //------- #lambda > 1 -----------------------------------------------//
	
		tempname Psi betas sbetas dof lambdamat0 lambdamat slambdamat0 slambdamat
		tempname l1norm sl1norm wl1norm swl1norm stdvec shat shat0
		mat `Psi' 			= r(Psi)
		mat `betas' 		= r(betas)
		mat `sbetas' 		= r(sbetas)
		mat `dof'			= r(dof)
		mat `shat'			= r(shat)
		mat `shat0'			= r(shat0)
		mat `lambdamat'		= r(lambdalist)
		mat `lambdamat0'	= r(lambdalist0)
		mat `slambdamat'	= r(slambdalist)
		mat `slambdamat0'	= r(slambdalist0)
		mat `l1norm'		= r(l1norm)
		mat `sl1norm'		= r(sl1norm)
		mat `wl1norm'		= r(wl1norm)
		mat `swl1norm'		= r(swl1norm)
		mat `stdvec'		= r(stdvec)
		if ("`holdout'"!="") {
			tempname mspe0
			mat `mspe0' = r(mspe)	
		}
		else {
			tempname rss ess tss rsq
			tempname aic aicc bic ebic saic saicc sbic sebic icstd
			mat `rss' = r(rss)
			mat `ess' = r(ess)
			mat `tss' = r(tss)
			mat `rsq' = r(rsq)
			mat `aic' = r(aic)
			mat `bic' = r(bic)
			mat `aicc' = r(aicc)
			mat `ebic' = r(ebic)
			scalar `icstd' = r(icstd)
			mat `saic' = r(saic)
			mat `sbic' = r(sbic)
			mat `saicc' = r(saicc)
			mat `sebic' = r(sebic)
			return scalar aicmin = r(aicmin)
			return scalar aiccmin = r(aiccmin)
			return scalar bicmin = r(bicmin)
			return scalar ebicmin = r(ebicmin)
			return scalar saicmin = r(saicmin)
			return scalar saiccmin = r(saiccmin)
			return scalar sbicmin = r(sbicmin)
			return scalar sebicmin = r(sebicmin)
			return scalar laicid = r(laicid)
			return scalar laiccid = r(laiccid)
			return scalar lbicid = r(lbicid)
			return scalar lebicid = r(lebicid)
			return scalar ebicgamma = r(ebicgamma)
		}
		
		return scalar N				= r(N)
		return scalar lcount 		= r(lcount)
		return scalar olsflag		= `olsflag'
		return scalar lmax 			= r(lmax)
		return scalar lmax0			= r(lmax0)
		return scalar lmin 			= r(lmin)
		return scalar lmin0			= r(lmin0)
		return matrix betas 		= `betas'
		return matrix sbetas 		= `sbetas'
		return matrix dof 			= `dof'
		return matrix shat 			= `shat'
		return matrix shat0			= `shat0'
		return matrix lambdalist 	= `lambdamat'
		return matrix lambdalist0 	= `lambdamat0'
		return matrix slambdalist 	= `slambdamat'
		return matrix slambdalist0 	= `slambdamat0'
		return matrix Psi			= `Psi'
		return matrix l1norm		= `l1norm'
		return matrix sl1norm		= `sl1norm'
		return matrix wl1norm		= `wl1norm'
		return matrix swl1norm		= `swl1norm'
		return matrix stdvec 		= `stdvec'
		if ("`holdout'"!="") {
			return matrix mspe 		= `mspe0'
		}
		else {
			return matrix rss 		= `rss' 
			return matrix ess 		= `ess' 
			return matrix tss		= `tss' 
			return matrix rsq 		= `rsq' 
			return matrix aic 		= `aic' 
			return matrix aicc 		= `aicc' 
			return matrix bic 		= `bic' 
			return matrix ebic 		= `ebic' 
			return matrix saic 		= `saic' 
			return matrix saicc 	= `saicc' 
			return matrix sbic 		= `sbic' 
			return matrix sebic 	= `sebic' 
		}
	}
	else if (`r(lcount)'==1) { 
	
		*** Via ReturnResults(.)
		*** the following code is based on _rlasso
		tempname b bOLS sb sbOLS Psi sPsi stdvec
		tempname bAll bAllOLS sbAll sbAllOLS
		tempname lambda slambda lambda0 rmse rmseOLS objfn sobjfn srmse srmseOLS r2

		// coefs are returned as column vectors
		// convert to row vectors (Stata standard)
		mat `b'				=r(b)'					//  in original units
		mat `bOLS'			=r(bOLS)'
		mat `bAll'			=r(bAll)'
		mat `bAllOLS'		=r(bAllOLS)'
		mat `sb'			=r(sb)'					//  in standardized units
		mat `sbOLS'			=r(sbOLS)'
		mat `sbAll'			=r(sbAll)'
		mat `sbAllOLS'		=r(sbAllOLS)'
		mat `Psi'			=r(Psi)
		mat `stdvec'		=r(stdvec)
		local selected0		`r(sel)'				//  selected variables INCL NOTPEN, EXCL CONSTANT
		local s0			=r(s)					//	number of selected vars INCL NOTPEN, EXCL CONSTANT; may be =0
		local k				=r(k)					//  number of all vars in estimated parameter vector INC CONSTANT; may be =0
		local niter			=r(niter)
		//*//local npsiiter		=r(npsiiter)
		local N				=r(N)
		local N_clust		=r(N_clust)
		scalar `lambda'		=r(lambda)				//  relates to depvar in original units
		scalar `slambda'	=r(slambda)				//  relates to standardized depvar
		scalar `rmse'		=r(rmse)				//  lasso rmse
		scalar `rmseOLS'	=r(rmsePL)				//  post-lasso rmse
		scalar `srmse'		=r(srmse)				//  standardized lasso rmse
		scalar `srmseOLS'	=r(srmsePL)				//  standardized post-lasso rmse
		scalar `r2'			=r(r2)
		scalar `objfn'		=r(objfn)				//  minimized objective function
		scalar `sobjfn'		=r(sobjfn)				//  standardized minimized objective function
		*

		// added to lasso2
		tempname df

		scalar `df'			=r(dof)

		*** check notpen is in selected0
		fvstrip `selected0'							// fvstrip removes b/n/o prefixes.
		local list1			`r(varlist)'
		fvstrip `notpen_o', dropomit				// use dropomit option here
		local list2		`r(varlist)'
		local missing		: list list2 - list1
		local missing_ct	: word count `missing'
		if `missing_ct' {
			di as err "internal _lassopath error - unpenalized `missing' missing from selected vars"
			di as err "set tolzero(.) or other tolerances smaller or use partial(.) option"
			exit 499
		}
		*
		
		*** conventions
		// k			= number of selected/notpen INCLUDING constant; k>0 means estimated coef vector is non-empty
		// s0			= number of selected/notpen EXCLUDING constant
		// s			= number of selected (penalized)
		// p0			= number of all variables in model INCLUDING constant
		// selected0	= varlist of selected and unpenalized EXCLUDING constant; has s0 elements
		// selected		= varlist of selected; has s elements
		// cons			= 1 if constant, 0 if no constant
		// notpen_ct	= number of unpenalized variables in the model EXCLUDING CONSTANT
		// coef vectors b, bOLS, sb, sbOLS have k elements
		// coef vectors bAll, bAllOLS have p0 elements
		// Psi, sPsi have p0-cons elements
		// stdvec has p0 - notpen_ct - cons elements (=#penalized)
		// Note that full set of regressors can including omitteds, factor base categories, etc.
		*

		*** fix colnames of beta vectors to include omitted "o." notation
		// includes post-lasso vector in case of collinearity => a selected variable was omitted
		// trick is to use _ms_findomitted utility but give it
		// diag(bAll) as vcv matrix where it looks for zeros
		// also build in fv info
		tempname tempvmat
		foreach bmat in `bOLS' `bAll' `bAllOLS' `sbAll' `sbAllOLS' {
			mat `tempvmat'	= diag(`bmat')
			_ms_findomitted	`bmat' `tempvmat'
			_ms_build_info	`bmat' if `toest'
		}
		*

		*** manipulate variable counts and lists
		// selected0 and s0 are selected+notpen (excluding constant)
		// selected and s are selected only
		local notpen_ct		: word count `notpen_o'			//  number of notpen EXCL CONSTANT
		local selected		: list selected0 - notpen_o		//  selected now has only selected/penalized variables
		local s				: word count `selected'			//  number of selected/penalized vars EXCL NOTPEN/CONSTANT
		*

		*** error checks
		// check col vector of b (selected) vs. k
		local col_ct		=colsof(`b')
		if colsof(`b')==1 & el(`b',1,1)==. {	//  coef vector is empty so k should be zero
			if `k'~=0 {
				di as err "internal _lassopath error - r(k)=`k' does not match number of selected vars/coefs=`col_ct'"
				exit 499
			}
		}
		else if `k'~=`col_ct' {					//  coef vector not empty so k should match col count
			di as err "internal _lassopath error - r(k)=`k' does not match number of selected vars/coefs=`col_ct'"
			exit 499
		}
		// check col vector of bAll vs. p0
		local col_ct		=colsof(`bAll')
		if `p0'~=`col_ct' {						// full coef vector not empty so k should match col count
			di as err "internal _lassopath error - p0=`p0' does not match number of model vars/coefs=`col_ct'"
			exit 499
		}
		*

		*** return results
		// coef vectors are row vectors (Stata standard)
		return matrix beta			=`b'
		return matrix betaOLS		=`bOLS'
		return matrix betaAll		=`bAll'
		return matrix sbeta			=`sb'
		return matrix sbetaOLS		=`sbOLS'
		return matrix betaAllOLS	=`bAllOLS'
		return matrix sbetaAll		=`sbAll'
		return matrix sbetaAllOLS	=`sbAllOLS'
		return matrix Psi			=`Psi'
		//*//return matrix sPsi		=`sPsi'
		return matrix stdvec		=`stdvec'				//  penalized vars only
		return scalar lambda		=`lambda'				//  relates to depvar in original units
		return scalar slambda		=`slambda'				//  relates to standardized depvar
		//*//return scalar lambda0	=`lambda0'				//  BCH definition of lambda
		return scalar rmse		=`rmse'						//  lasso rmse
		return scalar rmseOLS	=`rmseOLS'					//  post-lasso rmse
		return scalar srmse		=`srmse'					//  standardized lasso rmse
		return scalar srmseOLS	=`srmseOLS'					//  standardized post-lasso rmse
		return scalar r2		=`r2'
		return scalar objfn		=`objfn'
		return scalar sobjfn	=`sobjfn'					//  standardized objective function
		return local  selected0	`selected0'					//  all selected/notpen vars INCLUDING NOTPEN (but excl constant)
		return local  selected	`selected' 					//  all selected (penalized) vars EXCLUDING NOTPEN & CONS
		return scalar k			=`k'						//  number of all vars in sel/notpen parameter vector INCLUDING CONSTANT
		return scalar s0		=`s0'						//  number of vars selected INCLUDING NOTPEN (but excl constant)
		return scalar s			=`s'						//  number of vars selected EXCLUDING NOTPEN & CONS
		return scalar p0		=`p0'						//  number of all vars in original model including constant
		return scalar N_clust	=`N_clust'
		return scalar N			=`N'
		//*//return scalar center	=`center'
		return local  clustvar	`cluster'
		return local  robust	`robust'
		return scalar niter		=`niter'
		return scalar maxiter	=`maxiter'
		//*//return scalar npsiiter	=`npsiiter'
		//*//return scalar maxpsiiter=`maxpsiiter'
		return scalar lcount 	= 1
		return scalar olsflag	= `olsflag'
		
		// added to lasso2
		return scalar df		= `df'
	}
	else {
		di as err "internal _lassopath error - lcount=`lcount'"
		exit 499
	}
end

// subroutine for checking/converting ploadings
program define getPsi, rclass
	version 13
	syntax ,					///
		[						///
		ploadings(string)		/// string with name of Stata matrix
		p(int 0)				/// size of model should = #dim loadings
		]

	// Check that ploadings is indeed a matrix
	tempname Psi
	cap mat `Psi' = `ploadings'
	if _rc != 0 {
		di as err "invalid matrix `ploadings' in ploadings option"
		exit _rc
	}
	// check dimension of col vector
	if (colsof(`Psi')!=(`p')) | (rowsof(`Psi')!=1) {
		di as err "`ploadings' should be a row vector of length `p'=dim(X) (excl constant)."
		exit 503
	}
	// check for negative loadings
	tempname hasneg
	mata: st_numscalar("`hasneg'",rowsum(st_matrix("`Psi'") :< 0))
	if `hasneg' {
		di as err "invalid ploadings matrix `ploadings' - has negative values"
		exit 508
	}

	return matrix Psi	= `Psi'

end

// subroutine for partialling out
program define _partial, rclass sortpreserve
	version 13
	syntax anything ,							/// anything is actually a varlist but this is more flexible - see below
		touse(varlist numeric max=1)			/// required `touse' variable (so marksample unneeded)
		[										/// 
		toest(varlist numeric max=1)			/// optional `toest' variable (subset on which standardization is based)
		PARtial(string)							/// string so that fv operators aren't inserted
		tvarlist(string)						/// optional list of temp model vars - may be "unfilled" (NaNs)
		wvar(string)							/// optional weight variable
		dmflag(int 0)							/// =0 if cons and not demeaned; =1 if no cons, already demeaned or treat as demeaned
		psolver(string)							/// 
		*										/// to catch extraneous or ignored options
		]

	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename.
	local varlist `anything'
	*

	*** toest vs touse
	// touse = all data
	// toest = estimation sample
	// if toest is missing, set equal to touse
	if "`toest'"=="" {
		tempvar toest
		qui gen byte `toest' = `touse'
	}
	*

	*** Error check - if tvarlist provided, should have same number of elements as varlist
	local v_ct	: word count `varlist'
	local tv_ct	: word count `tvarlist'
	if `tv_ct' & (`v_ct'~=`tv_ct') {
		di as err "internal lassoutils partialling error - mismatched lists"
		exit 198
	}
	*

	*** error check - allowed solvers are svd, qr, lu, luxx, chol
	if "`psolver'"~="" {
		local optsolver			svd svdxx qr qrxx lu luxx chol
		local pos				: list posof "`psolver'" in optsolver
		if `pos' == 0 {
			di as err "Syntax error: solver `psolver' not allowed"
			exit 198
		}
	}
	else {
		// default is QR+quadcross
		local psolver			qrxx
	}
	*
	
	*** partial out
	mata: s_partial("`varlist'", 			/// y and X
					"`partial'",          	/// to be partialled out
					"`tvarlist'",			///
					"`touse'",				/// touse
					"`toest'",				/// toest
					"`wvar'",				/// optional weight variable
					`dmflag', 				/// treatment of constant/demeaned or not
					"`psolver'")			//  choice of solver (optional)
	return scalar rank	=`r(rank)'			//  rank of matrix P of partialled-out variables 
	return local dlist	`r(dlist)'			//  list of dropped collinear variables in P
	*

end


// subroutine for fe transformation
program define _fe, rclass sortpreserve
	version 13
	syntax anything ,							/// anything is actually a varlist but this is more flexible - see below
		touse(varlist numeric max=1)			/// required `touse' variable (so marksample unneeded)
		tvarlist(string)						/// optional list of temp vars - may be "unfilled" (NaNs)
		FE(varlist numeric min=1 max=1) 		/// fe argument is ivar
		[										///
		wvar(varlist numeric max=1)				/// optional weight variable
		toest(varlist numeric max=1)			/// optional `toest' variable (subset on which standardization is based)
		NOFTOOLS								///
		]

	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename.
	local varlist `anything'
	*

	*** toest vs touse
	// toest = estimation sample
	// touse = all data
	// if toest is missing, set equal to touse
	if ("`toest'"=="") {
		tempvar toest
		qui gen `toest'=`touse'
	}
	*

	*** ftools
	// do not use ftools if (a) specifcally not requested; (b) not installed
	// use "which ftools" to check - faster, plus compile was already triggered
	// in conditional load section at end of this ado
	if "`noftools'"=="" {
		cap which ftools
		if _rc {
			// fails check, likely not installed, so use (slower) Stata code
			local noftools "noftools"
		}
	}
	*
	
	*** error check - weights support requires ftools
	if "`noftools'"~="" & "`wvar'"~="" {
		di as err "error - fe option with weights requires " _c
		di as err `"{rnethelp "http://fmwww.bc.edu/RePEc/bocode/f/ftools.sthlp":ftools} package to be installed"'
		di as err `"see {rnethelp "http://fmwww.bc.edu/RePEc/bocode/f/ftools.sthlp":help ftools} for details"'
		exit 198
	}
	*
	
	*** Error check - if tvarlist provided, should have same number of elements as varlist
	local v_ct	: word count `varlist'
	local tv_ct	: word count `tvarlist'
	if (`v_ct'~=`tv_ct') {
		di as err "internal lassoutils FE error - mismatched lists"
		exit 198
	}
	*

	if "`noftools'"~="" {
		// timer on 1
		*** Within-transformation / demeaning
		// varlist should be all doubles so no recast needed
		// xtset is required for FEs so this check should never fail
		cap xtset
		if _rc {
			di as err "internal lassoutils xtset error"
			exit 499
		}
		// panelvar always exists; timevar may be empty
		// if data xtset by panelvar only, may not be sorted
		// if data xtset by panelvar and time var, will be sorted by both
		// resorting can be triggers by simple call to xtset
		local panelvar `r(panelvar)'
		local timevar `r(timevar)'
		// sort by panel and estimation sample; put latest observation last
		sort `panelvar' `toest' `timevar'
		// toest_m is indicator that this ob will have the mean
		// N_est is number of obs in panel and in estimation sample
		tempvar toest_m N_est
		// last ob in panel and estimation sample tagged to have mean
		qui by `panelvar' `toest' : gen `toest_m' = (_n==_N) & `toest'
		qui by `panelvar' `toest' : gen `N_est' = sum(`toest')
		// count is of panels used in estimation
		qui count if `toest_m'
		local N_g	=r(N)
		// if xtset by panel and time vars, restore sort
		cap xtset
		// create means for each variable
		foreach var of varlist `varlist' {
			tempvar `var'_m
			local mlist `mlist' ``var'_m'
			// use only training/estimation data to calculate mean (toest)
			qui by `panelvar' : gen double ``var'_m'=sum(`var')/`N_est' if `toest'
			qui by `panelvar' : replace ``var'_m' = . if ~`toest_m'
		}
		// sort so that last ob in each panel has the mean in it
		sort `panelvar' `toest_m'
		// and propagate to the rest of the panel
		foreach var of varlist `mlist' {
			qui by `panelvar' : replace `var' = `var'[_N] if _n<_N
		}
		// if xtset by panel and time vars, restore sort
		// need to do this if e.g. any time-series operators in use
		cap xtset
		// finally, demean data
		local i 1
		foreach var of varlist `varlist' {
			local tvar	: word `i' of `tvarlist'
			qui replace `tvar'=`var'-``var'_m' if `touse'
			local ++i
		}
		return scalar N_g = `N_g'
		*
		// timer off 1
	}
	else {
		// Mata routine; uses Sergio Correia's FTOOLS package.
		// timer on 2
		mata: s_fe("`varlist'", 				/// 
					"`tvarlist'",				///
					"`wvar'",					///
					"`fe'",         		 	///
					"`touse'",					///
					"`toest'")
		local N_g	=r(N_g)
		return scalar N_g = `N_g'
		// timer off 2
	}
	// indicate whether FTOOLS not used
	return local noftools `noftools'
end


// subroutine for standardizing in Stata
program define _std, rclass
	version 13
	syntax anything ,							/// anything is actually a varlist but this is more flexible - see below
		touse(varlist numeric max=1)			/// required `touse' variable (so marksample unneeded)
		[										///
		toest(varlist numeric max=1)			/// optional `toest' variable (subset on which standardization is based)
		tvarlist(string)						/// optional list of temp vars - may be "unfilled" (NaNs)
		consmodel(int 1)						/// =1 if constant in model (=> data are to be demeaned)
		dmflag(int 0)							/// =1 if data already demeaned or treat as demeaned
		NOChange								/// don't transform the data; just return the std vector
		]

	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename.
	local varlist `anything'
	*
	
	local transform = ("`nochange'"=="")
	
	*** toest vs touse
	// touse = all data
	// toest = estimation sample
	// if toest is missing, set equal to touse
	if "`toest'"=="" {
		tempvar toest
		qui gen byte `toest' = `touse'
	}
	*

	tempname stdvec mvec
	mata: s_std("`varlist'","`tvarlist'","`touse'","`toest'",`consmodel',`dmflag',`transform') 
	mat `stdvec' = r(stdvec)
	mat `mvec' = r(mvec)
	// these will be tempnames if std program called with tempvars
	mat colnames `stdvec' = `varlist'
	mat colnames `mvec' = `varlist'
	return matrix stdvec = `stdvec'
	return matrix mvec = `mvec'	
end


// subroutine for weighting in Stata
program define _wt, rclass
	version 13
	syntax anything ,							/// anything is actually a varlist but this is more flexible - see below
		touse(varlist numeric max=1)			/// required `touse' variable (so marksample unneeded)
		wvar(varlist numeric max=1)				/// required `wvar' variable
		tvarlist(string)						/// required list of temp vars - may be "unfilled" (NaNs)
		[										///
		NOChange ///
		]

	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename.
	local varlist `anything'
	*
	
	local transform = ("`nochange'"=="")
	
	mata: s_wt("`varlist'","`tvarlist'","`touse'","`wvar'",`transform')
	
end


// subroutine for recovering coefs of partialled-out vars
program define _unpartial, rclass sortpreserve
	version 13
	syntax anything ,							/// anything is actually a varlist but this is more flexible - see below
		touse(varlist numeric max=1)			/// required `touse' variable (so marksample unneeded)
		[										///
		beta(string)							///
		depvar(string)							///
		scorevars(string)						/// names of selected vars (in beta)
		wvar(string)							///
		names_t(string)							/// string so that fv operators aren't inserted
		names_o(string)							/// ditto
		PARtial(string)							/// ditto
		consmodel(int 1)						/// include constant when recovering coefs
		]

	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename.
	local varlist `anything'
	*

	local depvar `varlist'
	tempname b bpartial
	tempvar xb yminus

	local scorevar_ct	: word count `scorevars'
	local partial_ct	: word count `partial'
	
	// betaempty = 1 if not supplied
	local betaempty		= `scorevar_ct'==0

	if `scorevar_ct' {
		mat `b' = `beta'
		mat colnames `b' = `scorevars'
		qui mat score double `xb' = `b' if `touse'
		qui gen double `yminus' = `depvar' - `xb' if `touse'
	}
	else {
		qui gen double `yminus' = `depvar' if `touse'
	}

	// if beta has missings, yminus all missing, cannot unpartial
	qui count if `yminus'==. & `touse'
	local betamissing	= `r(N)' > 0
	if `betamissing' {
		// beta has missings so bpartial needs to be all missing
		// and must create b as all missings
		mat `b'				= J(1,`scorevar_ct',.)
		mat colnames `b'	= `scorevars'
		if `consmodel' & `partial_ct'==0 {
			// only constant was partialled-out
			mat `bpartial' = J(1,1,.)
			mat colnames `bpartial' = _cons
		}
		else if `partial_ct' & `consmodel'==0 {
			// partialled-out vars but no constant
			mat `bpartial' = J(1,`partial_ct',.)
			mat colnames `bpartial' = `partial'
		}
		else if `partial_ct' & `consmodel' {
			// both constant and partialled-out vars
			mat `bpartial' = J(1,`partial_ct'+1,.)
			mat colnames `bpartial' = `partial' _cons
		}
		else {
			di as err "internal lassoutils error - unpartialling"
			exit 499
		}
	}
	else {
		// unpartial
		// weights; note that since we use regress+weights, all vars here are unweighted
		if "`wvar'"~="" {
			local wexp [aw=`wvar']
		}
		//  partial uses _tnames; use _t names and then convert to _o names
		//  _t names will be FE-transformed or just the values of the original vars (without i/b/n etc.)
		if ~`consmodel' {
			local noconstant noconstant
		}
		qui reg `yminus' `partial' if `touse' `wexp', `noconstant'
		mat `bpartial' = e(b)
	}

	// replace temp names with original names
	local cnames_t	: colnames `bpartial'				//  may have _cons at the end already
	fvstrip `cnames_t'									//  regress output may have omitted vars
	local cnames_t	`r(varlist)'
	matchnames "`cnames_t'" "`names_t'" "`names_o'"
	local cnames_o	`r(names)'
	mat colnames `bpartial' = `cnames_o'
	// may be omitteds so need to reinsert o. notation
	// trick is to use _ms_findomitted utility but give it
	// diag(bAll) as vcv matrix where it looks for zeros
	tempname tempvmat
	mat `tempvmat'	= diag(`bpartial')
	_ms_findomitted `bpartial' `tempvmat'
	// build in fv info
	_ms_build_info `bpartial' if `touse'
	// attach to b matrix if b not empty
	if `betaempty' == 0 {
		mat `b' = `beta' , `bpartial'
	}
	else {
		mat `b' = `bpartial'
	}

	return matrix b			= `b'
	return matrix bpartial	= `bpartial'

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

********************************************************************************
*** Mata section															 ***
********************************************************************************

version 13
mata:

// data structure
struct dataStruct {
	pointer colvector y
	pointer matrix X
	string colvector nameX		// names of actual variables (can be tempvars)
	string colvector nameX_o	// original names of variables
	real scalar cons			// =1 if model also has a constant, =0 otherwise
	real scalar dmflag			// =1 if data have mean zero or should be treated as demeaned, =0 otherwise
	real scalar n				// number of observations
	real scalar p				// number of columns in X (may included constant)
	real rowvector sdvec		// vector of SDs of X
	real rowvector varvec		// vector of variances of X
	real scalar ysd				// standard deviation of y
	real scalar prestdflag		// flag indicating data have been pre-standardized
	real rowvector prestdx		// prestandaridization vector of SDs for X
	real scalar prestdy			// prestandaridization SD for y
	real rowvector mvec			// vector of means
	real scalar ymvec			// mean of y (not actually a vector)
	real scalar lglmnet			// =1 if glmnet lambda/alpha scaling used, =0 otherwise
	real matrix XX				// cross prod of all Xs
	real matrix Xy				// cross prod of all Xs and y
	real scalar TSS				// total sum of squares, used in various places
	real scalar dofminus		// degrees of freedom used by FEs
	real scalar sdofminus		// degrees of freedom used by partialled-out vars ("s"=small; includes cons)
	real scalar sqrtflag		// =1 if sqrt-lasso, =0 otherwise
// below used only by rlasso
	real matrix pihat			// used for partialling out with rlasso
	real matrix ypihat			// used for partialling out with rlasso
	pointer matrix Xp			// penalized Xs
	pointer matrix Xnp			// unpenalized Xs
	real scalar np				// number of unpenalized Xs; also used as flag
	string colvector nameXp		// names of actual variables (can be tempvars)
	string colvector nameXnp	// names of actual variables (can be tempvars)
	real scalar hetero			// =1 if het-robust penalty loadings, =0 otherwise
	real scalar center			// center x_i*e_i or x_ij*e_ij
	pointer colvector clustid	// cluster id
	pointer colvector clustid1	// cluster id1
	pointer colvector clustid2	// cluster id2
	pointer colvector clustid3	// cluster id3
	real scalar nclust			// number of clusters; 0 if no clustering
	real scalar nclust1			// number of clusters (2-way); 0 if no clustering
	real scalar nclust2			// number of clusters (2-way); 0 if no clustering
	real scalar nclust3			// number of clusters (intersection); 0 if no clustering
	real scalar nclust1flag		// use #nclust-1 instead of #nclust in cluster-lasso
	real scalar bw				// bandwidth (for HAC/AC); =0 implies not HAC/AC
	string scalar kernel		// kernel (for HAC/AC)
	real scalar maqflag			// specific to truncated kernel
	real scalar spectral		// =1 if spectral (use all all lag), =0 otherwise
	string scalar tindexname	// name of actual variable (will be a tempvar)
	real scalar tdelta			// tdelta for TS operators
	real scalar bsize			// block size for the HAC/AC multiplier bootstrap
	real rowvector selXp		// selection row vector for penalized vars
	real rowvector selXnp		// selection row vector for unpenalized vars
	real rowvector selindXp		// selection index row vector for penalized vars
	real rowvector selindXnp	// selection index row vector for unpenalized vars
	real rowvector mvecp		// vector of means of penalized Xs
	real rowvector mvecnp		// vector of means of unpenalized Xs
	real rowvector sdvecp		// sdvec of penalized vars after partialling-out unpenalized vars
	real rowvector sdvecpnp		// same as sdvecp except conformable with full set of Xs (p and np)
	real scalar ysdp			// SD of y after partialling out unpenalized vars
	real scalar psinegs			// number of times negative penalty loadings encountered
	string scalar psinegvars	// names of variables with neg penalty loadings
	}

struct outputStruct {
	real colvector beta			// beta in original units
	//real colvector sbeta		// beta in standardized units
	real colvector betaPL		// OLS (post-lasso) beta in std units 
	//real colvector sbetaPL
	real colvector betaAll		// full beta vector
	real colvector betaAllPL	// full OLS (post-lasso) beta vector
	real colvector beta_init
	real scalar cons			// flag =1 if cons in regression, =0 if no cons
	real scalar intercept		// estimated intercept
	real scalar interceptPL		// estimated OLS (post-lasso) intercept
	string colvector nameXSel	// vector of names of selected Xs
	real colvector index		// index of selected vars
	real rowvector Psi			// penalty loadings
	real rowvector sPsi			// standardized penalty loadings
	real rowvector ePsi			// estimated penalty loadings (rlasso only)
	real scalar prestdflag		// flag indicating data have been pre-standardized
	real colvector v			// residuals
	real colvector vPL			// OLS (post-lasso) residuals
	real scalar lambda			// penalty scalar
	real scalar slambda			// standardized
	real scalar lambda0			// lambda without estimate of sigma (rlasso only)
	real scalar c				// part of BCH lambda
	real scalar gamma			// part of BCH lambda
	real scalar gammad			// part of BCH lambda
	real scalar rmse			// rmse using estimated beta
	real scalar rmsePL			// rmse using OLS (post-lasso) beta
	real scalar r2				// r-sq using estimated beta
	real scalar n				// number of obs
	real scalar s				// number of selected vars
	real scalar nclust			// number of clusters (rlasso only)
	real scalar nclust1			// number of clusters (dimension 1)
	real scalar nclust2			// number of clusters (dimension 2)
	real scalar niter			// number of iterations
	real scalar npsiiter		// number of rlasso penalty loadings iterations
	real rowvector supscore		// sup-score stats: BCH stat, p-value, crit-value, signif; rlasso stat, p-value
	real scalar ssgamma			// gamma for sup-score conservative critical value
	real scalar objfn			// value of minimized objective function
	real scalar dof				// "effective" degrees of freedom
	}
// end outputStruct

struct outputStructPath {
	real colvector lambdalist0	// original (untruncated) lambdalist
	real colvector lambdalist
	real matrix betas
	real matrix sbetas
	real rowvector Psi
	real rowvector sdvec
	real rowvector sPsi
	real colvector dof			// "effective" degrees of freedom
	real colvector shat			// number of non-zero parameters
	real colvector shat0		// number of nonzero params including partialled-out/cons
	real scalar cons			// yes/no
	real rowvector intercept
	real rowvector interceptPL
	real scalar n
	real scalar nclust
	}
// end outputStructPath

struct betaStruct {
	real colvector beta			// estimated beta
	real scalar m				// number of iterations
	real scalar fobj			// value of objective function
	real scalar dev				// deviation (R-squared)
}

struct dataStruct scalar MakeData(	string scalar nameY,		/// #1
									string scalar nameX,		/// #2
									string scalar nameX_o,		/// #3
									string scalar touse,		/// #4
									real scalar cons,			/// #5
									real scalar dmflag,			/// #6
									real scalar dofminus,		/// #7  lost df from FEs
									real scalar sdofminus,		/// #8  lost df from partialling
									real scalar prestdflag,		/// #9
									real scalar sqrtflag,		/// #10
									string scalar stdymat,		/// #11
									string scalar stdxmat,		/// #12
									real scalar lglmnet,		/// #13
									|							/// optional arguments
									string scalar nameclustid,	/// #14 optional arguments - rlasso-specific
									real scalar hetero,			/// #15 =1 if het-robust penalty loadings, =0 otherwise
									real scalar center,			/// #16 center x_i*e_i or x_ij*e_ij in cluster-lasso
									real scalar nclust1flag,	/// #17 use #nclust-1 instead of #nclust in cluster-lasso
									string scalar nameP,		/// #18
									string scalar nameNP,		/// #19
									real scalar bw,				/// #20
									string scalar kernel,		/// #21
									real scalar maqflag,		/// #22
									real scalar spectral,		/// #23
									string scalar tindexname,	/// #24
									real scalar tdelta,			/// #25
									real scalar bsize,			/// #26
									string scalar psolver		/// #27
									)
{

	if (args()<12) {
		stdymat		= ""
		stdxmat		= ""
		nameclustid	= ""
		hetero		= 0
		center		= 0
		nclust1flag	= 0
	}
	if (args()<18) {
		nameP		= ""			//  default list is empty
		nameNP		= ""			//  default list is empty
	}
	if (args()<20) {
		bw			= 0
		kernel		= ""
		maqflag		= 0
		spectral	= 0
		tindexname	= ""			//  no time index provided
		tdelta		= 0
		bsize		= 0
		psolver		= ""
	}

	struct dataStruct scalar d
	

	// dep var
	st_view(y,.,nameY,touse)
	d.y			=&y

	// X vars
	st_view(X,.,nameX,touse)
	d.X			=&X
	d.nameX		=tokens(nameX)
	d.nameX_o	=tokens(nameX_o)
	d.n			=rows(X)
	d.p			=cols(X)
	d.cons		=cons			//  model has constant to be recovered by estimation code; also used in variable counts
	d.dmflag	=dmflag			//  treat data as zero mean
	d.dofminus	=dofminus
	d.sdofminus	=sdofminus
	d.prestdflag=prestdflag

	// sqrtflag
	d.sqrtflag	=sqrtflag
	
	// scaling
	d.lglmnet	=lglmnet
	
	// rlasso het-rob etc.
	d.hetero	=hetero
	d.center	=center

	// clustering
	clustertokens	= tokens(nameclustid)
	if (cols(clustertokens)==1) {
		// one-way clustering
		st_view(cid1,.,nameclustid,touse)
		d.clustid1	= &cid1
		info		= panelsetup(cid1, 1)
		d.nclust1	= rows(info)
		d.nclust2	= 0
		d.nclust3	= 0
		d.nclust	= d.nclust1
	}
	else if (cols(clustertokens)>1) {
		// two-way clustering
		// dimension 3 is insersection of dim 1 and dim 2
		// data have been sorted on dim 3 and dim 1 (so ordering preserved for both since 1 nests 3)
		// dimension 1:
		st_view(cid1,.,clustertokens[1,1],touse)
		d.clustid1	= &cid1
		info		= panelsetup(cid1, 1)
		d.nclust1	= rows(info)
		// dimension 2 (can't use panelsetup since not sorted on dim 2):
		st_view(cid2,.,clustertokens[1,2],touse)
		d.clustid2	= &cid2
		d.nclust2	= rows(uniqrows(cid2))
		// dimension 3:
		st_view(cid3,.,clustertokens[1,3],touse)
		d.clustid3	= &cid3
		info		= panelsetup(cid3, 1)
		d.nclust3	= rows(info)
		d.nclust	= min((d.nclust1,d.nclust2))
	}
	else {
		d.nclust	= 0
	}
	d.nclust1flag	= nclust1flag
	
	// time-series HAC/AC
	d.bw			= bw
	d.kernel		= kernel
	d.maqflag		= maqflag
	d.spectral		= spectral
	d.tindexname	= tindexname
	d.tdelta		= tdelta
	d.bsize			= bsize
	
	// misc
	d.psinegs		= 0
	d.psinegvars	= ""

	// mean vectors
	if (dmflag) {
		// already demeaned or treat as demeaned
		d.mvec		= J(1,d.p,0)
		d.ymvec		= 0
	}
	else {
		// calculate means
		d.mvec		= mean(*d.X)
		d.ymvec		= mean(*d.y)
	}

	// SDs used to prestandardize
	if (prestdflag) {
		// data are prestandardized so just store these
		d.prestdy	=st_matrix(stdymat)
		d.prestdx	=st_matrix(stdxmat)
	}
	else {
		// unit vectors
		d.prestdy	=1
		d.prestdx	=J(1,d.p,1)
	}

	// standardization vectors
	// nb: d.ysd may be unused in code
	if (prestdflag) {
		// data are prestandardized so all unit vectors
		d.ysd		=1
		d.sdvec		=J(1,d.p,1)
		d.varvec	=J(1,d.p,1)
	}
	else if (dmflag) {
		// already demeaned (mean zero) or treat as demeaned
		d.ysd		= sqrt(mean((*d.y):^2))
		d.varvec	= mean((*d.X):^2)
		d.sdvec		= sqrt(d.varvec)
	}
	else {
		// not mean zero so need to demean
		d.ysd		= sqrt(mean(((*d.y):-d.ymvec):^2))
		d.varvec	= mean(((*d.X):-d.mvec):^2)
		d.sdvec		= sqrt(d.varvec)
	}

	if (cons) {
		// unpenalized constant present so X'X is in mean-devation form
		d.XX	= quadcrossdev(*d.X,d.mvec,*d.X,d.mvec)
		d.Xy	= quadcrossdev(*d.X,d.mvec,*d.y,d.ymvec)
		d.TSS	= quadcrossdev(*d.y,d.ymvec,*d.y,d.ymvec)
	}	
	else {
		// either data are already zero-mean or no constant is present in the model
		d.XX	= quadcross(*d.X,*d.X)
		d.Xy	= quadcross(*d.X,*d.y)
		d.TSS	= quadcross(*d.y,*d.y)
	}

	// rlasso section
	// Note that rlasso code requires uses SDs for additional purposes, and the SDs are after partialling out.  
	// These are saved as d.ysdp, d.sdvecp, d.sdvecpnp (inserted into a vector conformable with X).
	// sdvecpnp has the SDs of the Xs after partialling out + 0s for the SDs of the partialled vars
	if (nameNP!="") {
		// unpenalized regressors in rlasso
		st_view(Xp,.,nameP,touse)
		st_view(Xnp,.,nameNP,touse)
		d.Xp		=&Xp
		d.Xnp		=&Xnp
		d.np		=cols(Xnp)
		d.nameXp	=tokens(nameP)
		d.nameXnp	=tokens(nameNP)
		selXp		= J(1,cols(*d.X),1)
		forbound	= cols(d.nameXnp)		//  faster
		for (i=1;i<=forbound;i++) {
			selXp = selXp - (d.nameX :== d.nameXnp[1,i])
		}
		d.selXp		=selXp
		d.selXnp	=1:-selXp
		d.selindXp	=selectindex(d.selXp)
		d.selindXnp	=selectindex(d.selXnp)

		if (cons) {
			// standard case - model has a constant, data are not demeaned
			// model has a constant (unpenalized) so cross-prods etc are in mean-dev form
			d.mvecp		= mean(*d.Xp)
			d.mvecnp	= mean(*d.Xnp)
			d.ypihat	= anysolver((*d.Xnp):-d.mvecnp,(*d.y):-d.ymvec,r=.,psolver)
			d.pihat		= anysolver((*d.Xnp):-d.mvecnp,(*d.Xp):-d.mvecp,r=.,psolver)
			d.ysdp		= sqrt(mean((((*d.y):-d.ymvec)-((*d.Xnp):-mean(*d.Xnp))*d.ypihat):^2))
			// std vector for just the Xp variables.
			d.sdvecp	= sqrt(mean((((*d.Xp):-d.mvecp)-((*d.Xnp):-mean(*d.Xnp))*d.pihat):^2))
		}
		else if (dmflag) {
			// zero-mean data and no constant is present in the model
			// means are zero vectors, cross-prods don't need demeaning, SDs don't need demeaning
			d.mvecp		= J(1,cols(*d.Xp),0)
			d.mvecnp	= J(1,cols(*d.Xnp),0)
			d.ypihat	= anysolver(*d.Xnp,*d.y,r=.,psolver)
			d.pihat		= anysolver(*d.Xnp,*d.Xp,r=.,psolver)
			d.ysdp		= sqrt(mean(((*d.y)-(*d.Xnp)*d.ypihat):^2))
			// std vector for just the Xp variables.
			d.sdvecp	= sqrt(mean(((*d.Xp)-(*d.Xnp)*d.pihat):^2))
		}
		else {
			// model has no constant but means may be nonzero
			// hence cross-prods don't need demeaning but SDs do
			d.mvecp		= mean(*d.Xp)
			d.mvecnp	= mean(*d.Xnp)
			d.ypihat	= anysolver(*d.Xnp,*d.y,r=.,psolver)
			d.pihat		= anysolver(*d.Xnp,*d.Xp,r=.,psolver)
			d.ysdp		= sqrt(	mean(																			///
								(centerpartial(d.y,d.ymvec,d.cons,d.Xnp,d.mvecnp,d.ypihat)						///
									:- mean(centerpartial(d.y,d.ymvec,d.cons,d.Xnp,d.mvecnp,d.ypihat))):^2		///
								 ) )
			// std vector for just the Xp variables.
			d.sdvecp	= sqrt( mean(																			///
								(centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat)						///
									:- mean(centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat))):^2		///
								) )
		}

		// Now create blank full-size std vector of zeros for all X variables.
		d.sdvecpnp					= J(1,cols(*d.X),0)
		// Then insert values in appropriate columns.
		d.sdvecpnp[1,(d.selindXp)]	= d.sdvecp
	}
	else {													//  no unpenalized regressors (in rlasso)
		d.Xp		= d.X									//  penalized = all
		d.nameXp	= d.nameX
		d.selXp		= J(1,cols(*d.X),1)
		d.selindXp	= selectindex(d.selXp)
		d.np		= 0
		d.pihat		= 0
		d.ypihat	= 0
		d.mvecp		= d.mvec
		d.ysdp		= d.ysd									//  can just reuse these
		d.sdvecp	= d.sdvec
		d.sdvecpnp	= d.sdvecp								//  since all penalized, pnp=p (penalized-not penalized = penalized)

	}

	return(d)
}
// end MakeData



// this function calculates lasso path for range of lambda values
struct outputStructPath DoLassoPath(									///
									struct dataStruct scalar d,			///
									real rowvector Psi,					/// vector of penalty loadings (L1 norm)
									real rowvector Psi2,				/// vector of penalty loadings (L1 norm)
									real rowvector lvec,				/// vector of lambdas (L1 norm)
									real rowvector lvec2,				/// vector of lambdas (L1 norm)
									real scalar post,					/// 
									real scalar verbose,				/// 
									real scalar optTol,					/// 
									real scalar maxIter,				/// 
									real scalar zeroTol,				/// 
									real scalar fdev,					/// (glmnet) minimum fractional change in deviance for stopping path
									real scalar devmax,					/// (glmnet) maximum fraction of explained deviance for stopping path
									real scalar alpha,					/// 
									real scalar lglmnet,				/// 
									real scalar nodevcrit)
{

		struct outputStructPath scalar	t
		struct betaStruct scalar		b

		if (verbose>=2) {
			printf("Lambda list: %s\n",invtokens(strofreal(lvec)))
		}
		
		lmax	= max(lvec)
		lcount	= cols(lvec)

		// initialize deviance (R-sq)
		dev		= 1
		lpath	= J(lcount,d.p,.) // create empty matrix which stores coef path
		if (verbose>=1) {
			printf("Iterating over lambda:")
		}
		for (k = 1;k<=lcount;k++) { // loop over lambda

			lambda	= lvec[1,k]
			lambda2	= lvec2[1,k]

			if (verbose>=1) {
				printf(" %f",lambda)
			}
			
			// set verbose to zero (no output for each lambda)
			if (k==1) {
				// DoShooting(.) handles default initial beta
				b			= DoShooting(d,Psi,Psi2,lambda,lambda2,0,optTol,maxIter,zeroTol,alpha)
			}
			else {
				// initial beta is previous beta on the path
				b			= DoShooting(d,Psi,Psi2,lambda,lambda2,0,optTol,maxIter,zeroTol,alpha,beta)
			}
			beta		= b.beta
			lpath[k,.]	= beta'
			m			= b.m
			fobj		= b.fobj
			dev0		= dev
			dev			= b.dev
			lcount1		= k					// count of betas/lambdas estimated
			if (nodevcrit==0) {												// check dev criterion for exit
				if ((abs(dev-dev0) < fdev) & (dev>0))		break			// don't exit if dev=R-sq=0
				if (dev >= devmax)							break
			}
		}
		if (verbose>=1) {
			printf("\n")
		}

		// small coefs set to zero IF NOT RIDGE
		// compares to STANDARDIZED coefs
		if (alpha~=0) {
			if (d.prestdflag) {
				lpath=edittozerotol(lpath, zeroTol)
			}
			else {
				// compare to standardized coefs
				lpath0	= edittozerotol(lpath :* d.sdvec / d.ysd, zeroTol) :/ d.sdvec * d.ysd
				// it is possible for the SD of a predictor to be zero, in which case lpath0 will have missings in it
				// so we set lpath to have zeros or missings where they appear in lpath0
				lpath	= lpath :* ( (lpath0 :~= 0) :& (lpath0 :~= .) )
			}
		}

		if (post) { 
			betasPL = J(lcount1,d.p,0)
			nonzero0 = J(1,d.p,0)
			for (k = 1;k<=lcount1;k++) {				// loop over lambda points estimated (lcount1 <= lcount)
				nonzero = lpath[k,.]:!=0				// 0-1 vector
				sk = sum(nonzero)			
				if ((0<sk) & (sk<d.n)) {				// only if 0<s<n
					if ((nonzero0==nonzero) & (k>=2)) {	// no change in active set
						betasPL[k,.] = betasPL[k-1,.]
					}
					else { 
						ix = selectindex(nonzero)		// index of non-zeros
						// obtain beta-hat
						if (d.cons==0) {
							// data are mean zero or there is no constant in the model
							betak=qrsolve(select((*d.X),nonzero),*d.y)
						}
						else {
							betak=qrsolve(select((*d.X),nonzero):-select(d.mvec,nonzero),((*d.y):-d.ymvec))
						}
						betasPL[k,ix] = betak'
						nonzero0=nonzero
					}
				}
			}
			t.betas 		= betasPL	
		}
		else {
			t.betas			= lpath[1..lcount1,.]
		}
		
		// insert intercept
		// no intercept if pre-standardized
		if (t.cons) {
			t.intercept		= mean(*d.y) :- mean(*d.X)*t.betas'
		}

		t.lambdalist0	= lvec'
		t.lambdalist	= lvec[1..lcount1]'
		t.Psi			= Psi
		t.cons 			= d.cons
		
		// degrees of freedom and dimension of the model (add constant)
		// either count all partialled out or none
		t.shat	= rowsum(t.betas:!=0)
		t.shat0	= rowsum(t.betas:!=0) :+ (d.sdofminus)

		// requires only lvec2 and Psi2 (for ridge/elastic net)
		t.dof		= getdf(d,t.betas',Psi2,lvec2,alpha,verbose)'
				
		return(t)		

}
// end DoLassoPath

// calculate effective degrees of freedom
// betas is a colvector => return scalar
// betas is a matrix (each col a beta vector) => return a rowvector
function getdf(											///
					struct dataStruct scalar d,			///
					real matrix betas,					/// can be a matrix of multiple betas; single beta is a COL vector
					real rowvector Psi2,				/// vector of penalty loadings (L2 norm)
					real rowvector lvec2,				/// vector of lambdas (L2 norm)
					real scalar alpha,					/// enet parameter
					real scalar verbose					/// reporting
					)

{

			// initialize
			df					= J(1,cols(betas),.)
			// not all lambdas used in estimation
			beta_ct				= cols(betas)

			if (alpha==1) { 
				// lasso dof
				df[.,.]			= colsum(betas:!=0)
			}
			else if (alpha==0) {
				// ridge dof  
				if (d.cons)		Xt	= (*d.X) :- d.mvec
				else			Xt	= (*d.X)

				for (k=1;k<=beta_ct;k++) { // loop over lambda points
					if (d.lglmnet) {
						df[1,k]		= trace((Xt)*invsym(quadcross(Xt,Xt):+d.n*lvec2[1,k]*diag(Psi2:^2))*(Xt)')
					}
					else {
						df[1,k]		= trace((Xt)*invsym(quadcross(Xt,Xt):+lvec2[1,k]*d.prestdy/2*diag(Psi2:^2))*(Xt)')
					}
				}
			}
			else {
				// elastic net dof
				if (d.cons)		Xt	= (*d.X) :- d.mvec
				else			Xt	= (*d.X)
				
				for (k=1;k<=beta_ct;k++) {			// loop over lambda points
					nonzero		= betas[.,k]:!=0	// 0-1 col vector corresponding to nonzero beta elements
					XA			= select((Xt),nonzero')
					Psi2A		= select((Psi2),nonzero')
					if (d.lglmnet) {
						df[1,k]		= trace((XA)*invsym(quadcross(XA,XA):+(1-alpha)*d.n*lvec2[1,k]*diag(Psi2A:^2))*(XA)')
					}
					else {
						df[1,k]		= trace((XA)*invsym(quadcross(XA,XA):+(1-alpha)*lvec2[1,k]*d.prestdy/2*diag(Psi2A:^2))*(XA)')
					}
				}
			}
			// need to add the constant / account for demeaning / partialled-out vars ("s" for "small")
			df				= df :+ d.sdofminus
			// dof cannot be greater than n
			df_trunc		= colmin( (df \ J(1,cols(df),d.n)) )
			if (df ~= df_trunc) {
				df			= df_trunc
				if ((verbose>=1) & (beta_ct==1)) {
					printf("{txt}Warning: effective degrees of freedom > sample size; truncating.\n")
				}
				else if (verbose>=1) {
					printf("{txt}Warning: effective degrees of freedom > sample size in some cases; truncating.\n")
				}

			}
		
			return(df)
}							

struct betaStruct DoShooting(									///
							struct dataStruct scalar d,			/// #1 
							real rowvector Psi,					/// #2  vector of penalty loadings (L1 norm)
							real rowvector Psi2,				/// #3  vector of penalty loadings (L2 norm)
							real scalar lambda,					/// #4  lambda (single value) (L1 norm)
							real scalar lambda2,				/// #5  lambda (single value) (L2 norm)
							real scalar verbose,				/// #6 
							real scalar optTol,					/// #7 
							real scalar maxIter,				/// #8 
							real scalar zeroTol,				/// #9 
							real scalar alpha,					/// #10
							|									///
							real colvector beta_init			/// #11 initial beta vector (optional)
							)

{

	// no initial beta supplied
	noinitflag		= args()<11

	// stores results
	struct betaStruct scalar	b

	// sqrt elastic net not supported - must be sqrt lasso
	if (d.sqrtflag)	alpha	= 1

	if (verbose>=1) {
		avgPsi=sum(abs(Psi))/sum(Psi:>0)
		printf("{txt}Lambda: {res}%f\n{txt}Average abs. loadings: {res}%f\n", lambda, avgPsi)
	}

	// initialize cross-product matrices
	if (d.lglmnet) {
		XX		= d.XX/d.n
		Xy		= d.Xy/d.n
	}
	else if (d.sqrtflag) {
		XX		= d.XX/d.n
		Xy		= d.Xy/d.n
	}
	else {
		XX		= (d.XX)*2
		Xy		= (d.Xy)*2
	}

	// ridge; also required if no initial beta matrix supplied
	if ((alpha==0) | (noinitflag)) {
		
		if (d.lglmnet) {
			beta_ridge	= anysolver(XX+lambda2*diag(Psi2),Xy, r=.)		// r=rank; default LU, use QR if rank-deficient
		}
		else if (d.sqrtflag) {
			beta_ridge	= anysolver(XX*d.n*2+lambda2*d.prestdy*diag(Psi2:^2),Xy*d.n*2, r=.)
		}
		else {
			beta_ridge	= anysolver(XX+lambda2*d.prestdy*diag(Psi2:^2),Xy, r=.)
		}
		if ((verbose>=1) & (r<cols(XX))) {
			printf("{txt}Note: collinearities encountered in obtaining ridge solution.\n")
		}
		if (verbose>=2) {
			printf("{txt}Ridge solution:\n")
			invtokens(strofreal(beta_ridge'))
		}
	}	// end ridge
	
	if (alpha==0) {
		// ridge - closed form solution above
		beta	= beta_ridge
		// set m so that #iterations is returned as nonmissing
		m		= 0
	}
	else if (d.sqrtflag==0) {
		// lasso and elastic net - coordinate descent

		// initialize beta with ridge if beta_init not provided
		if (noinitflag) {
			beta_init	= beta_ridge
			if (verbose>=1) {
				printf("{txt}Initial beta vector is ridge.\n")
			}
		}
		beta	= beta_init
	
		if (verbose==2){
			w_old = beta
			if (alpha==1) {
				// lasso verbose output
				printf("{txt}%8s %8s %10s %14s %14s\n","iter","shoots","n(w)","n(step)","f(w)")
			}
			else {
				// elastic net verbose output
				printf("{txt}%8s %8s %10s %14s %14s %14s %14s\n","iter","shoots","n1(w)","n1(step)","n2(w)","n2(step)","f(w)")
			}
		}
		
		m = 0
		while (m < maxIter)
		{
			beta_old = beta
			for (j = 1;j<=d.p;j++)
			{
				S0 = quadcolsum(XX[j,.]*beta) - XX[j,j]*beta[j] - Xy[j]
	
				if (alpha==1) {
					//  lasso
					if (d.lglmnet) {
						if (S0 > lambda*Psi[j]) {
							beta[j] = (lambda*Psi[j] - S0)/(XX[j,j])
						}
						else if (S0 < -lambda*Psi[j]) {
							beta[j] = (-lambda*Psi[j] - S0)/(XX[j,j])
						}
						else {
							beta[j] = 0
						}
					}
					else {
						if (S0 > lambda*Psi[j]) {
							beta[j] = (lambda*Psi[j] - S0)/(XX[j,j])
						}
						else if (S0 < -lambda*Psi[j]) {
							beta[j] = (-lambda*Psi[j] - S0)/(XX[j,j]) 
						}
						else {
							beta[j] = 0
						}
					}
				}								//  end lasso
				else if (alpha>0) {
					//  elastic net
					if (d.lglmnet) {
						if (S0 > lambda*Psi[j]*alpha) {
							beta[j] = (lambda*Psi[j]*alpha - S0)/(XX[j,j] + lambda2*Psi2[j]*(1-alpha))
						}
						else if (S0 < -lambda*Psi[j]*alpha) {
							beta[j] = (-lambda*Psi[j]*alpha - S0)/(XX[j,j] + lambda2*Psi2[j]*(1-alpha))
						}
						else {
							beta[j] = 0
						}
					}
					else {
						if (S0 > lambda*Psi[j]*alpha) {
							beta[j] = (lambda*Psi[j]*alpha - S0)/(XX[j,j] + lambda2*d.prestdy*Psi2[j]^2*(1-alpha))
						}
						else if (S0 < -lambda*Psi[j]*alpha) {
							beta[j] = (-lambda*Psi[j]*alpha - S0)/(XX[j,j] + lambda2*d.prestdy*Psi2[j]^2*(1-alpha))
						}
						else {
							beta[j] = 0
						}
					}
				}								//  end elastic net
				else if ((alpha==0)) {
					//  ridge - shooting not required for ridge since closed-formed solution exists
					errprintf("internal lassoshooting error\n")
					exit(499)
				}								//  end ridge
			}									//  end j loop over components of beta

			m++

			if (verbose==2) {
				// multiply by d.cons to zero out the mean if no constant
				mse			= mean( (((*d.y):-d.ymvec*d.cons) - ((*d.X):-d.mvec*d.cons)*beta):^2 )
				if (alpha==1) {
					// lasso
					if (d.lglmnet) {
						fobj	= 1/2*mse + lambda*Psi*abs(beta)*d.ysd
					}
					else {
						fobj	= mse + lambda/d.n*Psi*abs(beta)
					}
					printf("{res}%8.0g %8.0g %14.8e %14.8e %14.8e\n",					///
						m,																///
						m*d.p,															///
						colsum(abs(beta)),												///
						colsum(abs(beta-w_old)),										///
						fobj															///
						)
					w_old = beta
				}
				else {
					// elastic net
					if (d.lglmnet) {
						fobj	= 1/2*mse												///
									+     lambda * alpha   *Psi*abs(beta)*d.ysd			///
									+ 1/2*lambda2*(1-alpha)*(Psi2:^2)*(beta:^2)
					}
					else {
						fobj	= mse														///
									+     lambda/d.n *alpha    *Psi*abs(beta)				///
									+ 1/2*lambda2*d.prestdy/d.n*(1-alpha)*(Psi2:^2)*(beta:^2)
					}
					printf("{res}%8.0g %8.0g %14.8e %14.8e %14.8e %14.8e %14.8e\n",		///
						m,																///
						m*d.p,															///
						colsum(abs(beta)),												///
						colsum(abs(beta-w_old)),										///
						colsum(beta:^2),												///
						colsum(beta:^2-w_old:^2),										///
						fobj															///
						)
					w_old = beta
				
				}
			}

			// convergence should always be in standardized metric	    
			if (d.prestdflag) {
				if (quadcolsum(abs(beta-beta_old))<optTol)								break
			}
			else {
				if (quadcolsum( (abs(beta-beta_old)) :* d.sdvec' / d.ysd )<optTol)		break
			}
		}
		
		if ((verbose>=1) & (alpha~=0)) 	{
			printf("{txt}Number of iterations: {res}%g\n{txt}Total Shoots: {res}%g\n",m,m*d.p)
			if (m == maxIter) {
				printf("{txt}Warning: reached max shooting iterations w/o achieving convergence.\n")
			}
			else {
				printf("{txt}Convergence achieved.\n")
			}
		}
		
		if ((verbose==2) & (alpha~=0)) {
			printf("{txt}Initial beta and beta after estimation:\n")
			(beta_init'\beta')
		}	
	}
	else {	// square-root lasso

		// initialize beta with ridge if beta_init not provided
		if (noinitflag) {
			beta_init	= beta_ridge
			if (verbose>=1) {
				printf("{txt}Initial beta vector is ridge.\n")
			}
		}
		beta	= beta_init

		// d.cons=1 if there is a constant, =0 otherwise.
		// Demeaning below throughout is needed only if a constant is present,
		// hence we multiply means by d.cons.
		ERROR = ((*d.y):-(d.ymvec*d.cons)) - ((*d.X):-(d.mvec*d.cons))*beta
		Qhat = mean(ERROR:^2)
	
		if (verbose==2) {
			// lasso verbose output
			printf("{txt}%8s %8s %10s %14s %14s\n","iter","shoots","n(w)","n(step)","f(w)")
		}

		m = 0
		while (m < maxIter) {
			beta_old = beta
			for (j = 1;j<=d.p;j++) {
				S0 = quadcolsum(XX[j,.]*beta) - XX[j,j]*beta[j] - Xy[j]
	
				if ( abs(beta[j])>0 ) {
					ERROR = ERROR + ((*d.X)[.,j]:-(d.mvec*d.cons)[j])*beta[j]
					Qhat = mean(ERROR:^2)
				}
				
				if ( d.n^2 < (lambda * Psi[j])^2 / XX[j,j]) {
					beta[j] = 0
				}
				else if (S0 > lambda/d.n*Psi[j]*sqrt(Qhat))
				{
					beta[j]= ( ( lambda * Psi[j] / sqrt(d.n^2 - (lambda * Psi[j])^2 / XX[j,j] ) ) * sqrt(max((Qhat-(S0^2/XX[j,j]),0)))-S0)  / XX[j,j]
					ERROR = ERROR - ((*d.X)[.,j]:-(d.mvec*d.cons)[j])*beta[j]
				}
				else if (S0 < -lambda/d.n*Psi[j]*sqrt(Qhat))	
	       		{
					beta[j]= ( - ( lambda * Psi[j] / sqrt(d.n^2 - (lambda * Psi[j])^2 / XX[j,j] ) ) * sqrt(max((Qhat-(S0^2/XX[j,j]),0)))-S0)  / XX[j,j]
					ERROR = ERROR - ((*d.X)[.,j]:-(d.mvec*d.cons)[j])*beta[j]
	           	}
	       		else 
				{
	       			beta[j] = 0
				}
			}
	
			ERRnorm=sqrt( quadcolsum( (((*d.y):-(d.ymvec*d.cons))-((*d.X):-(d.mvec*d.cons))*beta):^2 ) )
			// sqrt-lasso objective function
			fobj = ERRnorm/sqrt(d.n) + (lambda/d.n)*Psi*abs(beta)
			
			if (ERRnorm>1e-10) {
				aaa = (sqrt(d.n)*ERROR/ERRnorm)
				dual = aaa'((*d.y):-(d.ymvec*d.cons))/d.n  - abs(lambda/d.n*Psi' - abs(((*d.X):-(d.mvec*d.cons))'aaa/d.n))'abs(beta)
			}
			else {
				dual = (lambda/d.n)*Psi*abs(beta)
			}
			
	       	m++
	
			if (verbose==2) {
					printf("{res}%8.0g %8.0g %14.8e %14.8e %14.8e\n",m,m*d.p,colsum(abs(beta)),colsum(abs(beta)),fobj)
					w_old = beta
			}
	
			// this block should be in standardized units
			if (d.prestdflag) {
				if (quadcolsum(abs(beta-beta_old))<optTol) {
					if (fobj - dual < 1e-6) {
						break
					}
				}
			}
			else {
				if (quadcolsum( (abs(beta-beta_old)) :* d.sdvec' / d.ysd )<optTol) {
					if ( (fobj - dual)/d.ysd < 1e-6) {
						break
					}
				}
			}
		}
		
		if (verbose>=1) {
			printf("{txt}Number of iterations: {res}%g\n{txt}Total Shoots: {res}%g\n",m,m*d.p)
			if (m == maxIter) {
				printf("{txt}Warning: reached max shooting iterations w/o achieving convergence.\n")
			}
			else {
				printf("{txt}Convergence achieved.\n")
			}
		}

	}
	
	// glmnet algo is in standardized metric, so need to unstandardize beta
	if ((d.lglmnet) & (d.prestdflag==0)) {
		beta	= beta :/ Psi' * d.ysd
	}

	// deviation (R-sq)
	ERROR	= ((*d.y):-(d.ymvec*d.cons)) - ((*d.X):-(d.mvec*d.cons))*beta

	RSS		= quadcolsum(ERROR:^2)
	RSQ		= 1-RSS:/d.TSS

	b.beta		= beta
	b.m			= m
	b.fobj		= fobj
	b.dev		= RSQ

	return(b)
}

	
void ReturnResultsPath(		struct outputStructPath scalar t,	/// #1
							struct dataStruct scalar d,			/// #2
							string scalar Xnames				/// #3
							)
{

		Xnamesall		= tokens(Xnames)
		betas			= t.betas
		// original (untruncated) and truncated list
		lambdalist0		= t.lambdalist0
		lambdalist		= t.lambdalist

		// standardize or unstandardize
		if (d.lglmnet) {
			if (d.prestdflag) {
				sbetas			= betas
				slambdalist0	= lambdalist0
				slambdalist		= lambdalist
				betas			= betas			:/ d.prestdx * d.prestdy
				lambdalist0		= lambdalist0	* d.prestdy
				lambdalist		= lambdalist	* d.prestdy
			}
			else {
				slambdalist0	= lambdalist0
				slambdalist		= lambdalist
				lambdalist0		= lambdalist0	* d.ysd
				lambdalist		= lambdalist	* d.ysd
			}
		}
		else if (d.prestdflag) {
			sbetas				= betas
			slambdalist0		= lambdalist0
			slambdalist			= lambdalist
			betas				= betas			:/ d.prestdx * d.prestdy
			if (d.sqrtflag==0) {
				// sqrt-lasso lambdas don't need unstandardizing (pivotal so doesn't depend on sigma/y)
				lambdalist0		= lambdalist0	* d.prestdy
				lambdalist		= lambdalist	* d.prestdy
			}
		}
		else {
			sbetas				= betas			:* d.sdvec / d.ysd
			slambdalist0		= lambdalist0	/ d.ysd
			slambdalist			= lambdalist	/ d.ysd
		}

		// "L1 norm" excluding constant and unpenalized vars
		// weighted L1 norm weights by penalty loadings (if prestandardized, by sd(X))
		l1norm				= rowsum(abs(betas :* (t.Psi :> 0)))
		if (d.prestdflag) {
			sl1norm				= l1norm / d.prestdy
		}
		else {
			sl1norm				= l1norm / d.ysd
		}
		if (d.prestdflag) {
			wl1norm			= rowsum(abs(betas) :* d.prestdx)
			swl1norm		= wl1norm / d.prestdy
		}
		else {
			wl1norm			= rowsum(abs(betas :* t.Psi))
			swl1norm		= wl1norm / d.ysd
		}

		pall				= cols(Xnamesall)

		// no intercept if pre-standardized
		// note that betas are in original metric
		// standardized betas do not have a constant so retain name list for these
		sXnamesall		= Xnamesall
		spall			= pall
		if (t.cons) {
			// note transposes
			intercept	= mean(*d.y) :- mean(*d.X)*betas'
			betas		= (betas , intercept')
			Xnamesall	= (Xnamesall, "_cons")
			pall		= pall+1
		}

		st_numscalar("r(lmax)", max(lambdalist))
		st_numscalar("r(lmax0)", max(lambdalist0))
		st_numscalar("r(lmin)", min(lambdalist))
		st_numscalar("r(lmin0)", min(lambdalist0))
		st_numscalar("r(lcount)",rows(lambdalist))
		st_matrix("r(lambdalist)",lambdalist)
		st_matrix("r(lambdalist0)",lambdalist0)
		st_matrix("r(slambdalist)",slambdalist)
		st_matrix("r(slambdalist0)",slambdalist0)
		st_matrix("r(l1norm)",l1norm)
		st_matrix("r(sl1norm)",sl1norm)
		st_matrix("r(wl1norm)",wl1norm)
		st_matrix("r(swl1norm)",swl1norm)
		st_matrix("r(betas)",betas)
		st_matrix("r(sbetas)",sbetas)
		st_matrix("r(Psi)",t.Psi)
		st_matrix("r(sPsi)",t.sPsi)
		st_matrix("r(shat)",t.shat)
		st_matrix("r(shat0)",t.shat0)
		st_matrix("r(stdvec)",d.sdvec)
		st_matrixcolstripe("r(betas)",(J(pall,1,""),Xnamesall'))
		st_matrixcolstripe("r(sbetas)",(J(spall,1,""),sXnamesall'))
		st_matrixcolstripe("r(lambdalist0)",("","Lambdas"))
		st_matrixcolstripe("r(lambdalist)",("","Lambdas"))
		st_matrixcolstripe("r(slambdalist0)",("","Lambdas"))
		st_matrixcolstripe("r(slambdalist)",("","Lambdas"))
		st_matrixcolstripe("r(l1norm)",("","L1norm"))
		st_matrixcolstripe("r(sl1norm)",("","sL1norm"))
		st_matrixcolstripe("r(wl1norm)",("","wL1norm"))
		st_matrixcolstripe("r(swl1norm)",("","swL1norm"))
}
// end ReturnResultsPath
	
void ReturnCVResults(real rowvector Lam, ///
					real rowvector MSPE,
					real rowvector SD, 
					real scalar minid,
					real scalar minseid)
{
	
	lnum=cols(Lam)
	
	printf("{txt}%10s{c |} {space 3} {txt}%10s {space 3} {txt}%10s {space 3} {txt}%10s\n","","Lambda","MSPE","st. dev.")
	printf("{hline 10}{c +}{hline 45}\n")
			
	for (j = 1;j<=lnum;j++) 	{
		
		if (j==minid) {
			marker="*"
		}
		else {
			marker=""
		}
		if (j==minseid) {
			marker=marker+"^"
		}
		printf("{txt}%10.0g{c |} {space 3} {res}%10.0g {space 3} {res}%10.0g {space 3} {res}%10.0g  %s\n",j,Lam[1,j],MSPE[1,j],SD[1,j],marker)
	
	}
}
// end ReturnCVResults


real scalar lambdaCalc(struct dataStruct scalar d,		///
						real scalar pminus,				/// #1  adjustment to number of Xs in model
						real scalar gamma,				/// #2  default is 0.1/log(N)
						real scalar c,					/// #3
						real scalar R,					/// #4
						real scalar xdep,				/// #5
						real colvector v,				/// #6
						real scalar rmse,				/// #7
						real rowvector ScoreStdVec,		/// #8
						real scalar lalt,				/// #9
						|								///
						real scalar newseed,			/// #10
						real scalar dotsflag,			/// #11
						real scalar verbose				/// #12
						) 
{

	if (args()<10)	newseed = -1
	if (args()<11)	dotsflag = 0
	if (args()<12)	verbose = 0

	// lasso gets a factor of 2c; sqrt-lasso gets a factor of c
	if (d.sqrtflag)	lassofactor = c
	else			lassofactor = 2*c

	// model may have partialled-out var with coeffs estimated separately; subtract in formulae below
	p = d.p - pminus

	if (xdep==0) {									// X-independent
		if (lalt==0) {								// standard lambda
			lambda	= lassofactor * sqrt(d.n)*invnormal(1-gamma/(2*p))
		}
		else {										//  alternative lambda
			lambda	= lassofactor * sqrt(d.n) * sqrt(2 * log(2*p/gamma))
		}
	}												//  end X-independent block
	else {											//  X-dependent
		multbs	= 0 // don't override
		sim		= SimMaxScoreDist(d,v,rmse,ScoreStdVec,multbs,verbose,R,newseed,dotsflag)
		// use quantile from simulated values to get lambda
		lambda	= lassofactor * d.n * mm_quantile(sim,1,1-gamma)
	}

	return(lambda)
}
// end lambdaCalc


real colvector InitialResiduals(	struct dataStruct scalar d,
									real scalar corrnumber) 
{ 
	// applies OLS to a reduced set of regressors exhibiting highest correlation with y
	// size of restricted set = corrnumber

	// in case corrnumber < dim(X)
	if (d.pihat==0) {
		dimX=cols(*d.X)
	}
	else {
		dimX=cols(*d.Xp)
	}
	corrnumber=min((corrnumber,dimX))

	// just return if corrnum = 0
	if (corrnumber <= 0) {
		// centerpartial(.) returns y centered and with Xnp partialled out.
		return(centerpartial(d.y,d.ymvec,d.cons,d.Xnp,d.mvecnp,d.ypihat))
	}
	else {
		dimZ=dimX+1
		// instead of official correlation(.), use m_quadcorr(.) defined below
		// m_quadcorr(.) accommodates case of zero-mean data or no constant
		Z=abs(m_quadcorr(																	///
							(centerpartial(d.y,d.ymvec,d.cons,d.Xnp,d.mvecnp,d.ypihat),		///
							 centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat))		///
							, d.cons	))
		z=Z[2..dimZ,1]
		ix=order(-z,1)
		ix=ix[1..corrnumber,1]

		if ((d.pihat==0) & (d.cons==0)) {				//  no notpen Xs, no constant in model or zero-mean data
			b	=qrsolve(								///
						(*d.X)[.,ix],					///
						 *d.y							///
						 )
			r	=(*d.y) :- ((*d.X)[.,ix])*b
		}
		else if (d.pihat==0) {							//  no notpen Xs, model has constant so demean in X'X
			b	=qrsolve(								///
						((*d.X)[.,ix]):-d.mvec[1,ix],	///
						(*d.y):-d.ymvec					///
						)
			r	=(*d.y) :- d.ymvec :- (((*d.X)[.,ix]):-d.mvec[1,ix])*b
		}
		else if (d.cons==0) {							//  notpen Xs, no constant in model or zero-mean data
			// ix has highest-correlation cols of Xp
			// replace with corresponding cols of X
			ix	= d.selindXp[1,ix]
			// and append np cols of X
			ix	= ix,d.selindXnp
			b	=qrsolve(								///
						 (*d.X)[.,ix],					/// selected cols of Xp plus all of Xnp
						 *d.y							///
						 )
			r	=(*d.y) :- ((*d.X)[.,ix])*b
		}
		else {											//  notpen Xs, model has constant so demean in X'X
			// ix has highest-correlation cols of Xp
			// replace with corresponding cols of X
			ix	= d.selindXp[1,ix]
			// and append np cols of X
			ix	= ix,d.selindXnp
			b	=qrsolve(								///
						 ((*d.X)[.,ix]):-d.mvec[1,ix],	/// selected cols of Xp plus all of Xnp
						 *d.y :- d.ymvec				///
						 )
			r	=(*d.y) :- d.ymvec :- (((*d.X)[.,ix]):-d.mvec[1,ix])*b
		}
		return(r)				// return residuals

	}
}
// end InitialResiduals


void EstimateRLasso(							///  Complete Mata code for RLasso.
				string scalar nameY,			///
				string scalar nameX,			///
				string scalar nameX_o,			///
				string scalar pen,				///
				string scalar notpen,			///
				string scalar touse,			///
				string scalar nameclustid, 		///
				string scalar stdymat,			///
				string scalar stdxmat,			///
				real scalar sqrtflag,			/// lasso or sqrt-lasso?
				real scalar hetero,				/// homosk or heteroskedasticity?
				real scalar bw,					/// bw for HAC; bw=0 means no HAC
				string scalar kernel,			///
				real scalar maqflag,			///
				real scalar spectral,			///
				string scalar tindexname,		/// time index for lags etc.
				real scalar tdelta,				///
				real scalar bsize,				///
				real scalar xdep,				/// X-dependent or independent?
				real scalar R,					/// number of simulations with xdep
				real scalar lassoPsiflag,		/// use lasso or post-lasso residuals for estimating penalty loadings?
				real scalar optTol,				///
				real scalar maxIter,			///
				real scalar zeroTol,			///
				real scalar maxPsiIter,			///
				real scalar PsiTol,				///
				real scalar verbose,			///
				real scalar c,					///
				real scalar c0,					///
				real scalar gamma,				///
				real scalar gammad,				///
				real scalar lambda0,			///
				real scalar lalt, 				///
				real scalar corrnumber,			///
				real scalar pminus,				///
				real scalar nclust1flag,		/// use #nclust-1 instead of #nclust in cluster-lasso
				real scalar center,				/// center x_i*e_i or x_ij*e_ij
				real scalar sigma,				/// user-specified sigma for standard lasso
				real scalar supscoreflag,		///
				real scalar ssnumsim,			///
				real scalar ssgamma,			///
				real scalar newseed,			/// rnd # seed; relevant for xdep and supscore
				real scalar dotsflag,			///
				real scalar cons,				///
				real scalar dmflag,				///
				real scalar dofminus,			///
				real scalar sdofminus,			///
				real scalar prestdflag,			///
				string scalar psolver)
{

// temporary measure
	lglmnet = 0

	struct dataStruct scalar d
	d = MakeData(nameY,nameX,nameX_o,touse,cons,dmflag,dofminus,sdofminus,prestdflag,sqrtflag,stdymat,stdxmat,lglmnet,nameclustid,hetero,center,nclust1flag,pen,notpen,bw,kernel,maqflag,spectral,tindexname,tdelta,bsize,psolver)

	struct outputStruct scalar OUT
	if (d.sqrtflag) {
		OUT	= RSqrtLasso(d,xdep,R,lassoPsiflag,optTol,maxIter,zeroTol,maxPsiIter,PsiTol,verbose,c,c0,gamma,gammad,lambda0,lalt,corrnumber,pminus,supscoreflag,ssnumsim,ssgamma,newseed,dotsflag)
	}
	else {
		OUT	= RLasso(d,xdep,R,lassoPsiflag,optTol,maxIter,zeroTol,maxPsiIter,PsiTol,verbose,c,c0,gamma,gammad,lambda0,lalt,corrnumber,pminus,sigma,supscoreflag,ssnumsim,ssgamma,newseed,dotsflag)
	}
	ReturnResults(OUT,d)		//  puts results into r(.) macros

}	// end EstimateRLasso


struct outputStruct scalar RLasso(							/// Mata code for BCH rlasso
							struct dataStruct scalar d,		/// data
							real scalar xdep,				/// X-dependent or independent?
							real scalar R,					/// number of simulations with xdep
							real scalar lassoPsiflag,		/// use lasso or post-lasso residuals for estimating penalty loadings?
							real scalar optTol,				///
							real scalar maxIter,			///
							real scalar zeroTol,			///
							real scalar maxPsiIter,			///
							real scalar PsiTol,				///
							real scalar verbose,			///
							real scalar c,					///
							real scalar c0,					///
							real scalar gamma,				///
							real scalar gammad,				///
							real scalar lambda0,			///
							real scalar lalt, 				///
							real scalar corrnumber,			///
							real scalar pminus,				///
							real scalar sigma,				/// user-specified sigma (default=0)
							real scalar supscoreflag,		///
							real scalar ssnumsim,			///
							real scalar ssgamma,			///
							real scalar newseed,			///
							real scalar dotsflag			///
							)
{

	struct outputStruct scalar betas

	alpha=1 // elastic net parameter. always one.
	if (gammad<=0) {													//  not user-provided, so set here
		if (d.nclust)						gammad=log(d.nclust)		//  cluster-lasso so #obs=nclust
		else if ((d.bw>0) & (d.maqflag))	gammad=log(d.n/(d.bw+1))	//  HAC case, truncated kernel
		else								gammad=log(d.n)				//  not cluster or truncated kerne so #obs=n 
	}
	if (gamma<=0) {														//  not user-provided, so set here
		gamma	= 0.1/gammad
	}
	
	if (sigma>0) {		// user-supplied sigma; assumes homoskedasticity

		if (verbose>=1) {
			printf("{txt}Estimation of penalty level/loadings with user-supplied sigma and i.i.d. data\n")
			printf("{txt}Obtaining lasso estimate...\n")
		}
		Psi		= d.sdvecpnp*sigma
		// note we use c and not c0
		lambda	= lambdaCalc(d,pminus,gamma,c,R,xdep,v,sigma,Psi,lalt,newseed,dotsflag,verbose)
		// initial estimator is ridge
		beta_ridge	= anysolver(d.XX*2+lambda*diag(d.sdvec:^2),d.Xy*2, r=.)
		if (verbose>=2) {
			printf("{txt}Initial estimator is ridge:\n")
			beta_ridge'
		}
		betas	= DoLasso(d, Psi, Psi, lambda, lambda, verbose, optTol, maxIter, zeroTol, alpha, beta_ridge)
		if (verbose>=1) {
			printf("{txt}Selected variables: {res}%s\n\n",invtokens(betas.nameXSel'))
		}

	}
	else {				// standard code block to estimate rmse/lambda/penalty loadings
	
		// initial residuals
		v		= InitialResiduals(d,corrnumber)
		s1		= sqrt(mean(v:^2))

		// Create the score standardization vector using initial residuals.
		// Dim = 1 x p with zeros for unpenalized Xs.
		Psi		= MakeScoreStdVec(d,v,s1)
		
		// initialize lambda
		if (lambda0) {
			// user-supplied lambda0
			lambda=lambda0
		}
		else {
			// lambda does not incorporate rmse, does incorporate lasso factor of 2
			// note we use c0 for the first lambda
			lambda	= lambdaCalc(d,pminus,gamma,c0,R,xdep,v,s1,Psi,lalt,newseed,dotsflag,verbose)
		}
		// "iteration 1" - get first lasso estimate based on initial lambda/loadings
		iter = 1
		if (verbose>=1) {
			printf("{txt}Estimation of penalty level/loadings: Step {res}%f.\n",iter)
			printf("{txt}Obtaining initial lasso estimate...\n")
		}
		// initial estimator is ridge
		beta_ridge	= anysolver(d.XX*2+lambda*diag(d.sdvec:^2),d.Xy*2, r=.)
		if (verbose>=2) {
			printf("{txt}Initial estimator is ridge:\n")
			beta_ridge'
		}

		betas	= DoLasso(d, Psi, Psi, lambda, lambda, verbose, optTol, maxIter, zeroTol, alpha, beta_ridge)
		if (verbose>=1) {
			printf("{txt}Selected variables: {res}%s\n\n",invtokens(betas.nameXSel'))
		}

		// initialize Delta
		Delta = 1e10
		
		while ((iter < maxPsiIter) & (Delta > PsiTol)) {

			s0	= s1
			// obtain residuals; based on betas(Psi)
			if (lassoPsiflag) {
				v	= betas.v		// lasso residuals
				s1	= betas.rmse 
			}
			else {
				v	= betas.vPL 	// post-lasso residuals
				s1	= betas.rmsePL
			}
			// change in RMSE used in loadings
			Delta = abs(s1-s0)
			// new loadings and lambda
			Psi	= MakeScoreStdVec(d,v,s1)
			// note we use c and not c0 from now on
			lambda=lambdaCalc(d,pminus,gamma,c,R,xdep,v,s1,Psi,lalt,newseed,dotsflag,verbose)

			// Reporting
			if (verbose>=1) {
				printf("{txt}Estimation of penalty level/loadings: Step {res}%f.\n",iter)
				printf("{txt}RMSE: {res}%f\n",s1)
				printf("{txt}Change in RMSE: {res}%f\n",Delta)
				printf("{txt}Obtaining new lasso estimate...\n")
			}
			// new lasso estimate
			betas = DoLasso(d, Psi, Psi, lambda, lambda, verbose, optTol, maxIter, zeroTol, alpha, beta_ridge)
			if (verbose>=1) {
				printf("{txt}Selected variables: {res}%s\n\n",invtokens(betas.nameXSel'))
			}
		
			iter++
		
		}
	
	}	// end of code to iterate penalty loadings/lambda

	if (verbose>=1) {
		printf("{txt}Number of penalty loading iterations: {res}%g\n",iter)
		if (iter == maxPsiIter) {
			printf("{txt}Warning: reached max penalty loading iterations w/o achieving convergence.\n")
		}
		else {
			printf("{txt}Penalty loadings (psi) convergence achieved.\n")
		}
	}

	// sup-score stat
	if (supscoreflag) {
		betas.supscore	= doSupScore(d, c, ssgamma, pminus, verbose, ssnumsim, newseed, dotsflag)
		betas.ssgamma	= ssgamma
	}

	// convention is lambda0 = penalty without rmse (similar to sqrt-lambda);
	//               lambda = penalty including rmse;
	//               Psi = penalty loadings excluding rmse but standardizing;
	//               sPsi = penalty loadings excluding rmse and standardization; =1 under homoskedasticity.
	betas.lambda0	= lambda
	betas.lambda	= lambda * s1
	betas.slambda	= lambda * s1 / d.ysdp
	betas.Psi		= betas.Psi / s1
	betas.sPsi		= betas.Psi :/ d.sdvecpnp				//  should be =1 under homosk.
	betas.sPsi		= editmissing(betas.sPsi,0)				//  in case any unpenalized (div by 0)

	// Misc
	betas.npsiiter	= iter
	betas.c			= c
	betas.gamma		= gamma
	betas.gammad	= gammad

	return(betas)
}
// end RLasso	

void EstimateSupScore(							/// 
				string scalar nameY,			///
				string scalar nameX,			///
				string scalar nameX_o,			///
				string scalar pen,				///
				string scalar notpen,			/// #5
				string scalar touse,			///
				string scalar nameclustid, 		///
				string scalar stdymat,			///
				string scalar stdxmat,			///
				real scalar sqrtflag,			/// #10; will always be 0
				real scalar hetero,				/// homosk or heteroskedasticity?
				real scalar bw,					/// bw for HAC; bw=0 means no HAC
				string scalar kernel,			///
				real scalar maqflag,			///
				real scalar spectral,			/// #15
				string scalar tindexname,		/// time index for lags etc.
				real scalar tdelta,				///
				real scalar bsize,				///
				real scalar verbose,			///
				real scalar R,					/// #20
				real scalar ssiidflag,			/// override use of multiplier bootstrap in iid case
				real scalar c,					///
				real scalar nclust1flag,		///
				real scalar center,				///
				real scalar ssgamma,			/// #25
				real scalar pminus,				///
				real scalar newseed,			///
				real scalar dotsflag,			///
				real scalar cons,				///
				real scalar dmflag,				/// #30
				real scalar dofminus,			///
				real scalar sdofminus,			///
				real scalar prestdflag,			///
				string scalar psolver)
{

// temporary measure
	lglmnet = 0

	struct dataStruct scalar d
	d = MakeData(nameY,nameX,nameX_o,touse,cons,dmflag,dofminus,sdofminus,prestdflag,sqrtflag,stdymat,stdxmat,lglmnet,nameclustid,hetero,center,nclust1flag,pen,notpen,bw,kernel,maqflag,spectral,tindexname,tdelta,bsize,psolver)

	struct outputStruct scalar OUT
	OUT.supscore	= doSupScore(d, c, ssgamma, pminus, verbose, R, newseed, dotsflag, ssiidflag)
	// Misc
	OUT.c			= c
	OUT.ssgamma		= ssgamma

	ReturnResults(OUT,d)		//  puts results into r(.) macros

}	// end EstimateRLasso

real colvector SimMaxScoreDist(									///					
								struct dataStruct scalar d,		/// #1
								real colvector v,				/// #2
								real scalar rmse,				/// #3
								real rowvector ScoreStdVec,		/// #4
								|								///
								real scalar multbs,				/// #5
								real scalar verbose,			/// #6
								real scalar R,					/// #7
								real scalar newseed,			/// #8
								real scalar dotsflag			/// #9
								)
{

	// defaults
	// default is multbs=0 => use multiplier bootstrap for i.n.i.d. or cluster but not for i.i.d.
	if (args()<5)	multbs		= 0
	if (args()<6)	verbose		= 0
	if (args()<7)	R			= 500
	if (args()<8)	newseed		= -1
	if (args()<9)	dotsflag	= 0
	//  set seed if requested
	if (newseed>-1)	rseed(newseed)

	if (dotsflag) {
		dotscmd = "_dots 0 0, title(Estimating score vector distribution using " + strofreal(R) + " repetitions)"
		stata(dotscmd)
	}
	
	// initialize
	sim=J(R,1,.)
	// initialize if needed for HAC/AC
	// gaps flag is set to 1 if we are using HAC/AC and there are gaps in the time series
	if (d.bw>0) {
		tnow		= st_data(., d.tindexname)
		gapsflag	= (rows(tnow) > rows(v))
	}
	else {
		gapsflag	= 0
	}
	// initialize n_mbs; will be updated to account for gaps in the HAC/AC loop below
	n_mbs			= d.n
	
	for (j=1; j<=R; j++) {
		if (dotsflag) {
			dotscmd = "_dots " + strofreal(j) + " 0"
			stata(dotscmd)
		}
		if ((d.nclust==0) & (d.bw==0)) {
			// standard non-clustered case - g is iid standard normal
			g=rnormal(d.n,1,0,1)
		}
		else if (d.nclust>0) {
			// g is iid by cluster and repeated within clusters
			info	= panelsetup(*d.clustid1, 1)
			// g is iid by cluster and repeated within clusters
			g=J(d.n,1,.)
			info = panelsetup(*d.clustid1, 1)
			for (i=1; i<=d.nclust1; i++) {
				g[info[i,1]..info[i,2],1] = J(info[i,2]-info[i,1]+1,1,rnormal(1,1,0,1))
			}
			// 2nd dimension of 2-way clustering
			if (d.nclust2>0) {
				g2=J(d.n,1,0)
				for (i=1; i<=d.nclust2; i++) {
					svar		= selectindex(*d.clustid2:==i)
					g2[svar,1]	= J(rows(svar),1,rnormal(1,1,0,1))
				}
				g	= g :* g2
			}
		}
		else {
			// HAC/AC case
			// g is iid by block and repeated within blocks
			g=J(d.n,1,.)
			// if the last block doesn't have a full set of observations, it's ignored
			numblocks	= floor(rows(tnow)/d.bsize)
			for (i=1; i<=numblocks; i++) {
				svar		= tnow[(1+(i-1)*d.bsize)..(i*d.bsize),1]
				svar		= select(svar,svar:<.)
				if (rows(svar)==d.bsize) {
					// bootstrap is based only on blocks with a full set of observations
					g[svar,1]	= J(rows(svar),1,rnormal(1,1,0,1))
				}
			}
		}
		
		if (gapsflag) {
			// since bootstrap is based only on full blocks, need an n which is #nonmissings
			// overwrites n_mbs that was initialized to d.n
			// note that we use n_mbs = d.n if there are no gaps, even though the last block
			// might not have a full set of observations (and hence g has missings)
			n_mbs	= colnonmissing(g)
		}

		// centerpartial(.) returns Xp demeaned and with Xnp partialled-out.
		if (d.sqrtflag==0) {
			// standard lasso or score
			if ((multbs) | (d.hetero) | (d.nclust)) {
				// use multiplier bootstrap
				ScoreVec		= 1/(n_mbs) * quadcolsum( centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat) :* (g:*v) )
			}
			else {
				// use i.i.d. version; rmse scales g so that it's the same scale as v
				ScoreVec		= 1/(n_mbs) * quadcolsum( centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat) :* (g*rmse) )
			}
		}
		else {
			// sqrt lasso
			if ((multbs) | (d.hetero) | (d.nclust)) {
				// use multiplier bootstrap
				ScoreVec		= 1/(n_mbs) * quadcolsum( centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat) :* (g:*v) * 1/rmse ) * 1/sqrt(mean(g:^2))
			}
			else {
				// pivotal so we don't need rmse, but we do need to normalize by sd of g (will be appx 1 anyway).
				ScoreVec		= 1/(n_mbs) * quadcolsum( centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat) :* g ) * 1/sqrt(mean(g:^2))
			}
		}

		// Now create blank full-size score vector for all X variables.
		FullScoreVec					= J(1,cols(*d.X),0)
		// Then insert values in appropriate columns.
		FullScoreVec[1,(d.selindXp)]	= ScoreVec
		sim[j,1]						= max(abs(FullScoreVec:/ScoreStdVec))
	}

	return(sim)
}

// Score standardization vector:
// 1. Standard lasso + homosked = SD(x)*rmse
// 2. Sqrt lasso + homosked = SD(x)
// Robust versions are in the same metric.
real rowvector MakeScoreStdVec(								///
								struct dataStruct scalar d,	///
								real colvector v,			///
								real scalar rmse			///
								)
{

	real rowvector ScoreStd									//  will be partialled-out version with zero penalties inserted

	if (d.nclust==0) {										//  no cluster dependence
		if ((d.hetero==0) & (d.bw==0)) {					//  homoskedastic
			if (d.sqrtflag) {
				// sqrt-lasso case
				ScoreStd = d.sdvecpnp
			}
			else {
				// standard case
				ScoreStd = d.sdvecpnp*rmse
			}
		}
		else if ((d.hetero) & (d.bw==0)) {
			// heteroskedastic case, independence
			St = J(1,cols(*d.X),0)
			if (d.center) {		
				// centerpartial(.) returns Xp demeaned and with Xnp partialled-out.
				centervec			= mean( centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat) :* v)
				St[1,(d.selindXp)]	= quadcolsum( ( (centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat) :* v) :- centervec ):^2 )
			}
			else {
				St[1,(d.selindXp)]	= quadcolsum( ( (centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat) :* v)              ):^2 )
			}
			if (d.sqrtflag) {
				// sqrt-lasso version
				ScoreStd = sqrt(St/d.n) * 1/rmse
			}
			else {
				ScoreStd = sqrt(St/d.n)
			}
		}
		else if ((d.hetero) & (d.bw>0)) {
			// HAC case
			St = J(1,cols(*d.X),0)

			// Gamma0; same as het-robust
			if (d.center) {
				centervec	= mean( centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat) :* v)
				Gamma0		= quadcolsum( ( (centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat) :* v) :- centervec ):^2 ) * 1/d.n
			}
			else {
				Gamma0		= quadcolsum( ( (centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat) :* v)              ):^2 ) * 1/d.n
			}

			St[1,(d.selindXp)] = Gamma0

			// initialize
			if (d.spectral)		TAU=d.n/d.tdelta-1
			else				TAU=d.bw

			for (s=1; s<=TAU; s++) {
				tau=s
				lstau = "L"+strofreal(tau)
				tnow=st_data(., d.tindexname)
				tlag=st_data(., lstau+"."+d.tindexname)
				tmatrix = tnow, tlag
				svar=(tnow:<.):*(tlag:<.)		// multiply column vectors of 1s and 0s
				tmatrix=select(tmatrix,svar)	// to get intersection, and replace tmatrix
				kw = m_calckw(tau, d.bw, d.kernel)

				if (rows(tmatrix)>0) {
					if (d.center) {
						Gamma = quadcolsum(																							///
							(((centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat)[tmatrix[.,1],.]) :* v[tmatrix[.,1],.]) :- centervec) :* 	///
							(((centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat)[tmatrix[.,2],.]) :* v[tmatrix[.,2],.]) :- centervec)		///
							) * 1/d.n
					}
					else {
						Gamma = quadcolsum(																							///
							(centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat)[tmatrix[.,1],.] :* v[tmatrix[.,1],.]) :* 	///
							(centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat)[tmatrix[.,2],.] :* v[tmatrix[.,2],.])		///
							) * 1/d.n
					}
					St[1,(d.selindXp)] = St[1,(d.selindXp)] + 2*kw*Gamma
				}
			}
			// handle case if some variances negative (possible in finite samples)

			if (sum(St :< 0)) {
				Stzeros			= (St :< 0)
				Stznames		= invtokens(select(d.nameX_o',Stzeros')')
				printf("Warning: negative penalty loading(s) encountered.\n")
				printf("Replacing gamma0+2*kw*gamma with gamma0+kw*gamma.\n")
				printf("Variables affected: %s\n\n",Stznames)
				d.psinegs		= d.psinegs + sum(St :< 0)
				d.psinegvars	= d.psinegvars + " " + Stznames
				// add gamma0 to (gamma0+2*kw*gamma) and then divide by 2
				St[1,selectindex(Stzeros)] = St[1,selectindex(Stzeros)] + (Gamma0[1,selectindex(Stzeros)])
				St[1,selectindex(Stzeros)] = St[1,selectindex(Stzeros)] * 0.5
			}

			if (d.sqrtflag) {
				// sqrt-lasso version
				ScoreStd = sqrt(St) * 1/rmse
			}
			else {
				ScoreStd = sqrt(St)
			}
		}
		else  {
			// AC case
			// Gamma0; same as iid case
			Gamma0 = (d.sdvecpnp*rmse) :^2
			St = Gamma0

			// initialize
			if (d.spectral)		TAU=d.n/d.tdelta-1
			else				TAU=d.bw

			for (s=1; s<=TAU; s++) {
				tau=s
				lstau = "L"+strofreal(tau)
				tnow=st_data(., d.tindexname)
				tlag=st_data(., lstau+"."+d.tindexname)
				tmatrix = tnow, tlag
				svar=(tnow:<.):*(tlag:<.)		// multiply column vectors of 1s and 0s
				tmatrix=select(tmatrix,svar)	// to get intersection, and replace tmatrix
				kw = m_calckw(tau, d.bw, d.kernel)

				if (rows(tmatrix)>0) {
					GammaX = quadcolsum(																///
						(centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat)[tmatrix[.,1],.]) :*	///
						(centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat)[tmatrix[.,2],.])		///
						) * 1/d.n
					GammaV = quadcross( v[tmatrix[.,1],.], v[tmatrix[.,2],.] )*1/d.n
					Gamma = GammaX * GammaV
					St[1,(d.selindXp)] = St[1,(d.selindXp)] + 2*kw*Gamma
				}
			}
			if (sum(St :< 0)) {
				Stzeros			= (St :< 0)
				Stznames		= invtokens(select(d.nameX_o',Stzeros')')
				printf("Warning: negative penalty loading(s) encountered.\n")
				printf("Replacing gamma0+2*kw*gamma with gamma0+kw*gamma.\n")
				printf("Variables affected: %s\n\n",Stznames)
				d.psinegs		= d.psinegs + sum(St :< 0)
				d.psinegvars	= d.psinegvars + " " + Stznames
				// add gamma0 to (gamma0+2*kw*gamma) and then divide by 2
				St[1,selectindex(Stzeros)] = St[1,selectindex(Stzeros)] + (Gamma0[1,selectindex(Stzeros)])
				St[1,selectindex(Stzeros)] = St[1,selectindex(Stzeros)] * 0.5
			}

			if (d.sqrtflag) {
				// sqrt-lasso version
				ScoreStd = sqrt(St) * 1/rmse
			}
			else {
				ScoreStd = sqrt(St)
			}
		}
	}
	else {													//  cluster dependence
		// relevant n is #clusters
		info = panelsetup(*d.clustid1, 1)
		// StTemp is n x p. Each row has the cluster-sum of x_it*e_it for cluster i.
		// First create blank n x p matrix of zeros to be populated.
		StTemp = J(d.nclust1,cols(*d.X),0)
		// Now populate correct columns, leaving other as zeros (unpenalized).
		// centerpartial(.) returns Xp demeaned and with Xnp partialled-out.
		StTemp[.,(d.selindXp)] = panelsum(centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat):*(v*J(1,cols(*d.Xp),1)),info)
		//  center StTemp
		if (d.center) {
			StTemp = StTemp - J(d.nclust1,1,1)*quadcolsum(StTemp)/d.nclust1
		}
		//  Default approach is simply to divide by nobs = nclust*T in a balanced panel.
		//  A finite-sample adjustment as in CBH's lassoCluster is to divide by (nclust-1)*T.
		//  In an unbalanced panel, we achieve this with 1/nobs * nclust/(nclust-1).
		//  In a balanced panel, 1/nobs * nclust/(nclust-1) = 1/(nclust*T) * nclust/(nclust-1) = 1/((nclust-1)*T).
		if (d.nclust1flag) {								//  override default: divide by (nclust-1)*T
			St = quadcolsum(StTemp:^2)/(d.n)*(d.nclust1)/(d.nclust1-1)
		}
		else {												//  default: simply divide by n = nobs*T
			St = quadcolsum(StTemp:^2)/(d.n)
		}

		if (d.nclust2>0) {

			// Dimension #2. Not sorted so need to use custom function.
			StTemp = J(d.nclust2,cols(*d.X),0)
			StTemp[.,(d.selindXp)] = uspanelsum(centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat):*(v*J(1,cols(*d.Xp),1)),*d.clustid2,d.nclust2)
			//  center StTemp
			if (d.center) {
				StTemp = StTemp - J(d.nclust2,1,1)*quadcolsum(StTemp)/d.nclust2
			}
			if (d.nclust1flag) {								//  override default: divide by (nclust-1)*T
				St2 = quadcolsum(StTemp:^2)/(d.n)*(d.nclust2)/(d.nclust2-1)
			}
			else {												//  default: simply divide by n = nobs*T
				St2 = quadcolsum(StTemp:^2)/(d.n)
			}

			// Dimension #3 is intersection of #1 and #2. panelsetup(.) works since sorted.
			info = panelsetup(*d.clustid3, 1)
			StTemp = J(d.nclust3,cols(*d.X),0)
			StTemp[.,(d.selindXp)] = panelsum(centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat):*(v*J(1,cols(*d.Xp),1)),info)
			//  center StTemp
			if (d.center) {
				StTemp = StTemp - J(d.nclust3,1,1)*quadcolsum(StTemp)/d.nclust3
			}
			if (d.nclust1flag) {								//  override default: divide by (nclust-1)*T
				St3 = quadcolsum(StTemp:^2)/(d.n)*(d.nclust3)/(d.nclust3-1)
			}
			else {												//  default: simply divide by n = nobs*T
				St3 = quadcolsum(StTemp:^2)/(d.n)
			}
			
			// 2-level clustering, completion (Cameron-Gelbach-Miller/Thompson)
			// Add 2 cluster variance matrices and subtract 3rd
			St	= St + St2 - St3

			// handle case if some variances negative (possible in finite samples)
			if (sum(St :< 0)) {
				Stzeros			= (St :< 0)
				Stznames		= invtokens(select(d.nameX_o',Stzeros')')
				printf("Warning: negative penalty loading(s) encountered.\n")
				printf("Replacing V1+V2-V3 with V1+V2-(1/2)V3.\n")
				printf("Variables affected: %s\n\n",Stznames)
				d.psinegs		= d.psinegs + sum(St :< 0)
				d.psinegvars	= d.psinegvars + " " + Stznames
				// add 1/2 V3
				St[1,selectindex(Stzeros)] = St[1,selectindex(Stzeros)] + 0.5*(St3[1,selectindex(Stzeros)])
			}

		}
		
		ScoreStd		= sqrt(St)
		if (d.sqrtflag) {
			// sqrt-lasso version
			ScoreStd	= ScoreStd * 1/rmse
		}

	}

	return(ScoreStd)
}
// end MakeScoreStdVec


real rowvector doSupScore(								///
							struct dataStruct scalar d,	/// #1
							real scalar c,				/// #2
							real scalar ssgamma,		/// #3
							real scalar pminus,			/// #4
							|							///
							real scalar verbose,		/// #5
							real scalar R,				/// #6
							real scalar newseed,		/// #7
							real scalar dotsflag,		/// #8
							real scalar ssiidflag		/// #9 undocumented option - override use of mult bootstrap
							)
{

	// defaults
	if (args()<5)	verbose = 0
	if (args()<6)	R = 500
	if (args()<7)	newseed = -1
	if (args()<8)	dotsflag = 0
	if (args()<9)	ssiidflag = 0

// ******* sup-score test stat ********** //
// Matlab code from BCH.
// Loop over jj values of aVec: assemble ynull (="eTemp") then calc SupScore for that jj ynull
// % Sup-Score Test
// aVec = (-.5:.001:1.5)';
// SupScore = zeros(size(aVec));
// for jj = 1:size(aVec,1)
//     aT = aVec(jj,1);
//     eTemp = My-Md*aT;
//     ScoreVec = eTemp'*Mz;
//     ScoreStd = sqrt((eTemp.^2)'*(Mz.^2));
//     ScaledScore = ScoreVec./(1.1*ScoreStd);
//     SupScore(jj,1) = max(abs(ScaledScore));
// end
//
// hdm rlasso R code:
//     object$supscore <- sqrt(n)*max(abs(colMeans(object$model*as.vector(object$dev))))

// ******* sup-score crit values ********** //
// Following rlasso code:
//   object$supscore <- sqrt(n)*max(abs(colMeans(object$model*as.vector(object$dev))))
//    R <- 500
//    stat <- vector("numeric", length=R)
//    for (i in 1:R) {
//      g <- rnorm(n)
//      dev.g <- as.vector(g*object$dev)
//      mat <- object$model*dev.g
//      stat[i] <- sqrt(n)*max(abs(colMeans(mat)))
//    }
//    object$pvalue <- sum(stat>object$supscore)/R

	// First generate v = epsilon under H0: all betas=0, after partialling out constant and unpenalized vars.
	// rmse is just mean v^2.
	if (d.cons) {
		// standard case - y is demeaned and standardized
		if (d.np) {
			v		= ((*d.y) :- d.ymvec) - (((*d.Xnp):-mean(*d.Xnp))*d.ypihat)
		}
		else {
			v		= (*d.y) :- d.ymvec
		}
		rmse	= d.ysdp
	}
	else if (d.dmflag) {
		// zero-mean data and no constant is present in the model
		// means are zero vectors, cross-prods don't need demeaning, SDs don't need demeaning
		if (d.np) {
			v		= (*d.y)-(*d.Xnp)*d.ypihat
		}
		else {
			v		= *d.y
		}
		rmse	= d.ysdp
	}
	else {
		// model has no constant but means may be nonzero
		if (d.np) {
			v		= centerpartial(d.y,d.ymvec,d.cons,d.Xnp,d.mvecnp,d.ypihat)						///
						:- mean(centerpartial(d.y,d.ymvec,d.cons,d.Xnp,d.mvecnp,d.ypihat))
		}
		else {
			v		= *d.y
		}
		rmse	= sqrt(mean(v:^2))
	}

	// *** supscore statistic *** //
	// NB: can't use quadcross or mean with code below because a column can have (all) missing (e.g. if std dev=0).
	//     and quadcross and mean use rowwise deletion (all missing => all rows dropped!)
	//     use 1/n * quadcolsum instead
	// Create the score standardization vector. Dim = 1 x p with zeros for unpenalized Xs.
	ScoreStdVec	= MakeScoreStdVec(d,v,rmse)
	// Create the unstandardized score vector. Unpenalized vars only.
	ScoreVec	= 1/(d.n) * quadcolsum( centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat) :* v )
	// Now create blank full-size score vector for all X variables.
	FullScoreVec					= J(1,cols(*d.X),0)
	// Then insert values in appropriate columns.
	FullScoreVec[1,(d.selindXp)]	= ScoreVec
	// Sup-score statistic
	supscore						= sqrt(d.n)*max(abs(FullScoreVec:/ScoreStdVec))

	// *** supscore pvalue/critical value *** //
	if (ssiidflag)	multbs	= 0		// override use of multiplier bootstrap for iid case
	else			multbs	= 1		// always use multiplier bootstrap for supscore unless overridden
	// simulate distribution of the maximal element of the score vector
	// returns a colvector with R obs from the simulated distribution
	sim = sqrt(d.n)*SimMaxScoreDist(d,v,rmse,ScoreStdVec,multbs,verbose,R,newseed,dotsflag)
	supscore_pvalue = sum(supscore:<sim)/R
	supscore_critvalue	=c*invnormal(1-ssgamma/(2*((d.p)-pminus)))

	res = (supscore, supscore_pvalue, supscore_critvalue, ssgamma)
	return(res)
}


struct outputStruct scalar RSqrtLasso(						/// Mata code for BCH sqrt rlasso
							struct dataStruct scalar d,		/// data
							real scalar xdep,				/// X-dependent or independent?
							real scalar R,					/// number of simulations with xdep
							real scalar lassoPsiflag,		/// use lasso or post-lasso residuals for estimating penalty loadings?
							real scalar optTol,				///
							real scalar maxIter,			///
							real scalar zeroTol,			///
							real scalar maxPsiIter,			///
							real scalar PsiTol,				///
							real scalar verbose,			///
							real scalar c,					///
							real scalar c0,					///
							real scalar gamma,				///
							real scalar gammad,				///
							real scalar lambda0,			///
							real scalar lalt, 				///
							real scalar corrnumber,			///
							real scalar pminus,				///
							real scalar supscoreflag,		///
							real scalar ssnumsim,			///
							real scalar ssgamma,			///
							real scalar newseed,			///
							real scalar dotsflag			///
							)
{

	struct outputStruct scalar betas

	alpha		= 1		// elastic net parameter. always one.

	if (gammad<=0) {										//  not user-provided, so set here
		if (d.nclust==0)	gammad=log(d.n)					//  not cluster-lasso so #obs=n
		else				gammad=log(d.nclust)			//  cluster-lasso so #obs=nclust
	}
	if (gamma<=0) {											//  not user-provided, so set here
		gamma	= 0.1/gammad
	}

	iter	= 1
	if (verbose>=1) {
		printf("Obtaining penalty level/loadings: Step %f.\n",iter)
	}

	if (d.hetero | d.nclust) {
		if (corrnumber>=0) {
			// corrnumber specified so use initial residuals for initial loadings
			v		= InitialResiduals(d,corrnumber)
			s1		= sqrt(mean(v:^2))
			Psi		= MakeScoreStdVec(d,v,s1)
			Psi		= (Psi :< d.sdvecpnp):*d.sdvecpnp + (Psi :>= d.sdvecpnp):*Psi
		}
		else {
			// corrnumber=-1 indicates use initial penalty loadings = colmax(abs(X)) as per BCW 2014 p. 769
			// Psi vector with zeros; cols corresponding to penalized vars updated below
			Psi		= J(1,cols(*d.X),0)
			if (d.nclust) {
				// initial loadings = colmax of panel means of abs(X)
				info					= panelsetup(*d.clustid, 1)
				psum					= panelsum(abs(centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat)),info)
				pmean					= psum :/ (info[.,2]-info[.,1]:+1)
				Psi[1,(d.selindXp)]		= colmax(pmean)
			}
			else {
				// initial loadings = colmax(abs(X))
				Psi[1,(d.selindXp)]		= colmax(abs(centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat)))
			}
			v		= .
			s1		= .
		}
	}
	else {
		// homoskedasticity: initial penalty loadings = vector of 1s after standardization
		Psi		= d.sdvecpnp
		v		= .
		s1		= .
	}

	// initialize lambda
	if (lambda0) {
		// user-supplied lambda0
		lambda=lambda0
	}
	else if ((xdep) & (d.hetero | d.nclust)) {
		// x-dep and inid => initial lambda is simple plug-in rather than simulated
		// 0 arg forces non-xdep.
		lambda=lambdaCalc(d,pminus,gamma,c, R,0,   v,s1,Psi,lalt,newseed,dotsflag,verbose)
	}
	else if (d.hetero | d.nclust) {
		// initial lambda if iterating penalty loadings; uses c0
	 	lambda=lambdaCalc(d,pminus,gamma,c0,R,xdep,v,s1,Psi,lalt,newseed,dotsflag,verbose)
	}
	else {
		// lambda for iid case
	 	lambda=lambdaCalc(d,pminus,gamma,c, R,xdep,v,s1,Psi,lalt,newseed,dotsflag,verbose)
	}

	// initial estimator is ridge (do not mult by 2 as in lasso since 2 not in lambda for sqrt lasso)
	beta_ridge	= anysolver(d.XX+lambda*diag(d.sdvec:^2),d.Xy, r=.)
	if (verbose>=2) {
		printf("{txt}Initial estimator is ridge:\n")
		beta_ridge'
	}

	if ((d.hetero==0) & (d.nclust==0) & (d.bw==0)) {
		// iid. No iteration necessary, even for x-dep.
		if (verbose>=1) {
			printf("Obtaining sqrt-lasso estimate...\n")
		}
		betas	= DoLasso(d, Psi, Psi, lambda, lambda, verbose, optTol, maxIter, zeroTol, alpha, beta_ridge)
		if (verbose>=1) {
			printf("Selected variables: %s\n\n",invtokens(betas.nameXSel'))
		}

	}
	else {
		// heteroskedasticity or clustering or HAC/AC

		// get first lasso estimate based on initial lambda/loadings
		if (verbose>=1) {
			printf("Obtaining initial sqrt-lasso estimate...\n")
		}
		betas	= DoLasso(d, Psi, Psi, lambda, lambda, verbose, optTol, maxIter, zeroTol, alpha, beta_ridge)
		if (verbose>=1) {
			printf("\n")
		}

		do {								// start loop; always iterated at least once
		
			iter++

			if (verbose>=1) {
				printf("Estimation of penalty level/loadings: Step %f.\n",iter)
			}

			// old rmse (will be missing in iteration #1)
			s0		= s1
			// obtain residuals; based on betas(Psi)
			if (lassoPsiflag) {
				v	= betas.v		// lasso residuals
				s1	= betas.rmse 
			}
			else {
				v	= betas.vPL		// post-lasso residuals
				s1	= betas.rmsePL
			}
			// change in RMSE				
			Delta	= abs(s1-s0)

			// Psi update
			// =max(Psi,1), see Alg 1 in Ann of Stat
			// but since X isn't pre-standardized, need to compare to standardization vector instead of unit vector
			// nb: this the the std vector after partialling and with zeros in place of the zero-pen variables
			Psi		= MakeScoreStdVec(d,v,s1)
			Psi		= (Psi :< d.sdvecpnp):*d.sdvecpnp + (Psi :>= d.sdvecpnp):*Psi
			// xdep lambda update
			if (xdep) {
				lambda	= lambdaCalc(d,pminus,gamma,c,R,xdep,v,s1,Psi,lalt,newseed,dotsflag,verbose)
			}
			
			// Reporting
			if (verbose>=1) {
				printf("RMSE: %f\n",s1)
				printf("Change in RMSE: %f\n",Delta)
				printf("Obtaining new sqrt-lasso estimate...\n")
			}

			// new lasso estimate
			betas	= DoLasso(d, Psi, Psi, lambda, lambda, verbose, optTol, maxIter, zeroTol, alpha, beta_ridge)
			if (verbose>=1) {
				printf("Selected variables: %s\n\n",invtokens(betas.nameXSel'))
			}

		} while ((iter < maxPsiIter) & (Delta > PsiTol))
		
		// when iteration concludes, final estimated betas consistent with (have used) final Psi and lambda.

		if (verbose>=1) {
			printf("Number of penalty loading iterations: %g\n",iter)
			if (iter == maxPsiIter) {
				printf("Warning: reached max penalty loading iterations w/o achieving convergence.\n")
			}
			else {
				printf("Penalty loadings (upsilon) convergence achieved.\n")
			}
		}

	}

	// sup-score stat
	if (supscoreflag) {
		betas.supscore	= doSupScore(d, c, gamma, pminus, verbose, ssnumsim, newseed, dotsflag)
	}

	// Misc
	betas.n			= d.n
	betas.nclust	= d.nclust
	betas.nclust1	= d.nclust1
	betas.nclust2	= d.nclust2
	betas.npsiiter	= iter
	betas.sPsi		= betas.Psi :/ d.sdvecpnp		//  should be =1 under homosk.
	betas.sPsi		= editmissing(betas.sPsi,0)		//  in case any unpenalized (div by 0)
	betas.lambda0	= lambda						//  no diff between lambda and lambda0 with sqrt lasso
	betas.slambda	= lambda						//  no diff between lambda and std lambda with sqrt lasso
	betas.c			= c
	betas.gamma		= gamma
	betas.gammad	= gammad
	
	return(betas)

}
// end RSqrtLasso	


void EstimateLassoPath(							///  Complete Mata code for lassopath
				string scalar nameY,			///
				string scalar nameX,			///
				string scalar nameX_o,			///
				string scalar notpen_o,			///
				string scalar notpen_t,			///
				string scalar toest,			///
				string scalar holdout, 			/// validation data
				real scalar cons,				///
				real scalar stdall,				///
				real scalar dmflag,				///
				real scalar dofminus,			///
				real scalar sdofminus,			///
				real scalar pminus,				///
				real scalar prestdflag,			///
				string scalar lambdamat,		/// L1 norm lambda single, list or missing (=> construct default list)
				string scalar lambda2mat,		/// L2 norm lambda (optional)
				real scalar lmax, 				///
				real scalar lcount,				///
				real scalar lminratio,			///
				real scalar lglmnet,			///
				string scalar PsiMat,			/// Optional L1 norm loadings
				string scalar Psi2Mat,			/// Optional L2 norm loadings
				string scalar stdymat,			/// 
				string scalar stdxmat,			///
				real scalar stdl, 				/// standardization loadings?
				real scalar sqrtflag,			/// lasso or sqrt-lasso?
				real scalar alpha,				///
				real scalar adaflag,			/// 
				string scalar AdaMat,			///
				real scalar theta,				/// adaptive lasso/elastic net parameter
				string scalar adafirst,			///
				real scalar post, 				///
				real scalar optTol,				///
				real scalar maxIter,			///
				real scalar zeroTol,			///
				real scalar fdev,				///
				real scalar devmax,				///
				real scalar verbose,			///
				real scalar noic, 				///
				real scalar ebicgamma,			///
				real scalar nodevcrit			///
				)
{

	struct dataStruct scalar d
	d			= MakeData(nameY,nameX,nameX_o,toest,cons,dmflag,dofminus,sdofminus,prestdflag,sqrtflag,stdymat,stdxmat,lglmnet)
	
	// Estimation accommodates pre-standardized data and standardization on-the-fly.
	// Standardization on-the-fly: standardization included in penalty loadings.
	
	if (adaflag==0) {
		// non-adaptive case
		// L1 norm loadings
		if (PsiMat!="") {						//  overall pen loadings supplied
			Psi		= st_matrix(PsiMat)
		}
		else if (stdl) {						//  std loadings - standardize on the fly
			Psi		= d.sdvec					//  L1 norm loadings
		}
		else {									//  (pre-standardized comes here)
			Psi = J(1,d.p,1)					//  default is 1
		}
		// L2 norm loadings - if not supplied, =L1 norm loadings
		if (Psi2Mat!="") {						//  overall L2 norm pen loadings supplied
			Psi2	= st_matrix(Psi2Mat)
		}
		else {									//  default
			Psi2	= Psi
		}
	}
	else {	
		// adaptive lasso/elastic net
		if (AdaMat=="") {
			// no adaptive coefficients provided
			// svd returns the minimum norm generalized solution; no rows/coefs dropped
			// but if x is all zeros (e.g. omitted base var) then coef is zero
			if ((d.n <= d.p) | (adafirst=="univariate")) {
				ada_init		= editmissing(d.Xy :/ diagonal(d.XX),0)'
				if (verbose>=1) {
					printf("{txt}Adaptive weights calculated using univariate OLS regressions.\n")
					r	=rank(d.XX)
					if (r<d.p) {
						printf("{txt}Warning: rank(X'X)=%f < number of regressors=%f (including any factor vars).\n",r,d.p)
					}
				}
				// missing adaptive coefficients set to zero for now
			}
			else {
				ada_init		= anysolver(d.XX/d.n,d.Xy/d.n, r=.,"svd")'
				if (verbose>=1) {
					printf("{txt}Adaptive weights calculated using OLS.\n")
					if (r<d.p) {
						printf("{txt}Warning: rank(X'X)=%f < number of regressors=%f (including any factor vars).\n",r,d.p)
						printf("{text}         generalized inverse used to obtain OLS coefficients.\n")
					}
				}
			}
		}
		else {
			// adaptive coefficients provided
			ada_init		= st_matrix(AdaMat)
		}
		if (verbose>=2) {
			printf("{txt}Initial adaptive coefficients:\n")
			invtokens(strofreal(ada_init))
		}
		ada_init	= ada_init / d.ysd
		if (verbose>=2) {
			printf("{txt}Initial adaptive coefficients after rescaling:\n")
			invtokens(strofreal(ada_init))
		}
		// do inversion; if not prestandardized, std accommodated by exp(theta-1).
		Psi		= ((1:/abs(ada_init)):^theta) :/ (d.sdvec:^(theta-1))
		// if any ada coefs were 0, Psi will be missing; temporarily set these to zero.
		Psi		= editmissing(Psi,0)
		// L2 loadings (not adaptive)
		if (d.prestdflag) {
			Psi2	= J(1,d.p,1)
		}
		else {
			Psi2	= d.sdvec
		}
	}

	//  need to set loadings of notpen vars = 0
	if (notpen_o~="") {	
		npnames=tokens(notpen_o)
		forbound = cols(npnames)	//  faster
		for (i=1; i<=forbound; i++) {
				Psi		=Psi  :* (1:-(d.nameX_o:==npnames[1,i]))
				Psi2	=Psi2 :* (1:-(d.nameX_o:==npnames[1,i]))
		}
	}

	// lglmnet parameterization - penalty loadings sum to p
	// note this accommodates zero penalty loadings by adjusting others upwards
	if (lglmnet) {
		if (verbose>=2) {
			printf("{txt}Rescaling loadings to sum to p (lglmnet convention).\n")
		}
		Psi		= Psi  * cols(Psi)  * 1/rowsum(Psi)
		Psi2	= Psi2 * cols(Psi2) * 1/rowsum(Psi2)
	}

	if (adaflag) {
		// Ada zero coefs should get a near-inf penalty after inversion.
		// Exception is when variable is all zeros (e.g. base of a factor var).
		ada_zero	= (d.mvec :== 0) :& (d.sdvec :== 0)
		ada_inf		= (ada_init :== 0) :& (ada_zero :== 0)
		ada_other	= 1 :- ada_zero :- ada_inf
		// max Stata float is 1.7 x 10^38 (nb: Mata's maxdouble() is 8.9885e+307)
		Psi			= (ada_other :* editmissing(Psi,0)) + (ada_inf) * 1e38
		if (verbose>=2) {
			printf("{txt}Adaptive weights after inversion:\n")
			invtokens(strofreal(Psi))
		}
	}

	if (lambdamat=="") {
		// no lambdas provided
		if (d.lglmnet) {
			if (lmax<=0) { // no lambda max given
				// see Friedman et al (J of Stats Software, 2010)  
				if (d.prestdflag) {
					lmax = max(abs(d.Xy))/d.n * 1/max((0.001,alpha))
				}
				else {
					lmax = max(abs((d.Xy)):/((Psi)'))/d.ysd/d.n * 1/max((0.001,alpha))
				}
			}
			else {		// lambda max provided
				// if lambda max has been provided in the standardized metric, then must first unstandardize
				if (stdall) {
					if (d.prestdflag) {
						lmax	=lmax * d.prestdy
					}
					else {
						lmax	=lmax * d.ysd
					}
				}
				if (d.prestdflag) {
						// lambda max provided, rescale it
						lmax = lmax / d.prestdy
				}
				else {
						// lambda max provided, rescale it
						lmax = lmax / d.ysd
				}
			}
			lmin	= lminratio*lmax
			lambda	=exp(rangen(log(lmax),log(lmin),lcount))'
		}
		else {
			// not glmnet
			if (lmax<=0) { // no lambda max given
				if (d.sqrtflag) {
					// sqrt lambda grid should be invariant to the scaling of y
					lmax = max(abs((d.Xy)):/((Psi)')) * 1/d.ysd
				}
				else {
					// see Friedman et al (J of Stats Software, 2010)  
					lmax = max(abs((d.Xy)):/((Psi)'))*2/max((0.001,alpha))
				}
			}
			else {		// lambda max provided
				// if lambda max has been provided in the standardized metric, then must first unstandardize
				if (stdall) {
					if (d.prestdflag) {
						lmax	=lmax * d.prestdy
					}
					else {
						lmax	=lmax * d.ysd
					}
				}
				if ((d.prestdflag) & (!d.sqrtflag)) {
						lmax = lmax / d.prestdy
				}
			}
			lmin = lminratio*lmax
			lambda=exp(rangen(log(lmax),log(lmin),lcount))'
		}
	}
	else {
		// lambdas provided
		lambda			=st_matrix(lambdamat)
		// if lambdas have been provided in the standardized metric, then must first unstandardize
		if (stdall) {
			if (d.prestdflag) {
				lambda	=lambda * d.prestdy
			}
			else {
				lambda	=lambda * d.ysd
			}
		}
		// now adjust lambdas for prestandardization as necessary
		if (d.lglmnet) {
			if (d.prestdflag) {
				lambda	=lambda * 1/(d.prestdy)
			}
			else {
				lambda	=lambda * 1/(d.ysd)
			}
		}
		else if ((d.prestdflag) & (!d.sqrtflag)) {		//  data have been pre-standardized, so adjust lambdas accordingly
			lambda		=lambda * 1/(d.prestdy)
		}
	}
	if (lambda2mat=="") {
		// no separate L2 norm lambda provided => default is L1 norm lambda
		lambda2			= lambda
	}
	else {
		// separate L2 norm lambda provided
		lambda2			=st_matrix(lambda2mat)
		// if lambdas have been provided in the standardized metric, then must first unstandardize
		if (stdall) {
			if (d.prestdflag) {
				lambda2	=lambda2 * d.prestdy
			}
			else {
				lambda2	=lambda2 * d.ysd
			}
		}
		// now adjust lambdas for prestandardization as necessary
		if (d.lglmnet) {
			if (d.prestdflag) {
				lambda2	=lambda2 * 1/(d.prestdy)
			}
			else {
				lambda2	=lambda2 * 1/(d.ysd)
			}
		}
		else if ((d.prestdflag) & (!d.sqrtflag)) {		//  data have been pre-standardized, so adjust lambdas accordingly
			lambda2		=lambda2 * 1/(d.prestdy)
		}
	}

	// check dimensions
	if (cols(lambda)!=cols(lambda2)) {
		errprintf("error - dimensions of L1 and L2 penalties must be the same\n")
		exit(198)
	}

	if ((cols(lambda)==1) & (!hasmissing(lambda))) {				//  one lambda
		struct outputStruct scalar OUT
		OUT = DoLasso(d,Psi,Psi2,lambda,lambda2,verbose,optTol,maxIter,zeroTol,alpha)
		ReturnResults(OUT,d)
	}
	else if ((cols(lambda)>1) & (!hasmissing(lambda))) {		//  lambda is a vector or missing (=> default list)
		struct outputStructPath scalar OUTPATH
		OUTPATH = DoLassoPath(d,Psi,Psi2,lambda,lambda2,post,verbose,optTol,maxIter,zeroTol,fdev,devmax,alpha,lglmnet,nodevcrit)
		ReturnResultsPath(OUTPATH,d,nameX_o)
		if (holdout!="") { // used for cross-validation
			getMSPE(OUTPATH,nameY,nameX,holdout,d)  
		}
		else if (!noic) { // calculate IC 
			getInfoCriteria(OUTPATH,d,ebicgamma,pminus)
		}
	}
}
// end EstimateLassoPath

struct outputStruct scalar DoLasso(								///
							struct dataStruct scalar d,			/// #1  data (y,X)
							real rowvector Psi,					/// #2  penalty loading vector (L1 norm)
							real rowvector Psi2,				/// #3  penalty loading vector (L2 norm)
							real scalar lambda,					/// #4  lambda (single value) (L1 norm)
							real scalar lambda2,				/// #5  lambda (single value) (L2 norm)
							real scalar verbose,				/// #6  reporting
							real scalar optTol,					/// #7  convergence of beta estimates
							real scalar maxIter,				/// #8  max number of shooting iterations
							real scalar zeroTol, 				/// #9  tolerance to set coefficient estimates to zero
							real scalar alpha, 					/// #10 elastic net parameter
							|									///
							real colvector beta_init			/// #11 initial beta vector (optional)
							)
{

	struct outputStruct scalar	t
	struct betaStruct scalar	b

	// no initial beta supplied
	noinitflag		= args()<11

	// no initial beta supplied (optional last argument) => will be initialized in DoShooting(.)
	// get estimates
	if (noinitflag) {
		b		= DoShooting(d,Psi,Psi2,lambda,lambda2,verbose,optTol,maxIter,zeroTol,alpha)
	}
	else {
		b		= DoShooting(d,Psi,Psi2,lambda,lambda2,verbose,optTol,maxIter,zeroTol,alpha,beta_init)
	}
	beta	= b.beta
	m		= b.m
	fobj	= b.fobj
	dev		= b.dev

	// save results in t struct
	t.niter = m

	// full vector
	t.betaAll = beta

	// compare initial beta vs estimated beta
	//(beta,beta_init)

	// should set coefs to zero IF NOT RIDGE
	// should compare STANDARDIZED coefs to zeroTol
	if (alpha==0) {
		// ridge, all selected
		t.index = beta*0 :+ 1
		s = rows(beta)
	}
	else if (d.prestdflag) {
		// prestandardized, so compare directly to zeroTol
		t.index = abs(beta) :> zeroTol
		s = sum(abs(beta) :> zeroTol) // number of selected vars, =0 if no var selected
	}
	else {
		// compare standardized coefs to zeroTol
		t.index = abs(beta :* d.sdvec' / d.ysd) :> zeroTol
		s = sum(abs(beta :* d.sdvec' / d.ysd) :> zeroTol)	
	}

	// if ridge, keep entire beta vector as-is
	// otherwise reduce beta vector to only non-zero coeffs
	if (alpha==0) {
		// ridge
		t.beta = beta
	}
	else if (s>0) {
		t.beta = select(beta,t.index)
	}
	else {
		t.beta = .
	}

	// obtain post-OLS estimates
	if ((s>0) & (s<d.n) & (d.cons==0)) {
		// data are zero mean or there is no constant
		betaPL			= qrsolve(select(*d.X,t.index'),*d.y)
	}
	else if ((s>0) & (s<d.n)) {
		// model has a constant
		betaPL			= qrsolve((select(*d.X,t.index'):-select(d.mvec,t.index')),(*d.y:-d.ymvec))
	}
	else if (s>0) {
		betaPL			= J(s,1,.)		// set post-OLS vector = missing if s-hat > n. 
	}
	else if (s==0) {
		betaPL			= .
	}
	t.betaPL = betaPL
	
	// intercept handled in ReturnResults() to accommodate standardization etc.
	
	// "All" version of PL coefs
	t.betaAllPL = J(rows(beta),1,0)
	// need to look out for missing betaPL; if s=0 betaAllPL will just be a vector of zeros
	if (s>0) {
		t.betaAllPL[selectindex(t.index),1] = betaPL
	}

	// other objects
	if (s>0) {
		t.nameXSel	= select(d.nameX_o',t.index)
	}
	else {
		t.nameXSel	=""
	}
	t.Psi		= Psi
	t.lambda	= lambda
	t.cons		= d.cons
	t.n			= d.n
	t.beta_init	= beta_init
	t.s 		= s

	// obtain residuals
	if ((s>0) & (d.cons==0)) {
		// data are zero mean or there is no constant
		t.v		= *d.y - select(*d.X,t.index')*t.beta
		t.vPL	= *d.y - select(*d.X,t.index')*t.betaPL
	}
	else if (s>0) {
		// model has a constant
		t.v		= (*d.y:-d.ymvec) - (select(*d.X,t.index'):-select(d.mvec,t.index'))*t.beta
		t.vPL	= (*d.y:-d.ymvec) - (select(*d.X,t.index'):-select(d.mvec,t.index'))*t.betaPL
	}
	else if (d.cons==0) {
		// data are zero mean or no constant in model; nothing selected; residual is just y
		t.v		= *d.y
		t.vPL	= *d.y
	}
	else {
		// nothing selected; residual is demeaned y
		t.v		= (*d.y:-d.ymvec) 
		t.vPL	= (*d.y:-d.ymvec)
	}
	
	// RMSE
	t.rmse		=sqrt(mean(t.v:^2))
	t.rmsePL	=sqrt(mean(t.vPL:^2))
	
	// R-squared. Same as "dev" (deviation)
	t.r2		=dev
	
	// rescaling minimized objective function
	if (d.lglmnet) {
		// use glmnet definitions
		// sqrt-lasso not supported by lglmnet option
		if (alpha==1) {
			// lasso
			t.objfn		= 1/2*(t.rmse)^2 + lambda*quadcross(Psi',abs(beta))*d.ysd
		}
		else if (alpha==0) {
			// ridge
			t.objfn		= 1/2*(t.rmse)^2 + 1/2*lambda2*quadcross(Psi2':^2,beta:^2)
		}
		else {
			// elastic net
			t.objfn		= 1/2*(t.rmse)^2												///
							+     lambda * alpha   *quadcross(Psi',abs(beta))*d.ysd			///
							+ 1/2*lambda2*(1-alpha)*quadcross(Psi2':^2,beta:^2)
		}
	}
	else {
		if (d.sqrtflag) {
			// sqrt-lasso - no rescaling needed
			t.objfn		= fobj
		}
		else if (alpha==1) {
			// lasso
			t.objfn		= (t.rmse)^2 + lambda/d.n * quadcross(Psi',abs(beta))
		}
		else if (alpha==0) {
			// ridge
			t.objfn		= (t.rmse)^2													///
							+ 1/2*lambda2*d.prestdy/d.n * quadcross((Psi2'):^2,beta:^2)
		}
		else {
			// elastic net
			t.objfn		= (t.rmse)^2														///
							+     lambda/d.n           * alpha   *quadcross(Psi',abs(beta))	///
							+ 1/2*lambda2*d.prestdy/d.n*(1-alpha)*quadcross((Psi2'):^2,beta:^2)
		}
	}

	if (verbose>=1) {
		printf("{txt}Minimized objective function: {res}%f\n",t.objfn)
	}

	// effective degrees of freedom
	t.dof			= getdf(d,beta,Psi,lambda,alpha,verbose)
	if (verbose>=1) {
		printf("{txt}Effective degrees of freedom: {res}%f\n",t.dof)
	}

	return(t)
}
// end DoLasso


void ReturnResults(		struct outputStruct scalar t,		/// #1
						struct dataStruct scalar d			/// #2
						)
{

	if (rows(t.betaAll)) {					// estimation results to insert

		// initialize from results struct
		s			= t.s					// all vars in param vector incl notpen but EXCL constant
		k			= t.s+t.cons			// all vars in param vector incl notpen PLUS constant
		betaAll		= t.betaAll
		betaAllPL	= t.betaAllPL
		Psi			= t.Psi
		ePsi		= t.Psi					// rlasso only
		sPsi		= t.sPsi
		rmse		= t.rmse
		rmsePL		= t.rmsePL
		r2			= t.r2
		lambda		= t.lambda
		objfn		= t.objfn
		// initialize from data struct
		AllNames	= d.nameX_o'			// d.nameX_o is a row vector; AllNames is a col vector (to match coef vectors)

		// standardized or unstandardize
		if (d.lglmnet) {
			if (d.prestdflag) {
				sbetaAll		= betaAll
				sbetaAllPL		= betaAllPL
				betaAll			= betaAll		:/ d.prestdx' * d.prestdy
				betaAllPL		= betaAllPL		:/ d.prestdx' * d.prestdy
				slambda			= lambda
				srmse			= rmse
				srmsePL			= rmsePL
				sobjfn			= objfn
				lambda			= lambda		* d.prestdy
				rmse			= rmse			* d.prestdy
				rmsePL			= rmsePL		* d.prestdy
				// objfn is pmse
				objfn			= objfn			* (d.prestdy)^2
			}
			else {
				sbetaAll		= betaAll		:* d.sdvec' / d.ysd
				sbetaAllPL		= betaAllPL		:* d.sdvec' / d.ysd
				slambda			= lambda
				lambda			= lambda		* d.ysd
				srmse			= rmse			/ d.ysd
				srmsePL			= rmsePL		/ d.ysd
				sobjfn			= objfn			/ (d.ysd)^2
			}
		}
		else if (d.prestdflag) {
			// standard case (no lglmnet)
			sbetaAll		= betaAll
			sbetaAllPL		= betaAllPL
			//  beta is col vector, prestdx is row vector, prestdy is scalar
			betaAll			= betaAll		:/ d.prestdx' * d.prestdy
			betaAllPL		= betaAllPL		:/ d.prestdx' * d.prestdy
			slambda			= lambda
			srmse			= rmse
			srmsePL			= rmsePL
			sobjfn			= objfn
			ePsi			= ePsi			:* d.prestdx				//  rlasso only; ePsi and prestsdx both row vectors; Psi does not change
			rmse			= rmse			* d.prestdy
			rmsePL			= rmsePL		* d.prestdy
			if (d.sqrtflag==0) {
				// lasso => objfn is pmse
				objfn		= objfn			* (d.prestdy)^2
			}
			else {
				// sqrt lasso => objfn is prmse
				objfn		= objfn			* d.prestdy
			}
			if (d.sqrtflag==0) {									//  sqrt-lasso lambdas don't need unstandardizing
				lambda		= lambda		* d.prestdy
			}
		}
		else {
			// standard case (no lglmnet), no prestd
			sbetaAll		= betaAll		:* d.sdvec' / d.ysd
			sbetaAllPL		= betaAllPL		:* d.sdvec' / d.ysd
			slambda			= lambda		/ d.ysd
			srmse			= rmse			/ d.ysd
			srmsePL			= rmsePL		/ d.ysd
			sobjfn			= objfn			/ (d.ysd)^2
		}
	
		// get intercept
		if (t.cons) {											//  pre-standardized means no constant
			intercept	= mean(*d.y) :- mean(*d.X)*betaAll
			interceptPL	= mean(*d.y) :- mean(*d.X)*betaAllPL
		}

		if (s>0) {												//  do here so that we select from std/unstd vector
			beta		= select(betaAll,t.index)
			betaPL		= select(betaAllPL,t.index)
			sbeta		= select(sbetaAll,t.index)
			sbetaPL		= select(sbetaAllPL,t.index)
			Names		= select(AllNames,t.index)
		}
		else {
			beta		= .
			betaPL		= .
			sbeta		= .
			sbetaPL		= .
			Names		= ""
		}
	
		if ((s>0) & (t.cons)) {									//  add constant to end of vectors (unstandardized only)
			beta		= (beta			\ intercept)		
			betaPL		= (betaPL		\ interceptPL)	
			betaAll		= (betaAll		\ intercept)		
			betaAllPL	= (betaAllPL	\ interceptPL)	
			NamesCons	= (Names		\ "_cons")
			AllNamesCons= (AllNames 	\ "_cons")
		}
		else if ((s>0) & (!t.cons)) {							//  coef vectors already OK, just need names with _cons
			NamesCons	= Names
			AllNamesCons= AllNames
		}
		else if ((s==0) & (t.cons)) {
			beta		= intercept						
			betaPL		= interceptPL					
			NamesCons	= "_cons"
			AllNamesCons= (AllNames 	\ "_cons")
			betaAll		= (betaAll		\ intercept)			//  will be all zeros + intercept at end
			betaAllPL	= (betaAllPL	\ interceptPL)	
		}
		else {
			beta		= .
			betaPL		= .
			sbeta		= .
			sbetaPL		= .
			NamesCons	= ""
			AllNamesCons= AllNames
		}

		st_rclear() 
		if ((s > 0) | (t.cons)) {
			// matrix stripes
			coln		=(J(rows(Names),1,""),Names)			//  Names is a col vector so need #rows
			colnCons	=(J(rows(NamesCons),1,""),NamesCons)	//  Names is a col vector so need #rows
			// row vector of names
			st_global("r(sel)",invtokens(Names'))				//  selected vars exclude cons but here may include notpen
			// column vectors
			st_matrix("r(b)",beta)
			st_matrix("r(bOLS)",betaPL)
			st_matrixrowstripe("r(b)",colnCons)
			st_matrixrowstripe("r(bOLS)",colnCons)
			if (s>0) {
				// standardized coefs exclude constant
				st_matrix("r(sb)",sbeta)
				st_matrix("r(sbOLS)",sbetaPL)
				st_matrixrowstripe("r(sb)",coln)
				st_matrixrowstripe("r(sbOLS)",coln)
			}
		}
		
		// matrix stripe
		AllcolnCons=(J(rows(AllNamesCons),1,""),AllNamesCons)
		Allcoln    =(J(rows(AllNames),1,""),    AllNames)
		// column vectors
		st_matrix("r(bAll)",betaAll)
		st_matrix("r(bAllOLS)",betaAllPL)
		st_matrixrowstripe("r(bAll)",AllcolnCons)
		st_matrixrowstripe("r(bAllOLS)",AllcolnCons)
		// standardized coefs exclude constant
		st_matrix("r(sbAll)",sbetaAll)
		st_matrix("r(sbAllOLS)",sbetaAllPL)
		st_matrixrowstripe("r(sbAll)",Allcoln)
		st_matrixrowstripe("r(sbAllOLS)",Allcoln)
	
		// matrix stripe
		coln=(J(rows(AllNames),1,""),AllNames)
		// row vectors
		st_matrix("r(stdvec)",d.sdvec)
		st_matrix("r(Psi)",Psi)
		st_matrix("r(ePsi)",ePsi)						//  rlasso only
		st_matrixcolstripe("r(stdvec)",coln)
		st_matrixcolstripe("r(Psi)",coln)
		st_matrixcolstripe("r(ePsi)",coln)				//  rlasso only
	
		if ((cols(sPsi)>0) & (cols(sPsi)<.)) {			//  "cols(t.sPsi)<." may be unnecessary - if missing, cols(.) = 0.
			st_matrix("r(sPsi)",sPsi)
			st_matrix("r(stdvecpnp)",d.sdvecpnp)
			st_matrixcolstripe("r(sPsi)",coln)
			st_matrixcolstripe("r(stdvecpnp)",coln)
		}
		
		st_matrix("r(beta_init)",t.beta_init)
	}

	// BCH stat, p-value, crit-value, signif; rlasso stat, p-value
	if (cols(t.supscore)) {
		st_matrix("r(supscore)",t.supscore)
		coln=(J(4,1,""),("CCK_ss" \ "CCK_p" \ "CCK_cv" \ "CCK_gamma"))
		st_matrixcolstripe("r(supscore)",coln)
	}
	// rlasso has its own standardized lambda
	if (t.slambda<.) {
		slambda	=t.slambda
	}

	// Can always return these; scalars will just be missing
	st_numscalar("r(lambda)",lambda)
	st_numscalar("r(rmse)",rmse)
	st_numscalar("r(rmsePL)",rmsePL)
	st_numscalar("r(srmse)",srmse)
	st_numscalar("r(srmsePL)",srmsePL)
	st_numscalar("r(r2)",r2)
	st_numscalar("r(objfn)",objfn)
	st_numscalar("r(sobjfn)",sobjfn)
	st_numscalar("r(s)",s)
	st_numscalar("r(k)",k)
	st_numscalar("r(lcount)",1)
	st_numscalar("r(lambda0)",t.lambda0)
	st_numscalar("r(slambda)",slambda)
	st_numscalar("r(c)",t.c)
	st_numscalar("r(gamma)",t.gamma)
	st_numscalar("r(gammad)",t.gammad)
	st_numscalar("r(niter)",t.niter)
	st_numscalar("r(npsiiter)",t.npsiiter)
	st_numscalar("r(dof)",t.dof)

	st_numscalar("r(N)",d.n)
	st_numscalar("r(N_clust)",d.nclust)
	st_numscalar("r(N_clust1)",d.nclust1)
	st_numscalar("r(N_clust2)",d.nclust2)
	
	// for HAC and 2-way cluster lasso (neg loadings possible)
	st_numscalar("r(psinegs)",d.psinegs)
	st_global("r(psinegvars)",d.psinegvars)
	
}
// end ReturnResults

real matrix getMinIC(real matrix IC)		//  =0 if nothing to partial out, =projection coefs if Zs.
{
		if (nonmissing(IC)) {		// at least one non-missing IC value
			licid=.
			minindex(IC,1,licid,.)	// returns index of lambda that minimises ic
			if (rows(licid)>1) {    // no unique lopt 
				licid=licid[1,1] 	
				icmin=IC[1,licid]
				licunique=0
			}
			else {
				icmin=IC[1,licid]
				licunique=1
			}
		}
		// address issue if all values are missing
		else {
			licid=1
			icmin=.
			licunique=0
		}

	R = (licid,icmin,licunique)

	return(R)
}
// end get minimum IC

void getInfoCriteria(struct outputStructPath scalar t,
 			struct dataStruct scalar d,
			real ebicgamma,
			real pminus)
{		
		// t.betas is lcount by p	
		// t.dof is 1 by lcount
		// need to check for constant
		// pminus is #missing vars (base or omitted variables)
		XB  = quadcross(((*d.X):-(d.mvec*d.cons))',(t.betas)')	// n by lcount

		TSS = quadcolsum(((*d.y):-(d.ymvec*d.cons)):^2)	// 1 by lcount
		ESS = quadcolsum(((XB) :-(d.ymvec*d.cons)):^2)
		// need to demean y here as well
		RSS = quadcolsum((((*d.y):-(d.ymvec*d.cons)) :-(XB)):^2)
		if (d.prestdflag) {	
			TSS=TSS*(d.prestdy)^2
			ESS=ESS*(d.prestdy)^2
			RSS=RSS*(d.prestdy)^2
		}
		RSQ =1:-RSS:/TSS
		// ebic parameter
		// default choice is based on P=n^kappa and gamma=1-1/(2*kappa)
		// see Chen & Chen (2008, p. 768, Section 5)
		// note we subtract pminus to account for non-variables (all zeros = base or omitted vars) 
		if ((ebicgamma<0) | (ebicgamma>1)) {
			ebicgamma = 1-log(d.n)/(2*log(d.p-pminus))
			ebicgamma = max((ebicgamma,0)) // ensures that ebicgamma are in [0,1]
			ebicgamma = min((ebicgamma,1))
		}

		// calculate aic and bic
		// For reference: in R the definitions would be
		// AIC = d.n + d.n*log(2*pi()) + d.n*log(RSS/d.n) + 2*(t.dof')
		// BIC = d.n + d.n*log(2*pi()) + d.n*log(RSS/d.n) + log(d.n)*(t.dof')
		AIC		= d.n*log(RSS/d.n) + (t.dof')*2 
		BIC 	= d.n*log(RSS/d.n) + (t.dof')*log(d.n)
		// note we subtract pminus to account for non-variables (all zeros = base or omitted vars) 
		EBIC 	= BIC :+ 2 * (t.dof') * log(d.p-pminus) * ebicgamma
		AICC	= d.n*log(RSS/d.n) + (t.dof')*2:*((d.n):/(d.n:-t.dof'))

		// obtain minimum IC and obtimal lambda
		AICinfo =  getMinIC(AIC)
		BICinfo = getMinIC(BIC)
		EBICinfo = getMinIC(EBIC)
		AICCinfo = getMinIC(AICC)
		laicid=AICinfo[1,1]
		aicmin=AICinfo[1,2]
		lbicid=BICinfo[1,1]
		bicmin=BICinfo[1,2]
		lebicid=EBICinfo[1,1]
		ebicmin=EBICinfo[1,2]
		laiccid=AICCinfo[1,1]
		aiccmin=AICCinfo[1,2]
		
		st_matrix("r(dof)",t.dof)
		st_matrixcolstripe("r(dof)",("","Df"))
		// return 
 		st_matrix("r(rsq)",RSQ')
 		st_matrix("r(tss)",TSS')
		st_matrix("r(ess)",ESS')
		st_matrix("r(rss)",RSS')
		/// aic
		st_matrix("r(aic)",AIC')
		st_numscalar("r(aicmin)",aicmin)
		st_numscalar("r(laicid)",laicid)
		/// bic
		st_matrix("r(bic)",BIC')
		st_numscalar("r(bicmin)",bicmin)
		st_numscalar("r(lbicid)",lbicid)	
		/// ebic
		st_matrix("r(ebic)",EBIC')
		st_numscalar("r(ebicmin)",ebicmin)
		st_numscalar("r(lebicid)",lebicid)
		st_numscalar("r(ebicgamma)",ebicgamma)
		/// aicc
		st_matrix("r(aicc)",AICC')
		st_numscalar("r(aiccmin)",aiccmin)
		st_numscalar("r(laiccid)",laiccid)
		// constant to remove from IC to standardize
		if (d.prestdflag) {	
			icstd	= 2*ln((d.prestdy))
		}
		else {
			icstd	= 2*ln((d.ysd))
		}
		st_numscalar("r(icstd)",icstd)
		st_matrix("r(saic)",AIC' :- icstd)
		st_matrix("r(sbic)",BIC' :- icstd)
		st_matrix("r(sebic)",EBIC' :- icstd)
		st_matrix("r(saicc)",AICC' :- icstd)
		st_numscalar("r(saicmin)",aicmin - icstd)
		st_numscalar("r(sbicmin)",bicmin - icstd)
		st_numscalar("r(sebicmin)",ebicmin - icstd)
		st_numscalar("r(saiccmin)",aiccmin - icstd)
}
// end 


void getMSPE(struct outputStructPath scalar t,
								string scalar varY,
								string scalar varX, 
								string scalar holdout, // marks validation data set
								struct dataStruct scalar d)
{		
		// get beta matrix
		bhat=t.betas 		// lcount by p	
	
		// get validation data
		st_view(y0,.,varY,holdout)
		st_view(X0,.,varX,holdout) 	// n by p
		
		// predicted values
		X0B=quadcross(X0',bhat') 	// n by lcount

		// add intercepts
		if (t.cons) { 	
			X0B=X0B :+ t.intercept 	// t.intercept is 1 by lcount
		}
	
		// mean squared prediction error
		//X0
		MSPE= mean((y0:-X0B):^2) 	// 1 by lcount vector
		
		if (d.prestdflag) (
			MSPE = MSPE :* (d.prestdy)^2
		)
		
		st_matrix("r(mspe)",MSPE)
}
// end getMSPE


// return the lambda id of the largest lambda at which the MSE
// is within one standard error of the minimal MSE.		
real scalar getOneSeLam(real rowvector mse, real rowvector sd, real scalar id) 
{
	minmse = mse[1,id] // minimal MSE
	minsd = sd[1,id]  // SE of minimal MSE
	criteria = mse[1,id]+sd[1,id] // max allowed MSE
	for (j=0; j<id; j++) {
		theid=id-j 
		thismspe= mse[1,theid]
		if (thismspe > criteria) { // if MSE is outside of interval, stop
				theid = id-j+1 // go back by one id and break
				break
		}
	} 
	return(theid)
}

// ********* s_vkernel (modified from livreg2) ******************* //
// Program checks whether kernel and bw choices are valid.
// s_vkernel is called from Stata.
// Arguments are the kernel name (req) and bandwidth (req), both strings.
// Returns results in r() macros.
// r(kernel) - name of kernel (string)
// r(spectral) - =1 if a spectral kernel (using all lags), =0 otherwise

void s_vkernel(string scalar kernel, string scalar bwstring)
{

	// Check bandwidth
	bw=strtoreal(bwstring)
	if (bw==.) {
		errprintf("bandwidth option bw() required for HAC-robust estimation\n")
		exit(102)
	}
	if (bw<=0) {
		errprintf("invalid bandwidth in option bw() - must be real > 0\n")
		exit(198)
	}

	// Check kernel
	// Valid kernel list is abbrev, full name, whether special case if bw=1
	// First in list is default kernel = Barlett
	vklist = 	(	("", "bartlett", "0")
				\	("bar", "bartlett", "0")
				\	("bartlett", "bartlett", "0")
				\	("par", "parzen", "0")
				\	("parzen", "parzen", "0")
				\	("tru", "truncated", "1")
				\	("truncated", "truncated", "1")
				\	("thann", "tukey-hanning", "0")
				\	("tukey-hanning", "tukey-hanning", "0")
				\	("thamm", "tukey-hamming", "0")
				\	("tukey-hamming", "tukey-hamming", "0")
				\	("qua", "quadratic spectral", "1")
				\	("qs", "quadratic spectral", "1")
				\	("quadratic-spectral", "quadratic spectral", "1")
				\	("quadratic spectral", "quadratic spectral", "1")
				\	("dan", "danielle", "1")
				\	("danielle", "danielle", "1")
				\	("ten", "tent", "1")
				\	("tent", "tent", "1")
			)
	kname=strltrim(strlower(kernel))
	pos = (vklist[.,1] :== kname)
	vkname		= strproper(select(vklist[.,2],pos))

	// Exit with error if not in list
	if (sum(pos)==0) {
		errprintf("invalid kernel\n")
		exit(198)
	}
	// Warn if kernel is type where bw=1 means no lags are used
	if (bw==1 & select(vklist[.,3],pos)=="0") {
		printf("{result}Note: kernel=%s", vkname)
		printf("{result} and bw=1 implies zero lags used.\n")
	}

	spectral	= ((vkname=="Quadratic Spectral") | (vkname=="Danielle") | (vkname=="Tent"))

	st_global("r(kernel)", vkname)
	st_numscalar("r(spectral)",spectral)

}  // end of program s_vkernel


// kernel weight (from livreg2)
real scalar m_calckw(	real scalar tau,
						real scalar bw,
						string scalar kernel) 
	{
				karg = tau / bw
				if (kernel=="Truncated") {
					kw=1
				}
				if (kernel=="Bartlett") {
					kw=(1-karg)
				}
				if (kernel=="Parzen") {
					if (karg <= 0.5) {
						kw = 1-6*karg^2+6*karg^3
					}
					else {
						kw = 2*(1-karg)^3
					}
				}
				if (kernel=="Tukey-Hanning") {
					kw=0.5+0.5*cos(pi()*karg)
				}
				if (kernel=="Tukey-Hamming") {
					kw=0.54+0.46*cos(pi()*karg)
				}
				if (kernel=="Tent") {
					kw=2*(1-cos(tau*karg)) / (karg^2)
				}
				if (kernel=="Danielle") {
					kw=sin(pi()*karg) / (pi()*karg)
				}
				if (kernel=="Quadratic Spectral") {
					kw=25/(12*pi()^2*karg^2) /*
						*/ * ( sin(6*pi()*karg/5)/(6*pi()*karg/5) /*
						*/     - cos(6*pi()*karg/5) )
				}
				return(kw)
	}  // end kw


// partial out program
void s_partial(	string scalar Ynames,
				string scalar Pnames,
				string scalar tYnames,
				string scalar touse,		//  indicates full sample
				string scalar toest,		//  indicates estimation subsample
				string scalar wvar,			//  optional weight variable
				scalar dmflag,				//  indicates data are already mean zero; no constant
				string scalar psolver)

{

// All varnames should be basic form, no FV or TS operators etc.
// Y = structural variables
// P = variables to be partialled out
// W = weight variable normalized to mean=1; vector of 1s if no weights
// touse = sample
// dmflag = 0 or 1
// psolver = 0, 1 or 2
// Strategy is to demean (numerically more stable in case of scaling problems)
// and then use A Mata solver.

	Ytokens=tokens(Ynames)
	Ptokens=tokens(Pnames)
	st_view(Y, ., Ytokens, touse)			//  full sample
	st_view(P, ., Ptokens, touse)
	st_view(Yest, ., Ytokens, toest)		//  estimation subsample
	st_view(Pest, ., Ptokens, toest)

	if (tYnames ~= "") {
		tYflag=1
		tYtokens=tokens(tYnames)
		st_view(tY, ., tYtokens, touse)		//  full sample
	}
	else {
		tYflag=0
	}
	
	if (wvar ~= "") {						//  weight variable provided
		st_view(W, ., wvar, touse)
	}
	else {									//  no weight variable => W is vector of 1s
		W = J(rows(Yest),1,1)
	}

	L = cols(P)

	// means are based on estimation subsample
	// weighted means using weight vector W; W is a vector of ones if unweighted
	// if dmflag=1, everything is already mean zero and there is no constant so means unneeded
	if ((!dmflag) & L>0) {					//  Vars to partial out including constant
		Ymeans = mean(Yest, W)
		Pmeans = mean(Pest, W)
	}
	else if (!dmflag) {						//  Only constant to partial out = demean
		Ymeans = mean(Yest, W)
	}

	//	Partial-out coeffs, incorporating weights in the projection matrix.
	//	Not necessary if no vars other than constant. r=rank.
	//  coef vector b is based on estimation subsample
	if ((!dmflag) & L>0) {
		b = anysolver((Pest :- Pmeans):*sqrt(W), (Yest :- Ymeans):*sqrt(W), r=., psolver)	//  partial out P + cons
	}
	else if (L>0) {
		b = anysolver(Pest:*sqrt(W), Yest:*sqrt(W), r=., psolver)							//  partial out P
	}
	else {
		r=1																					//  partial out cons (demean)
	}
	
	if (r==.) {
		_error(3352, "Singular matrix encountered; partialling out failed.")
	}
	else if (r<L) {
		printf("{text}Warning: reduced rank encountered when partialling out (may be caused by omitted vars).\n")
	}

	//	Replace with residuals.
	//  Use full sample (Y, P) and not estimation subsample (Yest, Pest).
	//  Use beta obtained with weighting but do not transform data.
	if ((!dmflag) & L>0) {					//  Vars to partial out including constant
		if (tYflag) {
			tY[.,.] = (Y :- Ymeans) - (P :- Pmeans)*b
		}
		else {
			Y[.,.] = (Y :- Ymeans) - (P :- Pmeans)*b
		}
	}
	else if (L>0) {							//  Vars to partial out NOT including constant
		if (tYflag) {
			tY[.,.] = Y - P*b
		}
		else {
			Y[.,.] = Y - P*b
		}

	}
	else {									//  Only constant to partial out = demean
		if (tYflag) {
			tY[.,.] = Y :- Ymeans
		}
		else {
			Y[.,.] = Y :- Ymeans
		}
	}
	
	// Return list of tokens of dropped collinear vars, if any
	if (r<cols(P)) {								//  something was dropped
		dsel = (b :== 0)							//  matrix, cols=#Y, rows=L, =1 if dropped
		dsel = (rowsum(dsel) :== cols(Y))			//  col vector, =1 if always dropped
		dlist = invtokens(select(Ptokens, dsel'))	//  string of names (tokens)
	}
	else {
		dlist = ""									//  empty list
	}
	st_global("r(dlist)",dlist)						//	Return list of dropped collinears.
	st_numscalar("r(rank)",r)						//  Don't include constant in rank

}  
//end program s_partial

// Mata utility to mimic panelsum(.) with unsorted data
function uspanelsum (	numeric matrix A,			///
						numeric matrix pid,			///
						|							///
						real scalar nclust			///
					)
{
			tempsum		= J(0,cols(A),.)
			for (i=1; i<=nclust; i++) {						// loop through clusters 1..nclust
				svar	= (pid:==i)							// mimics panelsubmatrix but doesn't require sorted data
				tempsum	= tempsum \ colsum(select(A,svar))
			}
			return(tempsum)

}

// Mata utility for sequential use of solvers
// Default is LU (faster) and if rank-deficient, use QR (drops columns).
// Can also specify lu, qr, svd or chol.
function anysolver (	numeric matrix A,			/// #1
						numeric matrix B,			/// #2
						|							///
						r,							/// #3 scalar, returned with rank of matrix
						string scalar psolver		/// #4
						)
{
	if (args()<4) psolver = ""

	real matrix C

	if (psolver=="") {
		// default behavior: LU, and if rank-deficient, then QR.
		// if nonsquare, need to use quadcross before calling lusolve
		if (rows(A) == cols(A)) {		//  square
			C = lusolve(A, B)			//  default is LU
		}
		else {							//  not square
			C = lusolve(quadcross(A,A),quadcross(A,B))
		}
		if (hasmissing(C)) {			//  not full rank, use QR
			C = qrsolve(A, B, r)		//  r = rank
		}
		else {
			r = cols(A)					//  full rank
		}
	}
	else if (psolver=="lu") {			//  override default, use LU no matter what
		if (rows(A) == cols(A)) {		//  square
			C = lusolve(A, B)
		}
		else {							//  not square
			C = lusolve(quadcross(A,A),quadcross(A,B))
		}
		if (!hasmissing(C)) {			//  if full rank return it, otherwise leave missing
			r = cols(A)
		}
	}
	else if (psolver=="luxx") {			//  override default, use LU+quadcross no matter what
		C = lusolve(quadcross(A,A),quadcross(A,B))
		if (!hasmissing(C)) {			//  if full rank return it, otherwise leave missing
			r = cols(A)
		}
	}
	else if (psolver=="qr") {			//  use QR no matter what
		C = qrsolve(A, B, r)			//  r = rank
	}
	else if (psolver=="qrxx") {			//  use QR+quadcross no matter what
		C = qrsolve(quadcross(A,A),quadcross(A,B),r)
	}
	else if (psolver=="svd") {			//  use SVD no matter what
		C = svsolve(A, B, r)
	}
	else if (psolver=="svdxx") {			//  use SVD+quadcross no matter what
		C = svsolve(quadcross(A,A),quadcross(A,B),r)
	}
	else if (psolver=="chol") {			//  use Cholesky no matter what
		// if nonsquare, need to use quadcross before calling cholsolve
		if (rows(A) == cols(A)) {		//  square
			C = cholsolve(A, B)
		}
		else {							//  not square
			C = cholsolve(quadcross(A,A),quadcross(A,B))
		}
		if (!hasmissing(C)) {			//  if full rank return it, otherwise leave missing
			r = cols(A)
		}
	}

	return(C)
}	// end anysolver


// Mata utility for standardizing; called from Stata.
// Behavior:
// "Standardize" means divide variable x by SD = mean((x-xbar)^2). If model has a constant, also demean x.
// dmflag=1 => treat data as demeaned. mvec set to 0 automatically.
//             means preweighted and demeaned data treated correctly (since weighted mean=0 even though mean is not)
// dmflag=0 => treatment depends on consmodel
//             consmodel=1 => mean is calculated and data are demeaned (hence not suitable for preweighted data)
//             consmodel=0 => just standardize without demeaning (so works with both unweighted and preweighted data)
void s_std(	string scalar Xnames,						//  names of original Stata variables
					string scalar tXnames,				//  names of optional Stata temp vars to initialize/modify
					string scalar touse,				//  full sample
					string scalar toest,				//  sample on which standardization is based
					real scalar consmodel,				//  =1 if constant in model (=> data to be demeaned)
					real scalar dmflag,					//  =1 if data already demeaned
					real scalar transform				//  flag to indicate that data are to be transformed
					)
{

// All varnames should be basic form, no FV or TS operators etc.
// Can include tempvars.

	Xtokens=tokens(Xnames)
	st_view(X, ., Xtokens, touse)		// full sample
	st_view(Xest, ., Xtokens, toest)	// estimation subsample
	
	if (tXnames ~= "") {
		tXflag=1
		tXtokens=tokens(tXnames)
		st_view(tX, ., tXtokens, touse)
	}
	else {
		tXflag=0
	}

	if (dmflag) {
		// treat data as already demeaned
		mvec			= J(1,cols(Xest),0)
		s				= sqrt(mean((Xest):^2))
		// if SD=0 (constant var), convention is to set SD to 1
		if (sum(s:==0)) {
			_editvalue(s,0,1)
		}
		if (transform) {
			if (tXflag) {
				tX[.,.]=X:/s
			}
			else {
				X[.,.]=X:/s
			}
		}
	}
	else {
		// mean and SD based on estimation subsample
		mvec			= mean(Xest)
		s				= sqrt(mean((Xest:-mvec):^2))

		// if SD=0 (constant var), convention is to set SD to 1
		if (sum(s:==0)) {
			_editvalue(s,0,1)
		}
		if (transform) {
			if ((tXflag) & (consmodel)) {
				tX[.,.]=(X:-mvec):/s			//  consmodel => demean temp vars before standardizing
			}
			else if (consmodel) {
				X[.,.]=(X:-mvec):/s				//  consmodel => demean orig vars before standardizing
			}
			else if (tXflag) {
				tX[.,.]=X:/s					//  no consmodel => just standardize temp vars
			}
			else {
				X[.,.]=X:/s						//  no consmodel => just standardize orig vars
			}
		}
	}

	st_matrix("r(stdvec)",s)
	st_matrix("r(mvec)",mvec)
	
}  
//end program s_std

// Mata utility for weighting; called from Stata.
// Weight by sqrt of weighting variable.
void s_wt(			string scalar Xnames,				//  names of original Stata variables
					string scalar tXnames,				//  names of Stata temp vars to initialize/modify
					string scalar touse,				//  full sample
					string scalar wvar,					//  weight variable
					real scalar transform				//  flag to indicate that data are to be transformed
					)
{

// All varnames should be basic form, no FV or TS operators etc.
// Can include tempvars.

	Xtokens=tokens(Xnames)
	st_view(X, ., Xtokens, touse)		// full sample
	
	tXtokens=tokens(tXnames)
	st_view(tX, ., tXtokens, touse)

	st_view(W, ., wvar, touse)
	
	tX[.,.]=X:*sqrt(W)

}  
//end program s_wt


// utility for rlasso
// returns X after centering (if cons present)
// and partialling-out (if any partialled-out vars present)
real matrix centerpartial(	pointer matrix X,		///
							real rowvector Xmvec,	///
							real scalar cons,		///
							pointer matrix Z,		///
							real matrix Zmvec,		///
							real matrix pihat)		//  =0 if nothing to partial out, =projection coefs if Zs.
{
	if ((pihat==0) & (cons==0)) {
		return(*X)
	}
	else if ((pihat==0) & (cons==1)) {
		return((*X) :- Xmvec)
	}
	else {
		return(((*X):-Xmvec) - (((*Z):-Zmvec)*pihat))
	}
}
// end program centerpartial



//********************* MODIFIED FROM OFFICIAL STATA/MATA *********************//
// quadcorrelation modified to allow for no intercept/means
// Based on official quadcorrelation(.):
// version 1.0.1  06jun2006
// version 9.0
// mata:

real matrix m_quadcorr(		real matrix X,			/// #1
							real scalar cons,		/// #2
							|						///
							real colvector w		/// #3
							)
{
		real rowvector  CP
		real rowvector  means
		real matrix	 res
		real scalar	 i, j
		real scalar	 n 

		if (args()<3) w = 1

		CP = quadcross(w,0, X,1)
		n  = cols(CP)
		if (cons) {
			means = CP[|1\n-1|] :/ CP[n]
			res = quadcrossdev(X,0,means, w, X,0,means)
		}
		else {
			res = quadcross(X,0, w, X,0)
		}

		for (i=1; i<=rows(res); i++) {
				res[i,i] = sqrt(res[i,i])
				for (j=1; j<i; j++) {
						res[i,j] = res[j,i] = res[i,j]/(res[i,i]*res[j,j])
				}
		}
		for (i=1; i<=rows(res); i++) res[i,i] = 1
		return(res)
}



//********************* FROM MOREMATA BY BEN JANN *********************//
// Based on:
// mm_quantile.mata
// version 1.0.8  20dec2007  Ben Jann

real matrix mm_quantile(real matrix X, | real colvector w,
 real matrix P, real scalar altdef)
{
	real rowvector result
	real scalar c, cX, cP, r, i

	if (args()<2) w = 1
	if (args()<3) P = (0, .25, .50, .75, 1)'
	if (args()<4) altdef = 0
	if (cols(X)==1 & cols(P)!=1 & rows(P)==1)
	 return(mm_quantile(X, w, P', altdef)')
	if (missing(P) | missing(X) | missing(w)) _error(3351)
	if (rows(w)!=1 & rows(w)!=rows(X)) _error(3200)
	r = rows(P)
	c = max(((cX=cols(X)), (cP=cols(P))))
	if (cX!=1 & cX<c) _error(3200)
	if (cP!=1 & cP<c) _error(3200)
	if (rows(X)==0 | r==0 | c==0) return(J(r,c,.))
	if (c==1) return(_mm_quantile(X, w, P, altdef))
	result = J(r, c, .)
	if (cP==1) for (i=1; i<=c; i++)
	 result[,i] = _mm_quantile(X[,i], w, P, altdef)
	else if (cX==1) for (i=1; i<=c; i++)
	 result[,i] = _mm_quantile(X, w, P[,i], altdef)
	else for (i=1; i<=c; i++)
	 result[,i] = _mm_quantile(X[,i], w, P[,i], altdef)
	return(result)
}

real colvector _mm_quantile(
 real colvector X,
 real colvector w,
 real colvector P,
 real scalar altdef)
{
	real colvector g, j, j1, p
	real scalar N

	if (w!=1) return(_mm_quantilew(X, w, P, altdef))
	N = rows(X)
	p = order(X,1)
	if (altdef) g = P*N + P
	else g = P*N
	j = floor(g)
	if (altdef) g = g - j
	else g = 0.5 :+ 0.5*((g - j):>0)
	j1 = j:+1
	j = j :* (j:>=1)
	_editvalue(j, 0, 1)
	j = j :* (j:<=N)
	_editvalue(j, 0, N)
	j1 = j1 :* (j1:>=1)
	_editvalue(j1, 0, 1)
	j1 = j1 :* (j1:<=N)
	_editvalue(j1, 0, N)
	return((1:-g):*X[p[j]] + g:*X[p[j1]])
}

real colvector _mm_quantilew(
 real colvector X,
 real colvector w,
 real colvector P,
 real scalar altdef)
{
	real colvector Q, pi, pj
	real scalar i, I, j, jj, J, rsum, W
	pointer scalar ww

	I  = rows(X)
	ww = (rows(w)==1 ? &J(I,1,w) : &w)
	if (altdef) return(_mm_quantilewalt(X, *ww, P))
	W  = quadsum(*ww)
	pi = order(X, 1)
	if (anyof(*ww, 0)) {
		pi = select(pi,(*ww)[pi]:!=0)
		I = rows(pi)
	}
	pj = order(P, 1)
	J  = rows(P)
	Q  = J(J, 1, .)
	j  = 1
	jj = pj[1]
	rsum = 0
	for (i=1; i<=I; i++) {
		rsum = rsum + (*ww)[pi[i]]
		if (i<I) {
			if (rsum<P[jj]*W) continue
			if (X[pi[i]]==X[pi[i+1]]) continue
		}
		while (1) {
			if (rsum>P[jj]*W | i==I) Q[jj] = X[pi[i]]
			else Q[jj] = (X[pi[i]] + X[pi[i+1]])/2
			j++
			if (j>J) break
			jj = pj[j]
			if (i<I & rsum<P[jj]*W) break
		}
		if (j>J) break
	}
	return(Q)
}

real colvector _mm_quantilewalt(
 real colvector X,
 real colvector w,
 real colvector P)
{
	real colvector Q, pi, pj
	real scalar i, I, j, jj, J, rsum, rsum0, W, ub, g

	W  = quadsum(w) + 1
	pi = order(X, 1)
	if (anyof(w, 0)) pi = select(pi, w[pi]:!=0)
	I  = rows(pi)
	pj = order(P, 1)
	J  = rows(P)
	Q  = J(J, 1, .)
	rsum = w[pi[1]]
	for (j=1; j<=J; j++) {
		jj = pj[j]
		if (P[jj]*W <= rsum) Q[jj] = X[pi[1]]
		else break
	}
	for (i=2; i<=I; i++) {
		rsum0 = rsum
		rsum = rsum + w[pi[i]]
		if (i<I & rsum < P[jj]*W) continue
		while (1) {
			ub = rsum0+1
			if (P[jj]*W>=ub | X[pi[i]]==X[pi[i-1]]) Q[jj] = X[pi[i]]
			else {
				g = (ub - P[jj]*W) / (ub - rsum0)
				Q[jj] = X[pi[i-1]]*g + X[pi[i]]*(1-g)
			}
			j++
			if (j>J) break
			jj = pj[j]
			if (i<I & rsum < P[jj]*W) break
		}
		if (j>J) break
	}
	return(Q)
}

// END MAIN MATA SECTION
end

*********** CONDITIONAL COMPILATION SECTION *****************
// Section is in Stata environment.

// FTOOLS
// Tell Stata to exit before trying to complile Mata function
// s_fe if required FTOOLS package is not installed.
// Note this runs only once, on loading, and if ftools has not
// been compiled, the check option will trigger compilation.

// ftools can flip matastrict to on, causing program to fail to load.
local mstrict `c(matastrict)'
cap ftools, check
if _rc {
	// fails check, likely not installed, so do not compile
	// _fe Stata program will use (slower) Stata code
	exit
}
else {
	// temporarily set matastrict off
	qui set matastrict off
}

// Compile Mata function s_fe.
version 13
mata:

// FE transformation.
// Uses Sergio Correia's FTOOLS package - faster and does not require the data to be sorted.
void s_fe(		string scalar Xnames,
				string scalar tXnames,
				string scalar Wname,
				string scalar fe,
				string scalar touse,
				string scalar toest)
{

	class Factor scalar F
	F = factor(fe, touse)
	F.panelsetup()

	Xtokens=tokens(Xnames)
	X = st_data( ., Xtokens, touse)
	tXtokens=tokens(tXnames)
	st_view(tX, ., tXtokens, touse)
	
	if (Wname~="") {
		st_view(Wvar, ., Wname, touse)
	}
	else {
		Wvar = J(rows(X),1,1)
	}

	w = F.sort(st_data(., toest, touse))				//  extract toest variable
	counts = panelsum(w:*Wvar, F.info)					//  weighted counts for just the toest subsample
	means = editmissing(panelsum(F.sort(X:*Wvar), w, F.info) :/ counts, 0)
	tX[.,.] = X - means[F.levels, .]

	N_g = F.num_levels
	st_numscalar("r(N_g)",N_g)

}

// End Mata section for s_fe
end

// reset matastrict to what it was prior to loading this file
qui set matastrict `mstrict'

// END CONDITIONAL COMPILATION SECTION

******************** END ALL PROGRAM CODE *******************
exit



**************** ADDITIONAL NOTES ***************************

// The code below is NOT implemented here but is for reference.
// Conditional coding section above works just once, for a single
// possibly uninstalled Mata package.
// To implement for multiple possibly uninstalled Mata packages,
// use the following comment-out trick (hat-tip to Sergio Correia).

// FTOOLS
// Set comment local to "*" before trying to complile Mata function
// s_fe if required FTOOLS package is not installed.
// Note this runs only once, on loading, and if ftools has not
// been compiled, the check option will trigger compilation.

// ftools can flip matastrict to on, causing program to fail to load.
local mstrict `c(matastrict)'
cap ftools, check
if _rc {
	// fails check, likely not installed, so do not complile
	// _fe Stata program will use (slower) Stata code
	loc c *
}
else {
	// temporarily set matastrict off
	qui set matastrict off
}

// Compile Mata function s_fe.
`c' version 13
`c' mata:

// FE transformation.
// Uses Sergio Correia's FTOOLS package - faster and does not require the data to be sorted.
`c' void s_fe(		string scalar Xnames,
`c' 				string scalar tXnames,
`c'					string scalar Wname,
`c' 				string scalar fe,
`c' 				string scalar touse,
`c' 				string scalar toest)
`c' {
`c' 
`c' 	class Factor scalar F
`c' 	F = factor(fe, touse)
`c' 	F.panelsetup()
`c' 
`c' 	Xtokens=tokens(Xnames)
`c' 	X = st_data( ., Xtokens, touse)
`c' 	tXtokens=tokens(tXnames)
`c' 	st_view(tX, ., tXtokens, touse)
`c' 
`c' 	w = F.sort(st_data(., toest, touse))				//  extract toest variable
`c' 	counts = panelsum(w, F.info)
`c' 	means = editmissing(panelsum(F.sort(X), w, F.info) :/ counts, 0)
`c' 	tX[.,.] = X - means[F.levels, .]
`c' 
`c' 	N_g = F.num_levels
`c' 	st_numscalar("r(N_g)",N_g)
`c' 
`c' }
`c' 
`c' // End Mata section for s_fe
`c' end

// reset matastrict
qui set matastrict `mstrict'

// Reset comment local for next package.
loc c

// ... further possibly uninstalled Mata packages ...

exit

