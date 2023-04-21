mata:
mata clear

void betaEstimatingEquation(todo, beta, r, pihat, xtemp, fullyObservedCovariates, outcome, pivars, lambda, v, g, H)
{
	n = st_nobs()

	epsilon = outcome - (xtemp, fullyObservedCovariates, J(n,1,1)) * beta'
	firstPart = ((epsilon :* r)) :* (xtemp, fullyObservedCovariates, J(n,1,1))
	
	numImputations = strtoreal(st_local("m"))
	st_view(ximpute, ., range(st_varindex("ximpute1"), st_varindex("ximpute1")+(numImputations-1), 1)')
	phi = - (outcome - (ximpute[.,1], fullyObservedCovariates, J(n,1,1)) * beta') :* (ximpute[.,1], fullyObservedCovariates, J(n,1,1))
	
	for (i=2; i<=numImputations; i++) {
		phi = phi :+ (- (outcome - (ximpute[.,i], fullyObservedCovariates, J(n,1,1)) * beta') :* (ximpute[.,i], fullyObservedCovariates, J(n,1,1)))
	}
	phi = phi :/ numImputations
			
	phialpha = (pivars, J(n,1,1))
	
	G_alpha = - (((pihat :* (1 :- pihat)) :* phialpha)' * phi ) :/ n
		
	//calculate tildephi
	tildephi = phi + phialpha * (luinv((((pihat :* (1 :- pihat)) :* phialpha)' * phialpha ) :/ n) * G_alpha)
	newSecondPart = (r - pihat) :* tildephi	
		
	estimatingFunction = firstPart + (r - pihat) :* (tildephi * lambda')
	v = -sum(colsum(estimatingFunction):^2)	
	
	// calculate derivative
	derivpart1 = -(r :* (xtemp, fullyObservedCovariates, J(n,1,1)))' * (xtemp, fullyObservedCovariates, J(n,1,1))
	deriv = derivpart1
	
	if (todo==1) {
		g = -2 :* colsum(estimatingFunction) * deriv
	}
}

real matrix sandwichVar(beta, r, pihat, xtemp, fullyObservedCovariates, outcome, pivars, lambda)
{
	n = st_nobs()
	
	G_theta = ( -(r :* (xtemp, fullyObservedCovariates, J(n,1,1)))' * (xtemp, fullyObservedCovariates, J(n,1,1)) ) :/ n
	
	epsilon = outcome - (xtemp, fullyObservedCovariates, J(n,1,1)) * beta'
	firstPart = ((epsilon :* r)) :* (xtemp, fullyObservedCovariates, J(n,1,1))
	
	numImputations = strtoreal(st_local("m"))
	st_view(ximpute, ., range(st_varindex("ximpute1"), st_varindex("ximpute1")+(numImputations-1), 1)')
	phi = - (outcome - (ximpute[.,1], fullyObservedCovariates, J(n,1,1)) * beta') :* (ximpute[.,1], fullyObservedCovariates, J(n,1,1))
	
	for (i=2; i<=numImputations; i++) {
		phi = phi :+ (- (outcome - (ximpute[.,i], fullyObservedCovariates, J(n,1,1)) * beta') :* (ximpute[.,i], fullyObservedCovariates, J(n,1,1)))
	}
	phi = phi :/ numImputations
	
	phialpha = (pivars, J(n,1,1))
	
	G_alpha = - (((pihat :* (1 :- pihat)) :* phialpha)' * phi ) :/ n
		
	//calculate tildephi
	tildephi = phi + phialpha * (luinv((((pihat :* (1 :- pihat)) :* phialpha)' * phialpha ) :/ n) * G_alpha)
	newSecondPart = (r - pihat) :* tildephi
		
	estimatingFunction = firstPart + (r - pihat) :* (tildephi * lambda')
	
	inner = (estimatingFunction' * estimatingFunction) :/ n
	
	sigma = (luinv(G_theta) * inner * luinv(G_theta)') :/ n
		
	return(sigma)
}

real matrix estimateLambda(beta, r, pihat, xtemp, fullyObservedCovariates, outcome, pivars)
{
	n = st_nobs()
	
	G_theta = ( -(r :* (xtemp, fullyObservedCovariates, J(n,1,1)))' * (xtemp, fullyObservedCovariates, J(n,1,1)) ) :/ n
	
	epsilon = outcome - (xtemp, fullyObservedCovariates, J(n,1,1)) * beta'
	firstPart = ((epsilon :* r)) :* (xtemp, fullyObservedCovariates, J(n,1,1))
	
	numImputations = strtoreal(st_local("m"))
	st_view(ximpute, ., range(st_varindex("ximpute1"), st_varindex("ximpute1")+(numImputations-1), 1)')
	phi = - (outcome - (ximpute[.,1], fullyObservedCovariates, J(n,1,1)) * beta') :* (ximpute[.,1], fullyObservedCovariates, J(n,1,1))
	
	for (i=2; i<=numImputations; i++) {
		phi = phi :+ (- (outcome - (ximpute[.,i], fullyObservedCovariates, J(n,1,1)) * beta') :* (ximpute[.,i], fullyObservedCovariates, J(n,1,1)))
	}
	phi = phi :/ numImputations
	
	phialpha = (pivars, J(n,1,1))
	
	G_alpha = - (((pihat :* (1 :- pihat)) :* phialpha)' * phi ) :/ n
		
	//calculate tildephi
	tildephi = phi + phialpha * (luinv((((pihat :* (1 :- pihat)) :* phialpha)' * phialpha ) :/ n) * G_alpha)
	newSecondPart = (r - pihat) :* tildephi
		
	//calculate lambda
	return(- (firstPart'*newSecondPart :/ n) * luinv(newSecondPart'*newSecondPart :/ n))
}
	

void estimatingBeta(string scalar r_name, string scalar pihat_name, string scalar x_name, string scalar fullyObservedCovariates_name, string scalar outcome_name, string scalar pivars_name)
{
	st_view(r, ., r_name)
	st_view(pihat, ., pihat_name)
	st_view(x, ., x_name)
	st_view(fullyObservedCovariates, ., fullyObservedCovariates_name)
	st_view(outcome, ., outcome_name)
	st_view(pivars, ., pivars_name)
		
	S = optimize_init()
	optimize_init_evaluator(S, &betaEstimatingEquation())
	optimize_init_evaluatortype(S, "d1")
	optimize_init_params(S, st_matrix("ccestimates"))
	
	optimize_init_argument(S, 1, r)
	optimize_init_argument(S, 2, pihat)
	optimize_init_argument(S, 3, x)
	optimize_init_argument(S, 4, fullyObservedCovariates)
	optimize_init_argument(S, 5, outcome)
	optimize_init_argument(S, 6, pivars)
	//optimize_init_tracelevel(S, "gradient")
	//optimize_init_conv_maxiter(S, 1)
	
	//calculate lambda
	lambda = estimateLambda(st_matrix("ccestimates"), r, pihat, x, fullyObservedCovariates, outcome, pivars)
	optimize_init_argument(S, 7, lambda)
	
	betahat = optimize(S)
	st_matrix("aipwestimates", betahat)
	
	st_matrix("aipwvar", sandwichVar(betahat, r, pihat, x, fullyObservedCovariates, outcome, pivars,lambda))
}

end


capture program drop augcca
program define augcca,eclass
version 11
syntax [varlist(default=none)], y(varname) x(varname) imps(varlist) pivars(varlist)
*generate missingness indicator
gen my_r=(`x'!=.)

di as result "Missingness model"
logistic my_r `pivars'
predict pihat, pr

*fit model to complete cases
di as result "Complete case analysis"
reg `y' `x' `varlist'
matrix ccestimates=e(b)

quietly gen xtemp = `x'
quietly replace xtemp = 0 if xtemp==.

local i=0
quietly foreach var of varlist `imps' {
	local i = `i' + 1
	gen ximpute`i' = `var'
}
local m = `i'

*estimate beta using augmented complete case analysis
di as result "Augmented complete case analysis"

mata: estimatingBeta("my_r", "pihat", "xtemp", "`varlist'", "`y'", "`pivars'")

return clear
ereturn clear
drop my_r pihat xtemp

quietly forvalues i=1(1)`m' {
	drop ximpute`i'
}

local expvarnames : colnames ccestimates
local outvarname: rownames ccestimates

matrix colnames aipwestimates = `expvarnames'
matrix rownames aipwestimates = `outvarname'
matrix colnames aipwvar = `expvarnames'
matrix rownames aipwvar = `expvarnames'

ereturn post aipwestimates aipwvar


ereturn display

end

