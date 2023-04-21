{smcl}
{cmd:help augcca}
{hline}

{title:Title}

{phang}
{bf:augcca} Augmented complete case analysis for partially observed covariates (augcca)

{title:Syntax}

{p 8 17 2} {cmd:augcca} [varlist], y(varname) x(varname) imps(varlist) pivars(varlist)

{synoptset 40 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}

{synopt :* {opt y(varname)}}the outcome variable in the linear regression model{p_end}
{synopt :* {opt x(varname)}}the partially observed covariate in the linear regression model{p_end}
{synopt :* {opt imps(varlist)}}the variables containing multiple imputations of x{p_end}
{synopt :* {opt pivars(varlist)}}the variables to be included in the logistic regression missingness model for x{p_end}

{synoptline}
{p2coldent :* denotes required option}{p_end}

{pstd}
[varlist] specifies the fully observed variables z in the linear regression outcome model.

{title:Description}

{pstd}
{cmd:augcca} estimates the parameters of a linear regression model for outcome y on covariates x and z, where x contains missing values,
using augmented complete case analysis. Estimation is performed assuming that missingness in x is independent of y, conditional on x and z. 
This is the same assumption required for validity of complete case analysis (CCA).

{title:Remarks}
{marker remarks}{...}

{pstd}
{cmd:augcca} gains efficiency over CCA by using information on y and z in the incomplete cases. It assumes that missingness in x,
conditional on y and z, follows a logistic regression model with covariates given by the pivars option. Correct specification of this
model is required for estimates to be consistent.

{pstd}
The command makes use of Monte-Carlo integration to approximate the expectation of the full data estimating function given y and z.
To do this the user must supply multiple imputations of x {bf:for all subjects} using the imp option. The command assumes
these imputations are stored in so called wide format. Using a larger number of imputations will result in somewhat more efficient 
inferences. The consistency of estimates does not rely on the imputation model used being correctly specified. Also note that the 
imputations need not be proper.

{pstd}
It is important to note that although the method fits a model for missingness in x given y and z, estimation is performed assuming
that missingness in x is independent of y, given z {bf: and x}.

{pstd}
Further details regarding the algorithm used by {bf:augcca} can be found in the article referenced below.

{title:Example}

{pstd}Assuming we have already created 10 imputations of x, stored in variables imp1-imp10, and that conditional on y and z, 
missingness in x follows a logistic regression with y and z as linear covariates, we can call {bf: augcca} by

{phang2}{cmd:augcca z, y(y) x(x) imps(imp1-imp10) pivars(y z)}{p_end}

{title:References}

{phang}Jonathan W. Bartlett, James R. Carpenter, Kate Tilling, Stijn Vansteelandt. {browse "http://dx.doi.org/10.1093/biostatistics/kxu023":Improving upon the efficiency of complete case analysis when covariates are MNAR.} Biostatistics (2014), 15; 719-730{p_end}

	
{title:Author}

    Jonathan Bartlett, London School of Hygiene & Tropical Medicine, UK
    jonathan.bartlett1@lshtm.ac.uk
	{browse "https://thestatsgeek.com/"}
