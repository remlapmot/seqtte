{smcl}
{* *! version 0.2.0  13may2026  Tom Palmer}{...}
{vieweralsosee "seqtte" "help seqtte"}{...}
{viewerjumpto "Syntax" "seqtte##syntax"}{...}
{viewerjumpto "Description" "seqtte##description"}{...}
{viewerjumpto "Options" "seqtte##options"}{...}
{viewerjumpto "Examples" "seqtte##examples"}{...}
{viewerjumpto "Stored results" "seqtte##results"}{...}
{viewerjumpto "References" "seqtte##references"}{...}
{viewerjumpto "Author" "seqtte##author"}{...}
{title:Title}

{phang}
{bf:seqtte} {hline 2} Sequential target trial emulation

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:seqtte} {it:outcomevar} {ifin}{cmd:,}
{cmdab:id(}{it:varname}{cmd:)}
{cmdab:time(}{it:varname}{cmd:)}
{cmdab:treatment(}{it:varname}{cmd:)}
[{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt id(varname)}}individual identifier variable{p_end}
{synopt:{opt time(varname)}}integer calendar time variable{p_end}
{synopt:{opt treatment(varname)}}binary treatment indicator (0/1){p_end}
{synoptline}
{syntab:Optional}
{synopt:{opt covariates(varlist)}}adjustment covariates for the outcome model{p_end}
{synopt:{opt estimator(string)}}{cmd:itt} (default) or {cmd:pp}{p_end}
{synopt:{opt wdenominator(varlist)}}denominator weight model covariates; required for {cmd:pp}{p_end}
{synopt:{opt wnumerator(varlist)}}numerator weight model covariates; if supplied, stabilized weights are used{p_end}
{synopt:{opt truncation(#)}}truncation threshold for cumulative weights; default 25{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:seqtte} estimates the causal effect of a sustained treatment strategy
using sequential target trial emulation
({help seqtte##HR2016:Hernán and Robins, 2016}).
Two estimators are available.

{pstd}
{ul:Intent-to-treat (ITT)} ({cmd:estimator(itt)}, the default).
Estimates the effect of being assigned to treatment at trial entry,
regardless of subsequent treatment changes.
No weighting is used.

{pstd}
{ul:Per-protocol (PP)} ({cmd:estimator(pp)}).
Estimates the effect of sustained adherence to the assigned treatment strategy.
Individuals are censored at the period in which they deviate from their
assigned treatment, and inverse probability of censoring weights (IPCW)
are used to adjust for informative censoring.

{pstd}
{ul:Input data.}
The data should be in long (person-period) format with one row per individual
per time period.
{it:outcomevar} should equal 1 only in the period the event first occurs and
0 otherwise.
The {it:time} variable should take consecutive integer values.

{pstd}
{ul:Algorithm.}
For each eligible person-period (i.e. periods in which the individual has
not yet received treatment), a trial is initiated.
The individual is then followed from that trial entry time to the end of
their observed follow-up.
Data are expanded so that each person can contribute to multiple trials.
A pooled logistic regression model is then fitted with quadratic polynomial
terms for follow-up time within trial and trial number, with standard errors
clustered by individual.

{pstd}
{ul:Weight models (PP only).}
Four logistic regression models are fitted on the pre-expansion data,
stratified by prior treatment status ({it:A_lag} = 0 or 1):
a denominator model including {cmd:wdenominator()} covariates
and a cubic polynomial in calendar time,
and (if {cmd:wnumerator()} is supplied) a numerator model including
{cmd:wnumerator()} covariates and the same time polynomial.
Unstabilized weights are used when {cmd:wnumerator()} is omitted;
stabilized weights (numerator/denominator) are used when it is supplied.
Cumulative products of per-period weights are formed within each trial,
then truncated at {cmd:truncation()}.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt id(varname)} specifies the variable identifying individuals.

{phang}
{opt time(varname)} specifies the integer calendar time variable.
Consecutive integer values represent consecutive periods.

{phang}
{opt treatment(varname)} specifies the binary treatment indicator
(0 = untreated, 1 = treated).

{dlgtab:Optional}

{phang}
{opt covariates(varlist)} specifies adjustment covariates for the
outcome model.
In the expanded dataset each person-trial record takes covariate
values from the trial entry period, so both time-fixed and
time-varying covariates are included at their trial-entry values.

{phang}
{opt estimator(string)} specifies the estimator: {cmd:itt} (default)
for the intent-to-treat effect, or {cmd:pp} for the per-protocol effect.

{phang}
{opt wdenominator(varlist)} specifies the covariates for the denominator
weight models.
These are fitted on the pre-expansion data and should include all
time-varying confounders of the treatment–outcome relationship.
Required when {cmd:estimator(pp)} is specified.

{phang}
{opt wnumerator(varlist)} specifies the covariates for the numerator
(stabilization) weight models.
These should be a subset of {cmd:wdenominator()}, typically restricted to
baseline (study-entry) values of covariates.
When omitted, unstabilized weights are used.

{phang}
{opt truncation(#)} specifies the upper truncation threshold applied
to the cumulative IPW weights.
Default is 25.

{marker examples}{...}
{title:Examples}

{pstd}Setup: generate a synthetic person-period dataset{p_end}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. set obs 500}{p_end}
{phang2}{cmd:. gen id = _n}{p_end}
{phang2}{cmd:. gen age = rnormal(50, 10)}{p_end}
{phang2}{cmd:. expand 15}{p_end}
{phang2}{cmd:. bysort id: gen time = _n - 1}{p_end}
{phang2}{cmd:. gen treatment = 0}{p_end}
{phang2}{cmd:. bysort id (time): replace treatment = 1 if time > 4 & id <= 250}{p_end}
{phang2}{cmd:. gen outcome = (runiform() < 0.05)}{p_end}
{phang2}{cmd:. bysort id (time): gen cumev = sum(outcome)}{p_end}
{phang2}{cmd:. drop if cumev > 1}{p_end}

{pstd}ITT estimator{p_end}

{phang2}{cmd:. seqtte outcome, id(id) time(time) treatment(treatment) covariates(age)}{p_end}

{pstd}PP estimator with unstabilized weights{p_end}

{phang2}{cmd:. seqtte outcome, id(id) time(time) treatment(treatment) covariates(age) estimator(pp) wdenominator(age)}{p_end}

{pstd}PP estimator with stabilized weights{p_end}

{phang2}{cmd:. seqtte outcome, id(id) time(time) treatment(treatment) covariates(age) estimator(pp) wdenominator(age) wnumerator(age)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:seqtte} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations in the pooled logistic regression{p_end}
{synopt:{cmd:e(N_indiv)}}number of individuals in the original data{p_end}
{synopt:{cmd:e(N_orig)}}number of observations in the original data{p_end}
{synopt:{cmd:e(N_exp)}}number of observations in the expanded dataset{p_end}
{synopt:{cmd:e(r2_p)}}pseudo R-squared{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:seqtte}{p_end}
{synopt:{cmd:e(estimator)}}{cmd:itt} or {cmd:pp}{p_end}
{synopt:{cmd:e(depvar)}}name of the outcome variable{p_end}
{synopt:{cmd:e(clustvar)}}name of the cluster variable{p_end}
{synopt:{cmd:e(vcetype)}}{cmd:Robust}{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector (log-odds scale){p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}

{marker references}{...}
{title:References}

{marker HR2016}{...}
{phang}
Hernán MA, Robins JM. 2016.
Using Big Data to Emulate a Target Trial When a Randomized Trial Is Not Available.
{it:American Journal of Epidemiology} 183(8): 758–764.

{phang}
Danaei G, Rodríguez LAG, Cantero OF, Logan R, Hernán MA. 2013.
Observational data for comparative effectiveness research:
an emulation of randomised trials of statins and primary prevention of coronary heart disease.
{it:Statistical Methods in Medical Research} 22(1): 70–96.

{phang}
Maringe C, Benitez Majano S, Exarchakou A, et al. 2020.
Reflections on modern methods: trial emulation in the presence of immortal-time bias.
Assessing the benefit of major surgery for elderly lung cancer patients using
observational data.
{it:International Journal of Epidemiology} 49(5): 1719–1729.

{marker author}{...}
{title:Author}

{phang}
Tom Palmer, University of Bristol, Bristol, UK.
{browse "mailto:remlapmot@hotmail.com":remlapmot@hotmail.com}

{phang}
Michalis Katsoulis, UCL, London, UK.

{phang}
Please report any bugs or feature requests at
{browse "https://github.com/remlapmot/seqtte/issues"}.
{p_end}
