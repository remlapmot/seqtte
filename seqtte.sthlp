{smcl}
{* *! version 0.1.0  29apr2026  Tom Palmer}{...}
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
[{cmd:covariates(}{it:varlist}{cmd:)}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt id(varname)}}individual identifier variable{p_end}
{synopt:{opt time(varname)}}integer calendar time variable{p_end}
{synopt:{opt treatment(varname)}}binary treatment indicator (0/1){p_end}
{synoptline}
{syntab:Optional}
{synopt:{opt covariates(varlist)}}adjustment covariates; values at trial entry are used{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:seqtte} estimates the intent-to-treat (ITT) effect of a sustained
treatment strategy using sequential target trial emulation
({help seqtte##HR2016:Hernán and Robins, 2016}).

{pstd}
The input data should be in long (person-period) format.
An individual is eligible to enter a trial at time {it:t} if they
have not received treatment in any earlier period.
For each eligible person-period, a trial is initiated; the person
is then followed from that trial entry time until the end of their
observed follow-up.

{pstd}
The ITT estimand compares outcomes between those who were assigned
(i.e., observed to initiate) treatment at trial entry versus those
who were not, regardless of subsequent treatment changes.

{pstd}
Internally, {cmd:seqtte} expands the dataset to create the pooled
person-trial observations, defines the binary outcome within each
trial, and fits a pooled logistic regression model with quadratic
polynomial terms for follow-up time within trial and trial number.
Standard errors are clustered by individual to account for the fact
that each person can contribute to multiple trials.

{pstd}
The {it:outcomevar} should equal 1 in the period the event first
occurs and 0 otherwise; the event period should be the last record
for each individual in the input data.

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
{opt covariates(varlist)} specifies adjustment covariates.
In the expanded dataset each person-trial record takes the covariate
values from the trial entry period, so both time-fixed and
time-varying covariates are included at their trial-entry values.

{marker examples}{...}
{title:Examples}

{pstd}
Generate a simple synthetic dataset and fit the ITT model:

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. set obs 200}{p_end}
{phang2}{cmd:. gen id = _n}{p_end}
{phang2}{cmd:. expand 10}{p_end}
{phang2}{cmd:. bysort id: gen time = _n - 1}{p_end}
{phang2}{cmd:. gen treatment = (id <= 100) * (time == 0)}{p_end}
{phang2}{cmd:. gen treatment_ever = treatment}{p_end}
{phang2}{cmd:. bysort id (time): replace treatment_ever = max(treatment_ever, treatment_ever[_n-1])}{p_end}
{phang2}{cmd:. gen treatment2 = treatment_ever}{p_end}
{phang2}{cmd:. gen u = runiform()}{p_end}
{phang2}{cmd:. gen outcome = (u < 0.05)}{p_end}
{phang2}{cmd:. bysort id (time): gen cumev = sum(outcome)}{p_end}
{phang2}{cmd:. drop if cumev > 1}{p_end}
{phang2}{cmd:. seqtte outcome, id(id) time(time) treatment(treatment2)}{p_end}

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
{synopt:{cmd:e(estimator)}}{cmd:itt}{p_end}
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
