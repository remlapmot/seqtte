cscript seqtte adofile seqtte

* ------------------------------------------------------------
* Generate a synthetic person-period dataset for testing
* ------------------------------------------------------------
* 200 individuals, 10 periods (0–9)
* Treatment follows a Markov chain so both A_lag==0 and A_lag==1
* strata have variation in treatment — required for PP weight models
* Covariate: age_grp (binary)

clear
set seed 42
set obs 200

gen id      = _n
gen age_grp   = mod(id, 2)
gen age_grp_0 = age_grp

expand 10
bysort id: gen time = _n - 1

* Markov treatment:
*   P(A=1 | A_lag=0) = 0.25   (initiation)
*   P(A=1 | A_lag=1) = 0.70   (continuation / some stopping)
gen treatment = .
bysort id (time): replace treatment = (runiform() < 0.25) if time == 0
bysort id (time): replace treatment = ///
    cond(treatment[_n-1] == 0, (runiform() < 0.25), (runiform() < 0.70)) ///
    if time > 0

* Period-specific outcome: 1 only in the period the event occurs
gen double u = runiform()
gen outcome = (u < 0.07 - 0.03 * treatment)

* Truncate at first event
bysort id (time): gen cumev = sum(outcome)
drop if cumev > 1
replace outcome = 0 if cumev == 1 & !outcome

drop cumev u
sort id time

* ------------------------------------------------------------
* Test 1: ITT, no covariates
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment)

assert e(N)       > 0
assert e(N_indiv) == 200
assert e(N_orig)  > 0
assert e(N_exp)  >= e(N_orig)
assert `"`e(cmd)'"'       == "seqtte"
assert `"`e(estimator)'"' == "itt"

* ------------------------------------------------------------
* Test 2: ITT, with covariate
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    covariates(age_grp)

assert e(N) > 0
assert `"`e(estimator)'"' == "itt"

* ------------------------------------------------------------
* Test 3: ITT, if/in restriction
* ------------------------------------------------------------
seqtte outcome if time <= 7, id(id) time(time) treatment(treatment)

assert e(N) > 0

* ------------------------------------------------------------
* Test 4: PP, unstabilized weights
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    covariates(age_grp) ///
    estimator(pp) wdenominator(age_grp)

assert e(N) > 0
assert `"`e(estimator)'"' == "pp"

* ------------------------------------------------------------
* Test 5: PP, stabilized weights (wnumerator = subset of wdenominator)
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    covariates(age_grp) ///
    estimator(pp) wdenominator(age_grp) wnumerator(age_grp_0)

assert e(N) > 0
assert `"`e(estimator)'"' == "pp"

* ------------------------------------------------------------
* Test 6: PP, custom truncation threshold
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    estimator(pp) wdenominator(age_grp) truncation(10)

assert e(N) > 0

* ------------------------------------------------------------
* Test 7: error if pp specified without wdenominator
* ------------------------------------------------------------
rcof `"seqtte outcome, id(id) time(time) treatment(treatment) estimator(pp)"' == 198

di as txt "seqtte cscript passed"
