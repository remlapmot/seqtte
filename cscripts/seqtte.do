cscript seqtte adofile seqtte

* ------------------------------------------------------------
* Generate a synthetic person-period dataset for testing
* ------------------------------------------------------------
* 200 individuals, up to 10 periods
* Treatment can be initiated at any period (time-varying)
* Covariate: age_grp (binary)

clear
set seed 42
set obs 200

gen id      = _n
gen age_grp = mod(id, 2)       // binary baseline covariate
gen age_grp_0 = age_grp        // study-baseline version (same here)

* Individuals initiate treatment at a random period (1–8); some never treated
gen trt_start = ceil(runiform() * 8) if runiform() > 0.3
replace trt_start = . if runiform() > 0.8  // ~20% never treated

expand 10
bysort id: gen time = _n - 1

* Treatment: initiated at trt_start and sustained thereafter
gen treatment = (time >= trt_start) if !missing(trt_start)
replace treatment = 0 if missing(trt_start)

* Period-specific outcome (lower hazard for treated)
gen double u = runiform()
gen outcome = (u < 0.06 - 0.03 * treatment)

* Keep only up to the first event (person exits after event)
bysort id (time): gen cumev = sum(outcome)
drop if cumev > 1
replace outcome = 0 if cumev == 1 & u >= (0.06 - 0.03 * treatment)

drop trt_start cumev u
sort id time

* ------------------------------------------------------------
* Test 1: ITT, no covariates
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment)

assert e(N)      > 0
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
* Test 5: PP, stabilized weights
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
rcof "seqtte outcome, id(id) time(time) treatment(treatment) estimator(pp)" == 198

di as txt "seqtte cscript passed"
