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
* Test 7: unweighted PP — censoring applied, no IPCW weights
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    estimator(pp)

assert e(N) > 0
assert `"`e(estimator)'"' == "pp"
* Post-censoring N must be <= expanded N (some controls are censored)
assert e(N) <= e(N_exp)

* ------------------------------------------------------------
* Test 8: factor-variable notation in covariates and weight models
* ------------------------------------------------------------
gen grp = mod(id, 3)

seqtte outcome, id(id) time(time) treatment(treatment) ///
    covariates(i.age_grp i.grp)

assert e(N) > 0
assert `"`e(estimator)'"' == "itt"

seqtte outcome, id(id) time(time) treatment(treatment) ///
    covariates(i.age_grp i.grp) ///
    estimator(pp) ///
    wdenominator(i.age_grp i.grp) wnumerator(i.age_grp)

assert e(N) > 0
assert `"`e(estimator)'"' == "pp"

* ------------------------------------------------------------
* Test 9: selection_random, ITT — runs and reduces dataset
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    selectionrandom selectionsample(0.5) seed(42)

assert e(N) > 0
assert e(N_sel) > 0
assert e(N_sel) <= e(N_exp)
assert `"`e(estimator)'"' == "itt"

* ------------------------------------------------------------
* Test 10: selection_random reproducibility — same seed → same N
* ------------------------------------------------------------
local n_rep1 = e(N)

seqtte outcome, id(id) time(time) treatment(treatment) ///
    selectionrandom selectionsample(0.5) seed(42)

assert e(N) == `n_rep1'

* ------------------------------------------------------------
* Test 11: selection_random with PP estimator
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    covariates(age_grp) estimator(pp) wdenominator(age_grp) ///
    selectionrandom selectionsample(0.3) seed(7)

assert e(N) > 0
assert `"`e(estimator)'"' == "pp"
assert e(N_sel) <= e(N_exp)
assert e(selection_sample) == .3

* ------------------------------------------------------------
* Test 12: selection_random invalid selectionsample
* ------------------------------------------------------------
rcof `"seqtte outcome, id(id) time(time) treatment(treatment) selectionrandom selectionsample(1.5)"' == 198
rcof `"seqtte outcome, id(id) time(time) treatment(treatment) selectionrandom selectionsample(0)"' == 198

* ------------------------------------------------------------
* Test 13: bootstrap, ITT — scalars and matrix returned
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    bootstrap(50) seed(42)

assert e(N) > 0
assert e(N_boot) > 0
assert e(N_boot) <= 50
assert !missing(e(bs_se))
assert !missing(e(bs_ll))
assert !missing(e(bs_ul))
assert e(bs_ll) < e(bs_ul)
assert rowsof(e(bs_b)) == 50
assert `"`e(estimator)'"' == "itt"

* ------------------------------------------------------------
* Test 14: bootstrap reproducibility — same seed → same SE
* ------------------------------------------------------------
local se_rep1 = e(bs_se)

seqtte outcome, id(id) time(time) treatment(treatment) ///
    bootstrap(50) seed(42)

assert reldif(e(bs_se), `se_rep1') < 1e-10

* ------------------------------------------------------------
* Test 15: bootstrap with PP estimator
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    covariates(age_grp) estimator(pp) wdenominator(age_grp) ///
    bootstrap(30) seed(7)

assert e(N) > 0
assert `"`e(estimator)'"' == "pp"
assert e(N_boot) > 0
assert !missing(e(bs_se))

* ------------------------------------------------------------
* Test 16: bootstrap combined with selection_random
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    selectionrandom selectionsample(0.5) bootstrap(30) seed(99)

assert e(N) > 0
assert e(N_sel) <= e(N_exp)
assert e(N_boot) > 0
assert !missing(e(bs_se))

* ------------------------------------------------------------
* Test 17: bootstrap when selectionrandom drops whole clusters
*   Regression test for r(498). n_indiv is counted before
*   selectionrandom, so when selection removes every (id, trial) pair
*   for some individuals the surviving cluster count falls below
*   n_indiv. bsample must resample the post-selection cluster count,
*   not n_indiv, or every replicate errors and 0 succeed.
* ------------------------------------------------------------
preserve

clear
set seed 12345
set obs 300
gen id = _n
* Odd ids: a single period, never treated (all control-arm), so
* selectionsample(0.5) drops roughly half of them entirely.
gen byte shortfu = mod(id, 2)
expand cond(shortfu, 1, 8)
bysort id: gen time = _n - 1
gen treatment = .
bysort id (time): replace treatment = 0 if time == 0
bysort id (time): replace treatment = ///
    cond(shortfu, 0, cond(treatment[_n-1] == 0, (runiform() < 0.3), 1)) ///
    if time > 0
gen double u = runiform()
gen outcome = (u < 0.1 - 0.04 * treatment)
bysort id (time): gen cumev = sum(outcome)
drop if cumev > 1
replace outcome = 0 if cumev == 1 & !outcome
drop cumev u

seqtte outcome, id(id) time(time) treatment(treatment) ///
    selectionrandom selectionsample(0.5) bootstrap(30) seed(99)

assert e(N) > 0
assert e(N_sel) < e(N_exp)
assert e(N_boot) == 30
assert !missing(e(bs_se))
assert !missing(e(bs_ll))
assert !missing(e(bs_ul))

restore

* ------------------------------------------------------------
* Test 18: bootstrap invalid (negative)
* ------------------------------------------------------------
rcof `"seqtte outcome, id(id) time(time) treatment(treatment) bootstrap(-1)"' == 198

* ------------------------------------------------------------
* Test 19: unweighted PP with covariates
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    covariates(age_grp) estimator(pp)

assert e(N) > 0
assert `"`e(estimator)'"' == "pp"

* ------------------------------------------------------------
* Test 20: unweighted PP with selection_random
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    covariates(age_grp) estimator(pp) ///
    selectionrandom selectionsample(0.5) seed(42)

assert e(N) > 0
assert `"`e(estimator)'"' == "pp"
assert e(N_sel) <= e(N_exp)

* ------------------------------------------------------------
* Test 21: unweighted PP with bootstrap
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    covariates(age_grp) estimator(pp) ///
    bootstrap(30) seed(7)

assert e(N) > 0
assert `"`e(estimator)'"' == "pp"
assert e(N_boot) > 0
assert !missing(e(bs_se))
assert !missing(e(bs_ll))
assert !missing(e(bs_ul))
assert e(bs_ll) < e(bs_ul)

* ------------------------------------------------------------
* Test 22: follow-up counts — ITT
*   - all four scalars returned, positive, non-missing
*   - nonunique totals partition the expanded dataset
*   - unique <= nonunique (can't have more individuals than intervals)
*   - unique <= N_indiv (can't exceed original sample size)
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment)

assert !missing(e(N_nonuniq_arm0)) & e(N_nonuniq_arm0) > 0
assert !missing(e(N_nonuniq_arm1)) & e(N_nonuniq_arm1) > 0
assert !missing(e(N_uniq_arm0))    & e(N_uniq_arm0)    > 0
assert !missing(e(N_uniq_arm1))    & e(N_uniq_arm1)    > 0

* intervals sum to the expanded N (= regression N for ITT)
assert e(N_nonuniq_arm0) + e(N_nonuniq_arm1) == e(N_exp)
assert e(N_nonuniq_arm0) + e(N_nonuniq_arm1) == e(N)

* unique individuals bounded by intervals and total sample
assert e(N_uniq_arm0) <= e(N_nonuniq_arm0)
assert e(N_uniq_arm1) <= e(N_nonuniq_arm1)
assert e(N_uniq_arm0) <= e(N_indiv)
assert e(N_uniq_arm1) <= e(N_indiv)

* ------------------------------------------------------------
* Test 23: follow-up counts — unweighted PP
*   nonunique totals must equal the post-censoring regression N
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    estimator(pp)

assert !missing(e(N_nonuniq_arm0)) & e(N_nonuniq_arm0) > 0
assert !missing(e(N_nonuniq_arm1)) & e(N_nonuniq_arm1) > 0
assert !missing(e(N_uniq_arm0))    & e(N_uniq_arm0)    > 0
assert !missing(e(N_uniq_arm1))    & e(N_uniq_arm1)    > 0

* intervals sum to post-censoring regression N (< N_exp due to censoring)
assert e(N_nonuniq_arm0) + e(N_nonuniq_arm1) == e(N)
assert e(N_nonuniq_arm0) + e(N_nonuniq_arm1) < e(N_exp)

assert e(N_uniq_arm0) <= e(N_nonuniq_arm0)
assert e(N_uniq_arm1) <= e(N_nonuniq_arm1)

* ------------------------------------------------------------
* Test 24: follow-up counts — ITT with selectionrandom
*   nonunique totals must equal the post-selection regression N
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    selectionrandom selectionsample(0.5) seed(42)

assert !missing(e(N_nonuniq_arm0)) & e(N_nonuniq_arm0) > 0
assert !missing(e(N_nonuniq_arm1)) & e(N_nonuniq_arm1) > 0
assert !missing(e(N_uniq_arm0))    & e(N_uniq_arm0)    > 0
assert !missing(e(N_uniq_arm1))    & e(N_uniq_arm1)    > 0

assert e(N_nonuniq_arm0) + e(N_nonuniq_arm1) == e(N)
assert e(N_nonuniq_arm0) + e(N_nonuniq_arm1) == e(N_sel)

assert e(N_uniq_arm0) <= e(N_nonuniq_arm0)
assert e(N_uniq_arm1) <= e(N_nonuniq_arm1)

* ------------------------------------------------------------
* Test 25: plot option — ITT
*   e(cif) matrix has 3 columns, CIF values in [0,1], non-decreasing
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) plot

local nr = rowsof(e(cif))
assert `nr' > 0
assert colsof(e(cif)) == 3

forvalues i = 1/`nr' {
    assert e(cif)[`i', 2] >= 0 & e(cif)[`i', 2] <= 1
    assert e(cif)[`i', 3] >= 0 & e(cif)[`i', 3] <= 1
}
forvalues i = 2/`nr' {
    assert e(cif)[`i', 2] >= e(cif)[`i'-1, 2] - 1e-10
    assert e(cif)[`i', 3] >= e(cif)[`i'-1, 3] - 1e-10
}

* without plot option, e(cif) should not be returned
seqtte outcome, id(id) time(time) treatment(treatment)
capture matrix list e(cif)
assert _rc != 0

* ------------------------------------------------------------
* Test 26: plot option — PP (unweighted)
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) ///
    estimator(pp) plot

local nr = rowsof(e(cif))
assert `nr' > 0
assert colsof(e(cif)) == 3

forvalues i = 1/`nr' {
    assert e(cif)[`i', 2] >= 0 & e(cif)[`i', 2] <= 1
    assert e(cif)[`i', 3] >= 0 & e(cif)[`i', 3] <= 1
}
forvalues i = 2/`nr' {
    assert e(cif)[`i', 2] >= e(cif)[`i'-1, 2] - 1e-10
    assert e(cif)[`i', 3] >= e(cif)[`i'-1, 3] - 1e-10
}

di as txt "seqtte cscript passed"
