cscript seqtte adofile seqtte

* ------------------------------------------------------------
* Generate a synthetic person-period dataset for testing
* ------------------------------------------------------------
* 200 individuals, 10 periods (0-9)
* First 100: treated from period 0 (ever-treated)
* Last 100:  never treated
* Outcome: approximately 8% per-period hazard, slightly lower for treated

clear
set seed 42
set obs 200

gen id = _n
gen trt_group = (id <= 100)  // 1 = treated, 0 = untreated
gen age_group = mod(id, 3)   // synthetic covariate

expand 10
bysort id: gen time = _n - 1

gen treatment = trt_group * (time == 0)
bysort id (time): replace treatment = max(treatment, treatment[_n-1])

* Period-specific outcome: 1 in the period the event occurs, 0 otherwise
* Slightly lower hazard for treated group
gen double u = runiform()
gen outcome_any = (u < 0.08 - 0.03 * trt_group)

* Keep only the first event per person (drop post-event rows)
bysort id (time): gen cumev = sum(outcome_any)
drop if cumev > 1
replace outcome_any = 0 if cumev == 1 & u >= (0.08 - 0.03 * trt_group)

* outcome = 1 only in the period the event first occurs
gen outcome = outcome_any

drop trt_group cumev u outcome_any

sort id time

* ------------------------------------------------------------
* Test 1: basic ITT run, no covariates
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment)

* Stored results exist
assert e(N) > 0
assert e(N_indiv) == 200
assert e(N_orig) > 0
assert e(N_exp) >= e(N_orig)
assert "`e(cmd)'" == "seqtte"
assert "`e(estimator)'" == "itt"

* ------------------------------------------------------------
* Test 2: with a covariate
* ------------------------------------------------------------
seqtte outcome, id(id) time(time) treatment(treatment) covariates(age_group)

assert e(N) > 0
assert "`e(cmd)'" == "seqtte"

* ------------------------------------------------------------
* Test 3: if/in restriction
* ------------------------------------------------------------
seqtte outcome if time <= 7, id(id) time(time) treatment(treatment)

assert e(N) > 0

di as txt "seqtte cscript passed"
