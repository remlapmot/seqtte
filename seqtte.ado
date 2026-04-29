*! version 0.1.0  29apr2026  Tom Palmer
program seqtte, eclass
    version 16

    syntax varlist(min=1 max=1 numeric) [if] [in], ///
        id(varname numeric) ///
        time(varname numeric) ///
        treatment(varname numeric) ///
        [covariates(varlist)]

    local outcome `varlist'

    preserve

    marksample touse
    markout `touse' `id' `time' `treatment' `covariates'
    keep if `touse'

    sort `id' `time'

    // Count original sample
    qui count
    local n_orig = r(N)
    qui levelsof `id'
    local n_indiv = r(r)

    // Step 1: Record event time per individual
    // outcome should be 1 only in the period the event occurs
    tempvar evttime evttime_max
    qui gen double `evttime' = .
    qui replace `evttime' = `time' if `outcome' == 1
    by `id': egen double `evttime_max' = max(`evttime')

    // Step 2: Eligibility — not treated in any earlier period
    tempvar A_lag eligible
    by `id': gen byte `A_lag' = `treatment'[_n-1]
    gen byte `eligible' = 1
    qui replace `eligible' = 0 if `A_lag' == 1
    by `id': replace `eligible' = 0 ///
        if `eligible'[_n-1] == 0 & `id' == `id'[_n-1]

    // Step 3: Trial variable — calendar time at trial entry
    tempvar trial
    qui gen long `trial' = `time'

    // Step 4: Expand each row to cover remaining follow-up
    tempvar max_t n_expand
    by `id': gen long `max_t' = `time'[_N]
    qui gen long `n_expand' = `max_t' - `time' + 1
    expand `n_expand'

    sort `id' `trial'

    // Step 5: Index within trial (1 = trial entry)
    tempvar time_in_trial
    by `id' `trial': gen long `time_in_trial' = _n

    // Step 6: Drop ineligible trial entries
    drop if `eligible' == 0

    // Step 7: Outcome — event occurs at this analysis time point
    // Calendar time = trial + time_in_trial - 1
    tempvar event
    qui gen byte `event' = 0
    qui replace `event' = 1 ///
        if `evttime_max' == `time_in_trial' + `trial' - 1 ///
        & !missing(`evttime_max')

    // Follow-up time within trial (0-indexed from trial entry)
    tempvar fu_time
    qui gen long `fu_time' = `time_in_trial' - 1

    // Count expanded sample
    qui count
    local n_exp = r(N)

    // Step 8: Pooled logistic regression (ITT)
    // Quadratic polynomials for follow-up time and trial number
    logistic `event' `treatment' ///
        c.`fu_time'##c.`fu_time' ///
        c.`trial'##c.`trial' ///
        `covariates', cluster(`id')

    restore

    ereturn scalar N_indiv = `n_indiv'
    ereturn scalar N_orig   = `n_orig'
    ereturn scalar N_exp    = `n_exp'
    ereturn local  estimator "itt"
    ereturn local  cmd       "seqtte"
end
