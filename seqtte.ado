*! version 0.5.0  20may2026  Tom Palmer
program seqtte, eclass
    version 16

    syntax varlist(min=1 max=1 numeric) [if] [in], ///
        id(varname numeric) ///
        time(varname numeric) ///
        treatment(varname numeric) ///
        [covariates(varlist fv) ///
         ESTIMator(string) ///
         wdenominator(varlist fv) ///
         wnumerator(varlist fv) ///
         TRUNCation(real 25) ///
         SELECTIONrandom ///
         SELECTIONsample(real 0.5) ///
         SEEd(integer -1) ///
         BOOTstrap(integer 0)]

    local outcome `varlist'

    // Default and validate estimator
    if "`estimator'" == "" local estimator "itt"
    local estimator = lower("`estimator'")
    if !inlist("`estimator'", "itt", "pp") {
        di as err "estimator() must be itt or pp"
        exit 198
    }
    if "`estimator'" == "pp" & "`wdenominator'" == "" {
        di as err "wdenominator() required for per-protocol estimation"
        exit 198
    }

    // Validate selection_random option
    if "`selectionrandom'" != "" {
        if `selectionsample' <= 0 | `selectionsample' > 1 {
            di as err "selectionsample() must be in (0, 1]"
            exit 198
        }
    }

    // Validate bootstrap option
    if `bootstrap' < 0 {
        di as err "bootstrap() must be a non-negative integer"
        exit 198
    }
    local do_bs = (`bootstrap' > 0)

    // Capture variable count before tempvars are created
    local k_orig = c(k)
    local seed_set 0

    preserve

    marksample touse
    // markout needs base variable names, not factor-variable notation
    local cov_base `covariates'
    if "`cov_base'" != "" {
        fvrevar `cov_base', list
        local cov_base `r(varlist)'
    }
    markout `touse' `id' `time' `treatment' `cov_base'
    if "`estimator'" == "pp" {
        local wden_base `wdenominator'
        local wnum_base `wnumerator'
        if "`wden_base'" != "" {
            fvrevar `wden_base', list
            local wden_base `r(varlist)'
        }
        if "`wnum_base'" != "" {
            fvrevar `wnum_base', list
            local wnum_base `r(varlist)'
        }
        markout `touse' `wden_base' `wnum_base'
    }
    keep if `touse'

    sort `id' `time'

    qui count
    local n_orig = r(N)
    qui levelsof `id'
    local n_indiv = r(r)

    di as txt _n "Original dataset: " %12.0fc `n_orig' " observations, " `k_orig' " variables"

    // Event time per individual
    tempvar evttime evttime_max
    qui gen double `evttime' = .
    qui replace `evttime' = `time' if `outcome' == 1
    by `id': egen double `evttime_max' = max(`evttime')

    // Eligibility: not treated in any earlier period
    tempvar A_lag eligible
    by `id': gen byte `A_lag' = `treatment'[_n-1]
    gen byte `eligible' = 1
    qui replace `eligible' = 0 if `A_lag' == 1
    by `id': replace `eligible' = 0 ///
        if `eligible'[_n-1] == 0 & `id' == `id'[_n-1]

    qui count if `eligible' == 1
    local n_elig = r(N)
    di as txt "Eligible observations: " %12.0fc `n_elig' " (" `n_indiv' " individuals)"

    // Trial = calendar time at trial entry
    tempvar trial
    qui gen long `trial' = `time'

    // ----- PP: weight models and per-period snapshots (before expand) -----
    if "`estimator'" == "pp" {

        di as txt _n "Fitting weight models..."

        qui sum `time'
        local t_min = r(min)
        local t_max = r(max)

        // Denominator: P(A | A_lag, time + time² + time³, wdenominator)
        // Guard each fit: if treatment does not vary within the A_lag stratum
        // (e.g. monotonic, non-reversible treatment, where everyone with
        // A_lag == 1 also has treatment == 1), logistic exits with r(2000)
        // "outcome does not vary". In that case the stratum's weight
        // contribution is 1, so set the predicted probability to the constant
        // value of treatment instead of fitting a model.
        tempvar p_d0 p_d1

        qui sum `treatment' if `A_lag' == 0, meanonly
        if r(N) == 0 | r(min) == r(max) {
            qui gen double `p_d0' = cond(r(N) == 0, 1, r(mean))
        }
        else {
            qui logistic `treatment' ///
                c.`time' c.`time'#c.`time' c.`time'#c.`time'#c.`time' ///
                `wdenominator' if `A_lag' == 0
            qui predict double `p_d0'
        }

        qui sum `treatment' if `A_lag' == 1, meanonly
        if r(N) == 0 | r(min) == r(max) {
            qui gen double `p_d1' = cond(r(N) == 0, 1, r(mean))
        }
        else {
            qui logistic `treatment' ///
                c.`time' c.`time'#c.`time' c.`time'#c.`time'#c.`time' ///
                `wdenominator' if `A_lag' == 1
            qui predict double `p_d1'
        }

        // Numerator: P(A | A_lag, time + time² + time³, wnumerator)
        local stabilized = ("`wnumerator'" != "")
        if `stabilized' {
            tempvar p_n0 p_n1

            qui sum `treatment' if `A_lag' == 0, meanonly
            if r(N) == 0 | r(min) == r(max) {
                qui gen double `p_n0' = cond(r(N) == 0, 1, r(mean))
            }
            else {
                qui logistic `treatment' ///
                    c.`time' c.`time'#c.`time' c.`time'#c.`time'#c.`time' ///
                    `wnumerator' if `A_lag' == 0
                qui predict double `p_n0'
            }

            qui sum `treatment' if `A_lag' == 1, meanonly
            if r(N) == 0 | r(min) == r(max) {
                qui gen double `p_n1' = cond(r(N) == 0, 1, r(mean))
            }
            else {
                qui logistic `treatment' ///
                    c.`time' c.`time'#c.`time' c.`time'#c.`time'#c.`time' ///
                    `wnumerator' if `A_lag' == 1
                qui predict double `p_n1'
            }
        }

        // IPW weights (reset to 1 at each person's first record)
        tempvar ipw
        if `stabilized' {
            qui gen double `ipw' = (1 - `p_n0') / (1 - `p_d0') if `treatment' == 0
            qui replace `ipw' = `p_n1' / `p_d1'                if `treatment' == 1
        }
        else {
            qui gen double `ipw' = 1 / (1 - `p_d0') if `treatment' == 0
            qui replace `ipw' = 1 / `p_d1'          if `treatment' == 1
        }
        qui replace `ipw' = 1 if `id' != `id'[_n-1]

        // Per-calendar-period censoring and weight snapshots
        local loop_start = `t_min' + 1
        forvalues i = `loop_start'/`t_max' {
            tempvar _c`i' _w`i'
            qui gen byte   `_c`i'' = 0
            qui gen double `_w`i'' = 1
            qui replace `_c`i'' = 1 ///
                if `trial' == `i' & `treatment' != `treatment'[_n-1] ///
                & `id' == `id'[_n-1]
            qui replace `_w`i'' = `ipw' if `trial' == `i'
        }

        // Aggregate snapshots to person level
        forvalues i = `loop_start'/`t_max' {
            tempvar _ca`i' _wa`i'
            by `id': egen byte   `_ca`i'' = max(`_c`i'')
            by `id': egen double `_wa`i'' = max(`_w`i'')
        }

        di as txt "Weight models fitted"
    }

    // ----- Expand each row to cover remaining follow-up -----
    di as txt _n "Expanding data..."

    tempvar max_t n_expand
    by `id': gen long `max_t' = `time'[_N]
    qui gen long `n_expand' = `max_t' - `time' + 1
    expand `n_expand'

    sort `id' `trial'

    tempvar time_in_trial
    by `id' `trial': gen long `time_in_trial' = _n

    drop if `eligible' == 0

    qui count
    local n_exp = r(N)
    di as txt "Expanded dataset: " %12.0fc `n_exp' " observations"

    // ----- PP: fill censoring and cumulative weights into expanded data -----
    if "`estimator'" == "pp" {
        tempvar censored wt
        qui gen byte   `censored' = 0
        qui gen double `wt'       = 1

        local loop_start = `t_min' + 1
        forvalues i = `loop_start'/`t_max' {
            local cond "`i' == `time_in_trial' + `trial' - 1 & `time_in_trial' > 1"
            qui replace `censored' = `_ca`i'' if `cond'
            qui replace `wt'       = `_wa`i'' if `cond'
        }

        // Propagate censoring forward within each (id, trial)
        sort `id' `trial' `time_in_trial'
        qui replace `censored' = 1 ///
            if `censored'[_n-1] == 1 ///
            & `id'    == `id'[_n-1] ///
            & `trial' == `trial'[_n-1]

        // Cumulative product of weights within (id, trial)
        tempvar wt_cum
        qui gen double `wt_cum' = `wt'
        by `id' `trial': ///
            replace `wt_cum' = `wt_cum' * `wt_cum'[_n-1] if _n > 1

        // Truncate extreme weights
        qui replace `wt_cum' = `truncation' ///
            if `wt_cum' > `truncation' & !missing(`wt_cum')

        qui count if `censored' == 0
        local n_pp = r(N)
        di as txt "Post-censoring dataset: " %12.0fc `n_pp' " observations"
    }

    // ----- Random selection of control-arm (id, trial) pairs -----
    if "`selectionrandom'" != "" {

        if `seed' != -1 & !`seed_set' {
            set seed `seed'
            local seed_set 1
        }

        // One random draw per (id, trial) pair; treatment is constant within pairs
        tempvar rnd_pair
        bysort `id' `trial' (`time_in_trial'): ///
            gen double `rnd_pair' = runiform() if _n == 1
        bysort `id' `trial' (`time_in_trial'): ///
            replace `rnd_pair' = `rnd_pair'[1]

        // Keep all treated-arm pairs; Bernoulli-sample control-arm pairs
        keep if `treatment' == 1 | (`treatment' == 0 & `rnd_pair' <= `selectionsample')

        qui count
        local n_sel = r(N)
        di as txt "After random selection (" %5.3f `selectionsample' ///
            " of control-arm pairs): " %12.0fc `n_sel' " observations"
    }

    // ----- Outcome and follow-up time -----
    tempvar event fu_time
    qui gen byte `event' = 0
    qui replace `event' = 1 ///
        if `evttime_max' == `time_in_trial' + `trial' - 1 ///
        & !missing(`evttime_max')
    qui gen long `fu_time' = `time_in_trial' - 1

    // ----- Readable names for the regression table -----
    // The data are preserved/restored, so these renames are discarded after the
    // fit. If a name is already in use in the caller's data we keep the
    // temporary name to avoid a collision.
    capture confirm new variable event
    if !_rc {
        rename `event' event
        local event event
    }
    capture confirm new variable followup
    if !_rc {
        rename `fu_time' followup
        local fu_time followup
    }
    capture confirm new variable trial
    if !_rc {
        rename `trial' trial
        local trial trial
    }

    // ----- Bootstrap -----
    if `do_bs' {

        if `seed' != -1 & !`seed_set' {
            set seed `seed'
            local seed_set 1
        }

        // Save fully processed dataset (weights, outcome, renaming applied)
        tempfile bsdata
        qui save `bsdata'

        di as txt _n "Running " `bootstrap' " bootstrap replicates..."

        tempname bs_b
        matrix `bs_b' = J(`bootstrap', 1, .)
        local bs_ok = 0
        // newid declared once; bsample creates it fresh each iteration after reload
        tempvar bs_newid

        forvalues b = 1/`bootstrap' {
            // reload from tempfile instead of nested preserve (r(621))
            qui use `bsdata', clear
            cap {
                bsample `n_indiv', cluster(`id') idcluster(`bs_newid')

                if "`estimator'" == "itt" {
                    qui logistic `event' `treatment' ///
                        c.`fu_time'##c.`fu_time' ///
                        c.`trial'##c.`trial' ///
                        `covariates', cluster(`bs_newid')
                }
                else {
                    qui logistic `event' `treatment' ///
                        c.`fu_time'##c.`fu_time' ///
                        c.`trial'##c.`trial' ///
                        `covariates' ///
                        if `censored' == 0 [pweight = `wt_cum'], cluster(`bs_newid')
                }

                matrix `bs_b'[`b', 1] = _b[`treatment']
                local bs_ok = `bs_ok' + 1
            }
        }

        // Reload processed data so the main regression has the correct dataset
        qui use `bsdata', clear

        // Bootstrap SE and percentile CI (95%) from log-OR distribution
        // svmat places values in obs 1..B; remaining obs are missing — both
        // sum and _pctile skip missing values automatically.
        svmat double `bs_b', names(_seqtte_bs_)
        qui sum _seqtte_bs_1
        local bs_se = r(sd)
        qui _pctile _seqtte_bs_1, p(2.5 97.5)
        local bs_ll = r(r1)
        local bs_ul = r(r2)
        drop _seqtte_bs_1

        di as txt "Bootstrap complete: " `bs_ok' "/" `bootstrap' " replicates succeeded"
        di as txt "Bootstrap SE (log-OR):      " %7.4f `bs_se'
        di as txt "Bootstrap 95% CI (log-OR): [" %7.4f `bs_ll' ", " %7.4f `bs_ul' "]"
    }

    // ----- Pooled logistic regression -----
    di as txt _n "Fitting " upper("`estimator'") " model..."

    if "`estimator'" == "itt" {
        logistic `event' `treatment' ///
            c.`fu_time'##c.`fu_time' ///
            c.`trial'##c.`trial' ///
            `covariates', cluster(`id')
    }
    else {
        logistic `event' `treatment' ///
            c.`fu_time'##c.`fu_time' ///
            c.`trial'##c.`trial' ///
            `covariates' ///
            if `censored' == 0 [pweight = `wt_cum'], cluster(`id')
    }

    restore

    ereturn scalar N_indiv   = `n_indiv'
    ereturn scalar N_orig    = `n_orig'
    ereturn scalar N_exp     = `n_exp'
    if "`selectionrandom'" != "" {
        ereturn scalar N_sel            = `n_sel'
        ereturn scalar selection_sample = `selectionsample'
    }
    if `do_bs' {
        ereturn scalar N_boot = `bs_ok'
        ereturn scalar bs_se  = `bs_se'
        ereturn scalar bs_ll  = `bs_ll'
        ereturn scalar bs_ul  = `bs_ul'
        ereturn matrix bs_b   = `bs_b'
    }
    ereturn local  estimator "`estimator'"
    ereturn local  cmd       "seqtte"
end
