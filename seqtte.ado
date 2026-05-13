*! version 0.2.0  13may2026  Tom Palmer
program seqtte, eclass
    version 16

    syntax varlist(min=1 max=1 numeric) [if] [in], ///
        id(varname numeric) ///
        time(varname numeric) ///
        treatment(varname numeric) ///
        [covariates(varlist) ///
         ESTIMator(string) ///
         wdenominator(varlist) ///
         wnumerator(varlist) ///
         TRUNCation(real 25)]

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

    preserve

    marksample touse
    markout `touse' `id' `time' `treatment' `covariates'
    if "`estimator'" == "pp" markout `touse' `wdenominator' `wnumerator'
    keep if `touse'

    sort `id' `time'

    qui count
    local n_orig = r(N)
    qui levelsof `id'
    local n_indiv = r(r)

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

    // Trial = calendar time at trial entry
    tempvar trial
    qui gen long `trial' = `time'

    // ----- PP: weight models and per-period snapshots (before expand) -----
    if "`estimator'" == "pp" {

        qui sum `time'
        local t_min = r(min)
        local t_max = r(max)

        // Denominator: P(A | A_lag, time + time² + time³, wdenominator)
        tempvar p_d0 p_d1
        qui logistic `treatment' ///
            c.`time' c.`time'#c.`time' c.`time'#c.`time'#c.`time' ///
            `wdenominator' if `A_lag' == 0
        qui predict double `p_d0'

        qui logistic `treatment' ///
            c.`time' c.`time'#c.`time' c.`time'#c.`time'#c.`time' ///
            `wdenominator' if `A_lag' == 1
        qui predict double `p_d1'

        // Numerator: P(A | A_lag, time + time² + time³, wnumerator)
        // If wnumerator not given, use unstabilized weights (numerator = 1)
        local stabilized = ("`wnumerator'" != "")
        if `stabilized' {
            tempvar p_n0 p_n1
            qui logistic `treatment' ///
                c.`time' c.`time'#c.`time' c.`time'#c.`time'#c.`time' ///
                `wnumerator' if `A_lag' == 0
            qui predict double `p_n0'

            qui logistic `treatment' ///
                c.`time' c.`time'#c.`time' c.`time'#c.`time'#c.`time' ///
                `wnumerator' if `A_lag' == 1
            qui predict double `p_n1'
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
    }

    // ----- Expand each row to cover remaining follow-up -----
    tempvar max_t n_expand
    by `id': gen long `max_t' = `time'[_N]
    qui gen long `n_expand' = `max_t' - `time' + 1
    expand `n_expand'

    sort `id' `trial'

    tempvar time_in_trial
    by `id' `trial': gen long `time_in_trial' = _n

    drop if `eligible' == 0

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
    }

    // ----- Outcome and follow-up time -----
    tempvar event fu_time
    qui gen byte `event' = 0
    qui replace `event' = 1 ///
        if `evttime_max' == `time_in_trial' + `trial' - 1 ///
        & !missing(`evttime_max')
    qui gen long `fu_time' = `time_in_trial' - 1

    qui count
    local n_exp = r(N)

    // ----- Pooled logistic regression -----
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
    ereturn local  estimator "`estimator'"
    ereturn local  cmd       "seqtte"
end
