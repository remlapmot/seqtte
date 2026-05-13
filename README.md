# seqtte

Stata package for sequential target trial emulation.

## Installation

```stata
net install seqtte, from("https://raw.githubusercontent.com/remlapmot/seqtte/main/") replace
```

To check for an update run

```stata
adoupdate
```

To update the package run

```stata
adoupdate seqtte, update
```

To uninstall the package run

```sh
ado uninstall seqtte
```

## Commands

| Command | Description |
|---------|-------------|
| `seqtte` | Intent-to-treat (ITT) and per-protocol (PP) estimators via sequential trial emulation |

## Usage

```stata
seqtte outcomevar, id(varname) time(varname) treatment(varname) [covariates(varlist)] [estimator(itt|pp)] [wdenominator(varlist)] [wnumerator(varlist)] [truncation(#)]
```

The input data should be in long (person-period) format with:
- one row per individual per time period
- `outcomevar` equal to 1 only in the period the event first occurs
- consecutive integer values for `time`
- binary `treatment` (0/1)

## For developers

### Running certification scripts

From the `cscripts/` directory in Stata:

```stata
do master
```

Or from top level of repo

```sh
just test
```

## Authors

* Tom Palmer, University of Bristol, Bristol, UK.
* Michalis Katsoulis, UCL, London, UK.

## References

Hernán MA, Robins JM (2016). Using Big Data to Emulate a Target Trial When a Randomized Trial Is Not Available. *American Journal of Epidemiology* 183(8): 758–764.
