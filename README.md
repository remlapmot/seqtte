# seqtte

Stata package for sequential target trial emulation.

## Installation

```stata
net install seqtte, from("https://raw.githubusercontent.com/remlapmot/seqtte/main/") replace
```

### If installation fails with a permission error

On managed machines Stata's default installation directory (`PLUS`) is sometimes not writable, and `net install` fails with an error such as `cannot write in directory ...` (`r(603)`). In that case redirect the installation to a directory you can write to using `net set ado`.

**Option 1 — install into your current working directory.** In Stata, `cd` to a folder you know you can write to (for example your Documents folder or your project folder), then run:

```stata
net set ado "`c(pwd)'"
net install seqtte, from("https://raw.githubusercontent.com/remlapmot/seqtte/main/") replace
```

Stata always searches the current working directory for programs, so `seqtte` will work whenever Stata's working directory is that folder. To use it from other folders as well, run `adopath + "<that folder>"` (per session, or add that line to your `profile.do`).

**Option 2 — install into your PERSONAL ado directory.** Run `sysdir` to see where `PERSONAL` points on your machine and create that folder if it does not exist (including any missing parent folders — note Stata's `mkdir` only creates one level at a time, so it may be easier in File Explorer/Finder). Then run:

```stata
// Find PERSONAL location with either
sysdir
// or
adopath
// or
di c(sysdir_personal)
// create the directory if it does not exist

net set ado PERSONAL
net install seqtte, from("https://raw.githubusercontent.com/remlapmot/seqtte/main/") replace
```

`PERSONAL` is on Stata's search path in every session, so the package will then be found automatically.

Note that `net set ado` only lasts for the current Stata session; subsequent `net install`/`adoupdate` runs in a new session will need it set again.

After installation has succeeded launch the helpfile with

```stata
help seqtte
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

```stata
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

Hernán MA, Robins JM (2016). Using Big Data to Emulate a Target Trial When a Randomized Trial Is Not Available. *American Journal of Epidemiology* 183(8): 758–764. <https://doi.org/10.1093/aje/kwv254>
