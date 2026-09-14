# Simple Quarto website for seqtte package

Install the R dependencies with

```r
if (!requireNamespace("Statamarkdown")) install.packages("Statamarkdown")
if (!requireNamespace("quarto")) install.packages("quarto")
```

To install the Stata dependencies run in Stata
```stata
do developer.do
```

To build the website, either open this directory as an RStudio project and then in the *Build* pane click *Build Website*, or run in R
```r
quarto::quarto_render()
```

Or in a shell
```sh
quarto render
```

Or from the top level of the repo, which also copies the rendered site into `docs/`
```sh
just render
```

To rerender a single file run say
```sh
quarto render index.qmd
```
