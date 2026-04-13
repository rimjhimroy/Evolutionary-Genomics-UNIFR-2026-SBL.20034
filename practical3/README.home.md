# Installing new R packages with renv

This image uses `renv` to manage packages for the project at `/home/rstudio/work`.

## Install packages with `renv::install()` (this installs into the project library and tracks versions):

### CRAN (examples not preinstalled)

```r
renv::install(c(
  "pcadapt",   # selection scans
  "poolfstat" # Pool-Seq statistics
))
```

### Bioconductor (examples not preinstalled)

```r
renv::install(c(
  "bioc::gdsfmt",
  "bioc::SNPRelate"
))
```

### GitHub (example not preinstalled)

```r
renv::install("grunwaldlab/poppr")
```

3) Update the lockfile so others (or future you) can reproduce the environment:

```r
renv::snapshot()
```


## If you prefer base R installers (not recommended with renv)

These install into the active library paths but won’t reliably update `renv.lock` unless you snapshot afterwards.

```r
install.packages("pcadapt")
BiocManager::install("SNPRelate")
remotes::install_github("grunwaldlab/poppr")
renv::snapshot()
```
