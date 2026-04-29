# Install renv + packages during image build.
# - renv cache/library are stored on the image under /home/rstudio/renv
# - the renv project lives in /home/rstudio/work (where users typically mount their project)

options(renv.consent = TRUE)
write("options(renv.consent = TRUE)", file = "/home/rstudio/.Rprofile", append = TRUE)
write(
  paste0(
    "\n",
    "# Ensure renv project library is on the default library paths\n",
    "local({\n",
    "  if (requireNamespace('renv', quietly = TRUE)) {\n",
    "    project <- Sys.getenv('RENV_PROJECT', '/home/rstudio/work')\n",
    "    lib <- tryCatch(renv::paths$library(project = project), error = function(e) NULL)\n",
    "    if (!is.null(lib) && dir.exists(lib)) .libPaths(c(lib, .libPaths()))\n",
    "  }\n",
    "})\n"
  ),
  file = "/home/rstudio/.Rprofile",
  append = TRUE
)

cran <- Sys.getenv("CRAN", "https://cloud.r-project.org")
options(repos = c(CRAN = cran))

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

project <- Sys.getenv("RENV_PROJECT", "/home/rstudio/work")
dir.create(project, recursive = TRUE, showWarnings = FALSE)

# Initialize the project. A brand-new directory may not contain a lockfile yet.
renv::init(project = project, force = TRUE)

lockfile <- file.path(project, "renv.lock")
if (file.exists(lockfile)) {
  renv::restore(project = project, confirm = FALSE)
}

cranPackages <- c(
  "optparse",
  "tidyverse",
  "cowplot",
  "data.table",
  "stringr",
  "locfit",
  "doMC",
  "gridExtra",
  "vcfR",
  "agricolae",
  "car",
  "lme4",
  "argparse",
  "adegenet",
  "ade4",
  "adegraphics",
  "pegas",
  "devtools",
  "vegan",
  "graph4lg",
  "reshape2",
  "rworldmap",
  "poolfstat"
)

renv::install(cranPackages, project = project)

biocPackages <- c(
  "bioc::ggplot2",
  "bioc::Biostrings"
)

renv::install(biocPackages, project = project)
renv::install("jgx65/hierfstat", project = project)
renv::install("pievos101/PopGenome", project = project)

# Make the installed set explicit and reusable.
# - Use type = "all" so we capture the packages we just installed even if the
#   project doesn't contain any R code yet (implicit snapshots can be nearly empty).
# - Snapshot writes renv.lock in the project directory.
# - Copy renv.lock to /home/rstudio/renv.lock so it remains accessible even when
#   /home/rstudio/work is bind-mounted.
renv::snapshot(project = project, prompt = FALSE, type = "all")

lockSrc <- file.path(project, "renv.lock")
lockDst <- "/home/rstudio/renv.lock"
if (file.exists(lockSrc)) {
  file.copy(lockSrc, lockDst, overwrite = TRUE)
}

lockDst2 <- "/opt/practical3/renv.lock"
dir.create(dirname(lockDst2), recursive = TRUE, showWarnings = FALSE)
if (file.exists(lockSrc)) {
  file.copy(lockSrc, lockDst2, overwrite = TRUE)
}
