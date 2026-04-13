#!/usr/bin/env bash
set -euo pipefail

# This script is executed during docker build.
# It prepares a renv project in /home/rstudio/work while keeping the renv cache/library
# outside the bind-mounted project directory (so packages remain on the image).

: "${CRAN:=https://packagemanager.posit.co/cran/2026-03-12}"
: "${RENV_PROJECT:=/home/rstudio/work}"

export RENV_CONFIG_REPOS_OVERRIDE="${CRAN}"

# Keep renv state on the image (NOT inside the bind-mounted project dir).
export RENV_PATHS_ROOT="/opt/renv"
export RENV_PATHS_CACHE="/opt/renv/cache"
export RENV_PATHS_LIBRARY_ROOT="/opt/renv/library"

mkdir -p /opt/renv /home/rstudio/work

R -q -f /tmp/install.R
