FROM rocker/rstudio:4.5.3

ARG DEBIAN_FRONTEND=noninteractive

# RStudio Server + renv with a pre-populated renv cache/library.
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    git \
    locales \
    build-essential \
    gfortran \
    pkg-config \
    cmake \
    gdal-bin \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libglpk-dev \
    libx11-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
  && rm -rf /var/lib/apt/lists/*

RUN sed -i 's/^# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
  && locale-gen
ENV LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

# Similar to the reference image (which pins CRAN), but simplified:
# use a Posit Package Manager CRAN snapshot by date.
ENV CRAN=https://packagemanager.posit.co/cran/2026-03-12
ENV RENV_CONFIG_REPOS_OVERRIDE=${CRAN}

# Keep renv state inside /home/rstudio (so it remains on the image even when /home/rstudio/work is bind-mounted).
# NOTE: Some platforms (e.g. Renku) mount /home/rstudio, which would hide anything stored
# under that path. Store renv state under /opt so installed packages remain available.
ENV RENV_PATHS_ROOT=/opt/renv
ENV RENV_PATHS_CACHE=/opt/renv/cache
ENV RENV_PATHS_LIBRARY_ROOT=/opt/renv/library
ENV RENV_PROJECT=/home/rstudio/work

# Ensure the renv project library is available in ALL R sessions (including when
# running R as root/non-rstudio). This avoids needing manual .libPaths() edits.
# NOTE: In rocker images, R reads the site profile from $R_HOME/etc/Rprofile.site.
RUN SITE_PROFILE=/usr/local/lib/R/etc/Rprofile.site; \
    if [ -f "$SITE_PROFILE" ] && grep -q "EVO_PRACTICAL3_RENV_LIBPATHS" "$SITE_PROFILE"; then \
      true; \
    else \
      mkdir -p "$(dirname "$SITE_PROFILE")"; \
      printf '%s\n' \
        '' \
        '# EVO_PRACTICAL3_RENV_LIBPATHS' \
        'local({' \
        '  if (requireNamespace("renv", quietly = TRUE)) {' \
        '    project <- Sys.getenv("RENV_PROJECT", "/home/rstudio/work")' \
        '    dir.create(project, recursive = TRUE, showWarnings = FALSE)' \
        '    lib <- tryCatch(renv::paths$library(project = project), error = function(e) NULL)' \
        '    if (!is.null(lib)) {' \
        '      dir.create(lib, recursive = TRUE, showWarnings = FALSE)' \
        '      .libPaths(c(lib, .libPaths()))' \
        '    }' \
        '  }' \
        '})' \
        >> "$SITE_PROFILE"; \
    fi

RUN R -q -e 'options(repos = c(CRAN = Sys.getenv("CRAN", "https://cloud.r-project.org"))); install.packages(c("renv","remotes","BiocManager"))'

RUN mkdir -p /home/rstudio/.cache/R/renv /home/rstudio/work /opt/renv /opt/practical3 \
  && chown -R rstudio:rstudio /home/rstudio/.cache /home/rstudio/work /opt/renv /opt/practical3

COPY --chown=root:root .renv_install.sh /tmp/.renv_install.sh
COPY --chown=root:root install.R /tmp/install.R

RUN chmod +x /tmp/.renv_install.sh \
  && su - rstudio -c 'bash /tmp/.renv_install.sh' \
  && rm -f /tmp/.renv_install.sh /tmp/install.R

# Put a short usage README in the user's home directory.
COPY --chown=rstudio:rstudio README.home.Rmd /home/rstudio/README.Rmd

# Also keep a copy outside /home (some platforms mount /home/rstudio).
COPY --chown=rstudio:rstudio README.home.Rmd /opt/practical3/README.Rmd

# Make renv-installed packages available on the default library path.
# (Symlinks are small and avoid duplicating the package trees.)
RUN RENV_LIB="$(find /opt/renv/library -mindepth 4 -maxdepth 4 -type d -name 'x86_64-pc-linux-gnu' | head -n 1)" \
  && test -n "$RENV_LIB" \
  && test -d "$RENV_LIB" \
  && for pkg_dir in "$RENV_LIB"/*; do \
       if [ -d "$pkg_dir" ] && [ -f "$pkg_dir/DESCRIPTION" ]; then \
         pkg_name="$(basename "$pkg_dir")"; \
         if [ ! -e "/usr/local/lib/R/site-library/$pkg_name" ]; then \
           ln -s "$pkg_dir" "/usr/local/lib/R/site-library/$pkg_name"; \
         fi; \
       fi; \
     done

EXPOSE 8787
