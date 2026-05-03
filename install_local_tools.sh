#!/usr/bin/env bash
set -euo pipefail

# Run from an Rmd with something like:
# system("bash ../install_local_tools.sh")

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_DIR="$ROOT_DIR/tools"
BIN_DIR="$TOOLS_DIR/bin"
MAMBA_DIR="$TOOLS_DIR/.micromamba"
MAMBA_BIN="$MAMBA_DIR/bin/micromamba"
ENV_PREFIX="$TOOLS_DIR/env"
PERL5_DIR="$TOOLS_DIR/perl5"
CLUMPP_URL="https://rosenberglab.stanford.edu/software/CLUMPP_Linux64.1.1.2.tar.gz"
POPOOLATION2_URL="https://sourceforge.net/projects/popoolation2/files/popoolation2_1201.zip/download"

log() {
  printf '[install_local_tools] %s\n' "$*"
}

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    printf 'Missing required command: %s\n' "$1" >&2
    exit 1
  fi
}

write_text_file() {
  local target="$1"
  shift
  cat >"$target" <<EOF
$*
EOF
}

install_micromamba() {
  if [[ -x "$MAMBA_BIN" ]]; then
    return
  fi

  require_cmd curl
  require_cmd tar

  log "Installing micromamba under $MAMBA_DIR"
  mkdir -p "$MAMBA_DIR/bin"
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xj -C "$MAMBA_DIR/bin" --strip-components=1 bin/micromamba
}

ensure_env() {
  install_micromamba

  if [[ ! -d "$ENV_PREFIX" ]]; then
    log "Creating local runtime environment in $ENV_PREFIX"
    "$MAMBA_BIN" create -y -p "$ENV_PREFIX" -c conda-forge -c bioconda \
      gsl perl perl-app-cpanminus openjdk structure
  fi
}

ensure_structure_package() {
  if [[ -x "$ENV_PREFIX/bin/structure" ]]; then
    return
  fi

  log "Installing Bioconda structure into $ENV_PREFIX"
  "$MAMBA_BIN" install -y -p "$ENV_PREFIX" -c conda-forge -c bioconda structure
}

install_perl_module() {
  mkdir -p "$PERL5_DIR"
  log "Installing Perl dependency for Popoolation2"
  export PATH="$ENV_PREFIX/bin:$PATH"
  export PERL5LIB="$PERL5_DIR/lib/perl5${PERL5LIB:+:$PERL5LIB}"
  "$ENV_PREFIX/bin/cpanm" -L "$PERL5_DIR" Text::NSP::Measures::2D::Fisher::twotailed >/dev/null
}

install_popoolation2() {
  local dest="$TOOLS_DIR/popoolation2"
  if [[ -d "$dest/popoolation2_1201" ]]; then
    return
  fi

  require_cmd curl
  require_cmd unzip

  log "Downloading Popoolation2"
  mkdir -p "$dest"
  curl -fsSL "$POPOOLATION2_URL" -o "$TOOLS_DIR/popoolation2.zip"
  unzip -q "$TOOLS_DIR/popoolation2.zip" -d "$dest"
  find "$dest" -type f \( -name '*.pl' -o -name '*.jar' \) -exec chmod 0755 {} +
  rm -f "$TOOLS_DIR/popoolation2.zip"
}

install_npstat() {
  local dest="$TOOLS_DIR/npstat"
  if [[ ! -d "$dest/.git" ]]; then
    require_cmd git
    log "Cloning npstat"
    git clone --depth 1 https://github.com/lucaferretti/npstat.git "$dest"
  fi

  require_cmd gcc

  log "Compiling npstat against local GSL"
  gcc -o "$dest/npstat" "$dest/NPStat-v1.c" \
    -I"$ENV_PREFIX/include" \
    -L"$ENV_PREFIX/lib" \
    -Wl,-rpath,"$ENV_PREFIX/lib" \
    -lgsl -lgslcblas -lm
  chmod 0755 "$dest/npstat"
}

install_clumpp() {
  local dest="$TOOLS_DIR/clumpp"
  if [[ -x "$dest/CLUMPP_Linux64.1.1.2/CLUMPP" ]]; then
    return
  fi

  require_cmd curl
  require_cmd tar

  log "Downloading CLUMPP"
  mkdir -p "$dest"
  curl -fsSL "$CLUMPP_URL" -o "$TOOLS_DIR/CLUMPP_Linux64.1.1.2.tar.gz"
  tar -xzf "$TOOLS_DIR/CLUMPP_Linux64.1.1.2.tar.gz" -C "$dest"
  chmod 0755 "$dest/CLUMPP_Linux64.1.1.2/CLUMPP"
  rm -f "$TOOLS_DIR/CLUMPP_Linux64.1.1.2.tar.gz"
}

install_structure_harvester() {
  local dest="$TOOLS_DIR/structureHarvester"
  if [[ -d "$dest/.git" ]]; then
    return
  fi

  require_cmd git

  log "Cloning structureHarvester"
  git clone --depth 1 https://github.com/dentearl/structureHarvester.git "$dest"
  chmod 0755 "$dest/structureHarvester.py"
}

write_wrapper_scripts() {
  mkdir -p "$BIN_DIR"

  cp "$TOOLS_DIR/npstat/npstat" "$BIN_DIR/npstat"
  cp "$TOOLS_DIR/clumpp/CLUMPP_Linux64.1.1.2/CLUMPP" "$BIN_DIR/CLUMPP"

  write_text_file "$BIN_DIR/structure" "#!/usr/bin/env bash
set -euo pipefail
tool_dir=\"\$(cd \"\$(dirname \"\${BASH_SOURCE[0]}\")/..\" && pwd)\"
binary=\"\$tool_dir/env/bin/structure\"
exec \"\$binary\" \"\$@\""

  write_text_file "$BIN_DIR/structureHarvester" "#!/usr/bin/env bash
set -euo pipefail
tool_dir=\"\$(cd \"\$(dirname \"\${BASH_SOURCE[0]}\")/..\" && pwd)\"
exec \"\$tool_dir/env/bin/python\" \"\$tool_dir/structureHarvester/structureHarvester.py\" \"\$@\""
  cp "$BIN_DIR/structureHarvester" "$BIN_DIR/structureHarvester.py"

  local pop_dir="$TOOLS_DIR/popoolation2/popoolation2_1201"
  local script_path
  for script_path in "$pop_dir"/*.pl "$pop_dir"/export/*.pl "$pop_dir"/indel_filtering/*.pl; do
    [[ -f "$script_path" ]] || continue
    write_text_file "$BIN_DIR/$(basename "$script_path")" "#!/usr/bin/env bash
set -euo pipefail
tool_dir=\"\$(cd \"\$(dirname \"\${BASH_SOURCE[0]}\")/..\" && pwd)\"
export PERL5LIB=\"\$tool_dir/perl5/lib/perl5\${PERL5LIB:+:\$PERL5LIB}\"
export PATH=\"\$tool_dir/env/bin:\$PATH\"
exec \"$script_path\" \"\$@\""
  done

  chmod 0755 "$BIN_DIR"/*
}

write_env_script() {
  write_text_file "$TOOLS_DIR/env.sh" "#!/usr/bin/env bash
tool_dir=\"\$(cd \"\$(dirname \"\${BASH_SOURCE[0]}\")\" && pwd)\"
export PATH=\"\$tool_dir/bin:\$tool_dir/env/bin:\$PATH\"
export LD_LIBRARY_PATH=\"\$tool_dir/env/lib\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}\"
export PERL5LIB=\"\$tool_dir/perl5/lib/perl5\${PERL5LIB:+:\$PERL5LIB}\""
  chmod 0755 "$TOOLS_DIR/env.sh"
}

main() {
  require_cmd bash
  require_cmd curl
  require_cmd tar
  require_cmd unzip
  require_cmd git
  require_cmd python3
  require_cmd perl
  mkdir -p "$TOOLS_DIR"

  ensure_env
  ensure_structure_package
  install_perl_module
  install_popoolation2
  install_npstat
  install_clumpp
  install_structure_harvester
  write_wrapper_scripts
  write_env_script

  log "Local tools installed under $TOOLS_DIR"
  log "From Rmd, use Sys.setenv(PATH = paste('$BIN_DIR', Sys.getenv('PATH'), sep = ':'), LD_LIBRARY_PATH = paste('$ENV_PREFIX/lib', Sys.getenv('LD_LIBRARY_PATH'), sep = ':'), PERL5LIB = paste('$PERL5_DIR/lib/perl5', Sys.getenv('PERL5LIB'), sep = ':'))"
  log "Or in a shell: source $TOOLS_DIR/env.sh"
}

main "$@"