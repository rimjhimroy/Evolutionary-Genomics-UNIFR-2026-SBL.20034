#!/usr/bin/env bash
# =============================================================================
# setup_practical2_tools.sh
#
# Installs the command-line tools needed for Practical 2:
#   - samtools      (SAM/BAM manipulation)
#   - bwa           (short-read alignment)
#   - java          (OpenJDK – required by popoolation2 .jar utilities)
#   - popoolation2  (pool-seq FST / sliding-window analyses)
#   - npstat        (pool-seq summary statistics, compiled from source)
#
# A lightweight micromamba environment (tools/env2) is created to supply
# samtools, bwa, OpenJDK, GSL (to compile npstat), and Perl + cpanm
# (needed by popoolation2).
#
# After running this script, source tools/env_practical2.sh (or restart your
# shell / RStudio session) to add the wrappers to your PATH.
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# This script lives in  <repo>/setup/ – derive the repo root from that.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TOOLS_DIR="$ROOT_DIR/tools"
BIN_DIR="$TOOLS_DIR/bin"
MAMBA_DIR="$TOOLS_DIR/.micromamba"
MAMBA_BIN="$MAMBA_DIR/bin/micromamba"
ENV_PREFIX="$TOOLS_DIR/env2"
PERL5_DIR="$TOOLS_DIR/perl5"
POPOOLATION2_URL="https://sourceforge.net/projects/popoolation2/files/popoolation2_1201.zip/download"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
log() { printf '\n[setup_practical2] %s\n' "$*"; }

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    printf 'ERROR: required command not found: %s\n' "$1" >&2
    exit 1
  fi
}

# ---------------------------------------------------------------------------
# 1. Micromamba (local, no sudo required)
# ---------------------------------------------------------------------------
install_micromamba() {
  if [[ -x "$MAMBA_BIN" ]]; then
    log "micromamba already present – skipping download"
    return
  fi
  require_cmd curl
  require_cmd tar
  log "Downloading micromamba into $MAMBA_DIR"
  mkdir -p "$MAMBA_DIR/bin"
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xj -C "$MAMBA_DIR/bin" --strip-components=1 bin/micromamba
}

# ---------------------------------------------------------------------------
# 2. Conda environment – samtools, bwa, OpenJDK, GSL + Perl
#    (samtools/bwa for alignment; java for popoolation2 jars;
#     GSL to compile npstat; Perl + cpanm for popoolation2 scripts)
# ---------------------------------------------------------------------------
ensure_env2() {
  install_micromamba
  if [[ -d "$ENV_PREFIX" ]]; then
    log "env2 already exists – skipping creation"
    return
  fi
  log "Creating conda environment env2 in $ENV_PREFIX"
  "$MAMBA_BIN" create -y -p "$ENV_PREFIX" \
    -c conda-forge -c bioconda \
    samtools bwa openjdk gsl perl perl-app-cpanminus
}

# ---------------------------------------------------------------------------
# 3. Perl module required by popoolation2
# ---------------------------------------------------------------------------
install_perl_module() {
  local marker="$PERL5_DIR/.installed"
  if [[ -f "$marker" ]]; then
    log "Perl module already installed – skipping"
    return
  fi
  log "Installing Perl module Text::NSP::Measures::2D::Fisher::twotailed"
  mkdir -p "$PERL5_DIR"
  export PATH="$ENV_PREFIX/bin:$PATH"
  export PERL5LIB="$PERL5_DIR/lib/perl5${PERL5LIB:+:$PERL5LIB}"
  "$ENV_PREFIX/bin/cpanm" --quiet -L "$PERL5_DIR" \
    Text::NSP::Measures::2D::Fisher::twotailed
  touch "$marker"
}

# ---------------------------------------------------------------------------
# 4. popoolation2
# ---------------------------------------------------------------------------
install_popoolation2() {
  local dest="$TOOLS_DIR/popoolation2"
  if [[ -d "$dest/popoolation2_1201" ]]; then
    log "popoolation2 already installed – skipping"
    return
  fi
  require_cmd curl
  require_cmd unzip
  log "Downloading popoolation2"
  mkdir -p "$dest"
  curl -fsSL "$POPOOLATION2_URL" -o "$TOOLS_DIR/popoolation2.zip"
  unzip -q "$TOOLS_DIR/popoolation2.zip" -d "$dest"
  find "$dest" -type f \( -name '*.pl' -o -name '*.jar' \) -exec chmod 0755 {} +
  rm -f "$TOOLS_DIR/popoolation2.zip"
}

# ---------------------------------------------------------------------------
# 5. npstat  (clone + compile)
# ---------------------------------------------------------------------------
install_npstat() {
  local dest="$TOOLS_DIR/npstat"
  if [[ ! -d "$dest/.git" ]]; then
    require_cmd git
    log "Cloning npstat"
    git clone --depth 1 https://github.com/lucaferretti/npstat.git "$dest"
  fi
  if [[ -x "$dest/npstat" ]]; then
    log "npstat binary already compiled – skipping"
    return
  fi
  require_cmd gcc
  log "Compiling npstat against GSL from env2"
  gcc -o "$dest/npstat" "$dest/NPStat-v1.c" \
    -I"$ENV_PREFIX/include" \
    -L"$ENV_PREFIX/lib" \
    -Wl,-rpath,"$ENV_PREFIX/lib" \
    -lgsl -lgslcblas -lm
  chmod 0755 "$dest/npstat"
}

# ---------------------------------------------------------------------------
# 6. Wrapper scripts in tools/bin
# ---------------------------------------------------------------------------
write_wrapper() {
  local target="$BIN_DIR/$1"
  cat >"$target" <<WRAPPER
#!/usr/bin/env bash
set -euo pipefail
tool_dir="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")/.." && pwd)"
${2}
WRAPPER
  chmod 0755 "$target"
}

write_wrappers_practical2() {
  mkdir -p "$BIN_DIR"

  # npstat – standalone binary, copy directly
  cp "$TOOLS_DIR/npstat/npstat" "$BIN_DIR/npstat"
  chmod 0755 "$BIN_DIR/npstat"

  # popoolation2 Perl scripts – one wrapper per script
  local pop_dir="$TOOLS_DIR/popoolation2/popoolation2_1201"
  local script_path
  for script_path in \
    "$pop_dir"/*.pl \
    "$pop_dir"/export/*.pl \
    "$pop_dir"/indel_filtering/*.pl; do
    [[ -f "$script_path" ]] || continue
    local name
    name="$(basename "$script_path")"
    write_wrapper "$name" \
      "export PERL5LIB=\"\$tool_dir/perl5/lib/perl5\${PERL5LIB:+:\$PERL5LIB}\"
export PATH=\"\$tool_dir/env2/bin:\$PATH\"
exec \"$script_path\" \"\$@\""
  done

  log "Wrappers written to $BIN_DIR"
}

# ---------------------------------------------------------------------------
# 7. env_practical2.sh  (source this to activate the environment)
# ---------------------------------------------------------------------------
write_env_script() {
  cat >"$TOOLS_DIR/env_practical2.sh" <<'ENVSCRIPT'
#!/usr/bin/env bash
# Source this file to add Practical 2 tools to the current shell session:
#   source tools/env_practical2.sh
_tool_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="$_tool_dir/bin:$_tool_dir/env2/bin:$PATH"
export LD_LIBRARY_PATH="$_tool_dir/env2/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PERL5LIB="$_tool_dir/perl5/lib/perl5${PERL5LIB:+:$PERL5LIB}"
unset _tool_dir
ENVSCRIPT
  chmod 0755 "$TOOLS_DIR/env_practical2.sh"
}

# ---------------------------------------------------------------------------
# 8. Persistent PATH hooks  (~/.bashrc, ~/.profile, ~/.Renviron)
# ---------------------------------------------------------------------------
upsert_block() {
  local target="$1" marker="$2" content="$3"
  local start="# >>> $marker >>>" end="# <<< $marker <<<"
  mkdir -p "$(dirname "$target")"
  touch "$target"
  if grep -Fq "$start" "$target"; then
    awk -v s="$start" -v e="$end" -v r="$content" '
      BEGIN{i=0} $0==s{print s; print r; print e; i=1; next}
      $0==e{i=0; next} !i{print}
    ' "$target" >"$target.tmp" && mv "$target.tmp" "$target"
  else
    { printf '\n%s\n%s\n%s\n' "$start" "$content" "$end"; } >>"$target"
  fi
}

install_startup_hooks() {
  local shell_block
  shell_block="$(printf 'if [ -f "%s/env_practical2.sh" ]; then\n  . "%s/env_practical2.sh"\nfi' \
    "$TOOLS_DIR" "$TOOLS_DIR")"
  local renviron_block
  renviron_block="$(printf 'PATH=%s:%s/bin:%s\nLD_LIBRARY_PATH=%s/lib:%s\nPERL5LIB=%s/lib/perl5:%s' \
    "$BIN_DIR" "$ENV_PREFIX" "${PATH:-}" \
    "$ENV_PREFIX" "${LD_LIBRARY_PATH:-}" \
    "$PERL5_DIR" "${PERL5LIB:-}")"

  upsert_block "$HOME/.bashrc"    "EVO_PRACTICAL2_TOOLS" "$shell_block"
  upsert_block "$HOME/.profile"   "EVO_PRACTICAL2_TOOLS" "$shell_block"
  upsert_block "$HOME/.Renviron"  "EVO_PRACTICAL2_TOOLS" "$renviron_block"
  log "Startup hooks written to ~/.bashrc, ~/.profile, ~/.Renviron"
}

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
main() {
  require_cmd bash
  require_cmd curl
  require_cmd git
  require_cmd gcc
  require_cmd unzip

  mkdir -p "$TOOLS_DIR"

  ensure_env2
  install_perl_module
  install_popoolation2
  install_npstat
  write_wrappers_practical2
  write_env_script
  install_startup_hooks

  log "Done. Tools installed under $TOOLS_DIR"
  log "  samtools      →  $ENV_PREFIX/bin/samtools  (on PATH via env_practical2.sh)"
  log "  bwa           →  $ENV_PREFIX/bin/bwa        (on PATH via env_practical2.sh)"
  log "  java          →  $ENV_PREFIX/bin/java        (on PATH via env_practical2.sh)"
  log "  popoolation2  →  $(ls "$TOOLS_DIR/popoolation2/popoolation2_1201/"*.pl 2>/dev/null | wc -l) .pl scripts wrapped in $BIN_DIR"
  log "  npstat        →  $BIN_DIR/npstat"
  log ""
  log "To activate in the current shell, run:"
  log "  source $TOOLS_DIR/env_practical2.sh"
}

main "$@"