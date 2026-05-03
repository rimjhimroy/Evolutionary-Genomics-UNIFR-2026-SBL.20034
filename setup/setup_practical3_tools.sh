#!/usr/bin/env bash
# =============================================================================
# setup_practical3_tools.sh
#
# Installs the command-line tools needed for Practical 3:
#   - structure          (Bayesian clustering via Bioconda)
#   - CLUMPP             (cluster alignment, standalone binary)
#   - structureHarvester (post-processing script, Python)
#
# A lightweight micromamba environment (tools/env3) is created to host
# structure (from Bioconda) and a Python interpreter for structureHarvester.
# CLUMPP is downloaded as a pre-compiled Linux binary.
#
# After running this script, source tools/env_practical3.sh (or restart your
# shell / RStudio session) to add the wrappers to your PATH.
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TOOLS_DIR="$ROOT_DIR/tools"
BIN_DIR="$TOOLS_DIR/bin"
MAMBA_DIR="$TOOLS_DIR/.micromamba"
MAMBA_BIN="$MAMBA_DIR/bin/micromamba"
ENV_PREFIX="$TOOLS_DIR/env3"
CLUMPP_URL="https://rosenberglab.stanford.edu/software/CLUMPP_Linux64.1.1.2.tar.gz"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
log() { printf '\n[setup_practical3] %s\n' "$*"; }

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
# 2. Conda environment – structure + python (for structureHarvester)
# ---------------------------------------------------------------------------
ensure_env3() {
  install_micromamba
  if [[ -d "$ENV_PREFIX" ]]; then
    log "env3 already exists – skipping creation"
    return
  fi
  log "Creating conda environment env3 in $ENV_PREFIX"
  "$MAMBA_BIN" create -y -p "$ENV_PREFIX" \
    -c conda-forge -c bioconda \
    structure python
}

# ---------------------------------------------------------------------------
# 3. CLUMPP  (pre-compiled Linux binary)
# ---------------------------------------------------------------------------
install_clumpp() {
  local dest="$TOOLS_DIR/clumpp"
  if [[ -x "$dest/CLUMPP_Linux64.1.1.2/CLUMPP" ]]; then
    log "CLUMPP already installed – skipping"
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

# ---------------------------------------------------------------------------
# 4. structureHarvester  (Python script, cloned from GitHub)
# ---------------------------------------------------------------------------
install_structure_harvester() {
  local dest="$TOOLS_DIR/structureHarvester"
  if [[ -d "$dest/.git" ]]; then
    log "structureHarvester already cloned – skipping"
    return
  fi
  require_cmd git
  log "Cloning structureHarvester"
  git clone --depth 1 https://github.com/dentearl/structureHarvester.git "$dest"
  chmod 0755 "$dest/structureHarvester.py"
}

# ---------------------------------------------------------------------------
# 5. Wrapper scripts in tools/bin
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

write_wrappers_practical3() {
  mkdir -p "$BIN_DIR"

  # structure – run the binary from env3
  write_wrapper "structure" \
    "exec \"\$tool_dir/env3/bin/structure\" \"\$@\""

  # CLUMPP – standalone binary
  cp "$TOOLS_DIR/clumpp/CLUMPP_Linux64.1.1.2/CLUMPP" "$BIN_DIR/CLUMPP"
  chmod 0755 "$BIN_DIR/CLUMPP"

  # structureHarvester – Python script via env3 interpreter
  write_wrapper "structureHarvester" \
    "exec \"\$tool_dir/env3/bin/python\" \"\$tool_dir/structureHarvester/structureHarvester.py\" \"\$@\""
  # also expose as .py alias
  cp "$BIN_DIR/structureHarvester" "$BIN_DIR/structureHarvester.py"
  chmod 0755 "$BIN_DIR/structureHarvester.py"

  log "Wrappers written to $BIN_DIR"
}

# ---------------------------------------------------------------------------
# 6. env_practical3.sh  (source this to activate the environment)
# ---------------------------------------------------------------------------
write_env_script() {
  cat >"$TOOLS_DIR/env_practical3.sh" <<'ENVSCRIPT'
#!/usr/bin/env bash
# Source this file to add Practical 3 tools to the current shell session:
#   source tools/env_practical3.sh
_tool_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="$_tool_dir/bin:$_tool_dir/env3/bin:$PATH"
export LD_LIBRARY_PATH="$_tool_dir/env3/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
unset _tool_dir
ENVSCRIPT
  chmod 0755 "$TOOLS_DIR/env_practical3.sh"
}

# ---------------------------------------------------------------------------
# 7. Persistent PATH hooks  (~/.bashrc, ~/.profile, ~/.Renviron)
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
  shell_block="$(printf 'if [ -f "%s/env_practical3.sh" ]; then\n  . "%s/env_practical3.sh"\nfi' \
    "$TOOLS_DIR" "$TOOLS_DIR")"
  local renviron_block
  renviron_block="$(printf 'PATH=%s:%s/bin:%s\nLD_LIBRARY_PATH=%s/lib:%s' \
    "$BIN_DIR" "$ENV_PREFIX" "${PATH:-}" \
    "$ENV_PREFIX" "${LD_LIBRARY_PATH:-}")"

  upsert_block "$HOME/.bashrc"   "EVO_PRACTICAL3_TOOLS" "$shell_block"
  upsert_block "$HOME/.profile"  "EVO_PRACTICAL3_TOOLS" "$shell_block"
  upsert_block "$HOME/.Renviron" "EVO_PRACTICAL3_TOOLS" "$renviron_block"
  log "Startup hooks written to ~/.bashrc, ~/.profile, ~/.Renviron"
}

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
main() {
  require_cmd bash
  require_cmd curl
  require_cmd tar
  require_cmd git

  mkdir -p "$TOOLS_DIR"

  ensure_env3
  install_clumpp
  install_structure_harvester
  write_wrappers_practical3
  write_env_script
  install_startup_hooks

  log "Done. Tools installed under $TOOLS_DIR"
  log "  structure           →  $BIN_DIR/structure  (env3: $ENV_PREFIX/bin/structure)"
  log "  CLUMPP              →  $BIN_DIR/CLUMPP"
  log "  structureHarvester  →  $BIN_DIR/structureHarvester"
  log ""
  log "To activate in the current shell, run:"
  log "  source $TOOLS_DIR/env_practical3.sh"
}

main "$@"