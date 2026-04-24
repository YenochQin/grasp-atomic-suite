#!/usr/bin/env sh

set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)

DRY_RUN=0

usage() {
  cat <<'EOF'
Usage: ./scripts/clean-build-artifacts.sh [--dry-run]

Clean build outputs for both workflows:
  - run `make cleanall` when available
  - remove top-level `build/` and `build-debug/`
  - remove non-hidden files under `bin/` and `lib/`

Options:
  --dry-run   Show what would be removed without changing anything
  -h, --help  Show this help message
EOF
}

run_cmd() {
  printf '+'
  for arg in "$@"; do
    printf ' %s' "$arg"
  done
  printf '\n'

  if [ "$DRY_RUN" -eq 0 ]; then
    "$@"
  fi
}

clean_dir_contents() {
  dir_path=$1

  if [ ! -d "$dir_path" ]; then
    return 0
  fi

  if find "$dir_path" -mindepth 1 -maxdepth 1 ! -name '.*' | grep -q .; then
    run_cmd find "$dir_path" -mindepth 1 -maxdepth 1 ! -name '.*' -exec rm -rf {} +
  else
    printf 'No removable entries under %s\n' "$dir_path"
  fi
}

remove_dir_if_exists() {
  dir_path=$1

  if [ -d "$dir_path" ]; then
    run_cmd rm -rf "$dir_path"
  else
    printf 'Directory not present: %s\n' "$dir_path"
  fi
}

while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run)
      DRY_RUN=1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'Unknown option: %s\n\n' "$1" >&2
      usage >&2
      exit 1
      ;;
  esac
  shift
done

printf 'Repository root: %s\n' "$REPO_ROOT"

if [ -f "$REPO_ROOT/Makefile" ]; then
  if [ "$DRY_RUN" -eq 0 ]; then
    printf '+ make -C %s cleanall\n' "$REPO_ROOT"
    make -C "$REPO_ROOT" cleanall || printf 'Warning: `make cleanall` failed; continuing with directory cleanup.\n' >&2
  else
    printf '+ make -C %s cleanall\n' "$REPO_ROOT"
  fi
fi

clean_dir_contents "$REPO_ROOT/bin"
clean_dir_contents "$REPO_ROOT/lib"
remove_dir_if_exists "$REPO_ROOT/build"
remove_dir_if_exists "$REPO_ROOT/build-debug"
