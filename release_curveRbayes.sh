#!/usr/bin/env bash
# =============================================================================
# release_curveRbayes.sh — regenerate docs, then cut a versioned, tagged release
# of the curveRbayes R package.
#
# Usage (from the curveRbayes repo root, in Git Bash / WSL / macOS / Linux):
#   ./release_curveRbayes.sh                    # version 0.4.0, tag v0.4.0
#   ./release_curveRbayes.sh 0.4.0 v0.4.0       # explicit version + tag
#   ./release_curveRbayes.sh 0.4.1              # version 0.4.1, tag v0.4.1
#
# Env toggles:
#   RSCRIPT=/path/to/Rscript.exe   # if Rscript isn't on PATH (common on Windows)
#   NO_DOCS=1                      # skip doc regeneration (docs already current)
#   NO_SITE=1                      # run devtools::document() but skip pkgdown
#   NO_NEWS=1                      # skip the NEWS.md edit
#   NO_PUSH=1                      # do everything locally, don't push / release
#
# Order:
#   1. Preflight (repo, package == curveRbayes, refactor present, tag free).
#   2. Bump DESCRIPTION Version.
#   3. Prepend NEWS.md section.
#   4. Regenerate docs: devtools::document() [required] + pkgdown::build_site()
#      [optional, non-fatal] — so the committed docs match this version.
#   5. Show diff + PAUSE for confirmation.
#   6. Commit, annotated tag, push.
#   7. GitHub Release via `gh` if available, else print the URL.
#
# Never force-pushes, never deletes tags. Safe to read before running.
# =============================================================================
set -euo pipefail

VERSION="${1:-0.4.0}"
TAG="${2:-v${VERSION}}"
PKG_EXPECTED="curveRbayes"

# ── 1. Preflight ─────────────────────────────────────────────────────────────
command -v git >/dev/null || { echo "ERROR: git not found on PATH."; exit 1; }
git rev-parse --is-inside-work-tree >/dev/null 2>&1 \
  || { echo "ERROR: not inside a git repository. cd to the curveRbayes repo first."; exit 1; }
cd "$(git rev-parse --show-toplevel)"
[ -f DESCRIPTION ] || { echo "ERROR: no DESCRIPTION here — is this the package root?"; exit 1; }

PKG_NAME="$(awk -F': *' '/^Package:/{print $2; exit}' DESCRIPTION | tr -d '\r')"
if [ "$PKG_NAME" != "$PKG_EXPECTED" ]; then
  echo "WARNING: DESCRIPTION Package is '$PKG_NAME', expected '$PKG_EXPECTED'."
  read -r -p "Continue anyway? [y/N] " ans; [ "${ans:-N}" = "y" ] || exit 1
fi

# Guard: make sure the refactored source is actually in the tree.
if [ -f R/fit_calibration_bayes.R ] && ! grep -q "include_measurement_error" R/fit_calibration_bayes.R; then
  echo "WARNING: R/fit_calibration_bayes.R has no 'include_measurement_error'."
  echo "         Did you copy in the refactored predict_bayes.R / fit_calibration_bayes.R?"
  read -r -p "Continue anyway? [y/N] " ans; [ "${ans:-N}" = "y" ] || exit 1
fi

BRANCH="$(git rev-parse --abbrev-ref HEAD)"
echo "Repo:    $(git rev-parse --show-toplevel)"
echo "Package: $PKG_NAME"
echo "Branch:  $BRANCH"
echo "Version: $VERSION   Tag: $TAG"
echo

if git rev-parse -q --verify "refs/tags/$TAG" >/dev/null; then
  echo "ERROR: tag $TAG already exists locally. Aborting."; exit 1
fi
if git ls-remote --exit-code --tags origin "$TAG" >/dev/null 2>&1; then
  echo "ERROR: tag $TAG already exists on origin. Aborting."; exit 1
fi

# Resolve Rscript (needed for doc regeneration). PATH first, then $RSCRIPT, then
# common Windows install locations.
resolve_rscript() {
  if [ -n "${RSCRIPT:-}" ]; then echo "$RSCRIPT"; return; fi
  if command -v Rscript >/dev/null 2>&1; then echo "Rscript"; return; fi
  local cand
  for cand in "/c/Program Files/R"/R-*/bin/x64/Rscript.exe \
              "/c/Program Files/R"/R-*/bin/Rscript.exe; do
    [ -x "$cand" ] && { echo "$cand"; return; }
  done
  echo ""
}
RS="$(resolve_rscript)"

# ── 2. Bump DESCRIPTION Version ──────────────────────────────────────────────
if [ -n "$RS" ] && \
   "$RS" -e 'quit(status = as.integer(!requireNamespace("desc", quietly=TRUE)))' >/dev/null 2>&1; then
  "$RS" -e "desc::desc_set_version('${VERSION}'); cat('DESCRIPTION Version ->', as.character(desc::desc_get_version()), '\n')"
else
  echo "(desc/Rscript not available — editing DESCRIPTION with awk)"
  awk -v v="$VERSION" 'BEGIN{done=0}
    /^Version:/ && !done {print "Version: " v; done=1; next} {print}' \
    DESCRIPTION > DESCRIPTION.tmp && mv DESCRIPTION.tmp DESCRIPTION
fi
DESC_VER="$(awk -F': *' '/^Version:/{print $2; exit}' DESCRIPTION | tr -d '\r')"
[ "$DESC_VER" = "$VERSION" ] || { echo "ERROR: DESCRIPTION Version is '$DESC_VER', expected '$VERSION'."; exit 1; }

# ── 3. NEWS.md ───────────────────────────────────────────────────────────────
NEWS_BODY="$(cat <<'EOF'
* **Unified the precision-profile variance definition** across the grid and the
  per-sample back-calculation via a new `include_measurement_error` argument to
  `fit_calibration_bayes()` (default `TRUE`), threaded identically into
  `predict_grid_bayes()` and `predict_samples_bayes()`. This closes the gap where
  Bayesian sample points floated ~9.6 %CV below the precision profile with a
  concentration-dependent slope: with the default they now lie on the
  measurement/CDAN profile.
* `predict_grid_bayes()` and `predict_samples_bayes()` now share a single
  `.obs_noise_sigma()` noise definition, so the grid and sample paths cannot
  diverge again.
* `include_measurement_error = FALSE` yields a curve/parameter-uncertainty-only
  profile; in that mode the grid inverts a fixed posterior-mean reference response
  across draws (non-degenerate) to match the sample path.
* The active mode is recorded in the `noise_mode` column and in the result `meta`.
EOF
)"

if [ "${NO_NEWS:-0}" != "1" ]; then
  DATE="$(date +%Y-%m-%d)"
  HEADER="# ${PKG_NAME} ${VERSION} (${DATE})"
  if [ -f NEWS.md ] && grep -q "^# ${PKG_NAME} ${VERSION}\b" NEWS.md; then
    echo "(NEWS.md already has a ${VERSION} section — leaving it as is)"
  else
    { printf '%s\n\n%s\n\n' "$HEADER" "$NEWS_BODY"; [ -f NEWS.md ] && cat NEWS.md; } > NEWS.md.tmp \
      && mv NEWS.md.tmp NEWS.md
    echo "Prepended ${VERSION} section to NEWS.md"
  fi
fi

# ── 4. Regenerate documentation ──────────────────────────────────────────────
if [ "${NO_DOCS:-0}" = "1" ]; then
  echo "NO_DOCS=1 — skipping doc regeneration."
elif [ -z "$RS" ]; then
  echo "WARNING: Rscript not found; cannot regenerate docs automatically."
  echo "  Run  devtools::document()  (and pkgdown::build_site()) in R, then re-run"
  echo "  with NO_DOCS=1, or set RSCRIPT=/path/to/Rscript.exe and re-run."
  read -r -p "Continue without regenerating docs? [y/N] " ans; [ "${ans:-N}" = "y" ] || exit 1
else
  echo "Regenerating roxygen docs (NAMESPACE + man/) ..."
  "$RS" -e 'devtools::document()'            # required: updates NAMESPACE + man/*.Rd
  if [ "${NO_SITE:-0}" = "1" ]; then
    echo "NO_SITE=1 — skipping pkgdown site build."
  else
    echo "Building pkgdown site (non-fatal) ..."
    if ! "$RS" -e 'pkgdown::build_site(preview = FALSE)'; then
      echo "WARNING: pkgdown::build_site() failed (docs site only; the package is unaffected)."
      echo "  Common cause: a newly exported topic missing from _pkgdown.yml's reference:"
      echo "  index — add it there (or @keywords internal to hide it), then rebuild."
      read -r -p "Continue the release without a rebuilt site? [y/N] " ans; [ "${ans:-N}" = "y" ] || exit 1
    fi
  fi
fi

# ── 5. Review & confirm ──────────────────────────────────────────────────────
git add -A
echo
echo "──────── staged for release ────────"
git --no-pager status --short
echo "────────────────────────────────────"
git --no-pager diff --cached -- DESCRIPTION NEWS.md NAMESPACE || true
echo
read -r -p "Commit, tag ${TAG}, and push? [y/N] " ans
[ "${ans:-N}" = "y" ] || { echo "Aborted before commit. (Edits are staged but uncommitted.)"; exit 1; }

# ── 6. Commit + annotated tag ────────────────────────────────────────────────
git commit -m "Release ${TAG}: unify grid/sample precision via include_measurement_error; regenerate docs"
git tag -a "$TAG" -m "${PKG_NAME} ${VERSION}

${NEWS_BODY}"
echo "Created annotated tag ${TAG}."

# ── 7. Push ──────────────────────────────────────────────────────────────────
if [ "${NO_PUSH:-0}" = "1" ]; then
  echo "NO_PUSH=1 set — skipping push. To push later:"
  echo "  git push origin ${BRANCH} && git push origin ${TAG}"
  exit 0
fi
git push origin "$BRANCH"
git push origin "$TAG"
echo "Pushed ${BRANCH} and ${TAG} to origin."

# ── 8. GitHub Release (optional) ─────────────────────────────────────────────
if command -v gh >/dev/null && gh auth status >/dev/null 2>&1; then
  printf '%s\n' "$NEWS_BODY" > .release_notes.tmp
  gh release create "$TAG" \
    --title "${PKG_NAME} ${VERSION}" \
    --notes-file .release_notes.tmp \
    --target "$BRANCH"
  rm -f .release_notes.tmp
  echo "GitHub Release ${TAG} created."
else
  REMOTE_URL="$(git remote get-url origin 2>/dev/null || echo '')"
  echo "gh CLI not available/authenticated — the tag is pushed; create the Release here:"
  case "$REMOTE_URL" in
    git@github.com:*)     SLUG="${REMOTE_URL#git@github.com:}"; SLUG="${SLUG%.git}";;
    https://github.com/*) SLUG="${REMOTE_URL#https://github.com/}"; SLUG="${SLUG%.git}";;
    *) SLUG="";;
  esac
  [ -n "$SLUG" ] && echo "  https://github.com/${SLUG}/releases/new?tag=${TAG}"
fi
