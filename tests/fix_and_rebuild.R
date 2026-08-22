# ============================================================
# Fix curveRbayes pkgdown logo and rebuild the site
# ============================================================
#
# Run from the curveRbayes project root (folder with DESCRIPTION).
# The problem is almost certainly one of three things:
#   1. man/figures/logo_small.png is not committed to git
#   2. pkgdown also needs man/figures/logo.png for the home-page badge
#   3. Favicons are stale (built before the logo existed)
#
# This script fixes all three, then rebuilds and pushes.


# ── Step 1: Place the logo files ────────────────────────────
#
# pkgdown needs TWO copies of the logo in man/figures/:
#
#   logo_small.png  — referenced in _pkgdown.yml navbar (already set)
#   logo.png        — auto-detected by pkgdown for the home-page badge
#                     and used as the source image for favicon generation
#
# Copy logo_small.png to both filenames:

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
file.copy("logo_small.png", "man/figures/logo_small.png", overwrite = TRUE)
file.copy("logo_small.png", "man/figures/logo.png",       overwrite = TRUE)

# Verify both files are present:
stopifnot(file.exists("man/figures/logo_small.png"))
stopifnot(file.exists("man/figures/logo.png"))
message("✔ Logo files placed in man/figures/")


# ── Step 2: Regenerate favicons ─────────────────────────────
#
# pkgdown generates favicons from man/figures/logo.png.
# Must be done before build_site() so the favicon files land in pkgdown/
# and are then copied to docs/.

pkgdown::build_favicons(overwrite = TRUE)
message("✔ Favicons rebuilt")


# ── Step 3: Rebuild the full site ───────────────────────────

pkgdown::build_site()
message("✔ Site built")


# ── Step 4: Verify the logo landed in docs/ ─────────────────
#
# pkgdown copies man/figures/ into docs/reference/figures/.
# The navbar image and home-page badge both reference the file from there.

logo_in_docs <- file.exists("docs/reference/figures/logo_small.png")
message(if (logo_in_docs) "✔ logo_small.png found in docs/" else
  "✖ logo_small.png NOT found in docs/ — check the build output above")


# ── Step 5: Commit and push ─────────────────────────────────
#
# From the terminal (Git Bash or RStudio terminal):
#
#   git add man/figures/logo_small.png
#   git add man/figures/logo.png
#   git add pkgdown/favicon/                  # new favicon files
#   git add docs/
#   git commit -m "Add logo and rebuild pkgdown site"
#   git push origin main
#
# IMPORTANT: man/figures/ must be committed — not just present locally.
# If git status shows the logo files as untracked, add them explicitly.
# Check that .gitignore does not contain man/ or *.png.


# ── Troubleshooting ─────────────────────────────────────────
#
# Logo appears in docs/ locally but not on GitHub Pages
#   → The man/figures/ files were not committed. Run git add man/figures/
#
# Logo shows in the navbar but not on the home page
#   → man/figures/logo.png is missing. Step 1 above creates it.
#
# Favicon is the generic pkgdown icon, not the hex sticker
#   → Rerun build_favicons(overwrite = TRUE) then build_site()
#
# Navbar shows a broken image icon
#   → The image: path in _pkgdown.yml does not match the filename.
#     Confirm it reads:  image: man/figures/logo_small.png
#
# build_favicons() errors with "no logo found"
#   → man/figures/logo.png must exist before running build_favicons()
