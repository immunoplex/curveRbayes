# ============================================================
# Build and deploy the curveRbayes pkgdown site
# ============================================================
#
# Run from the curveRbayes project root (the folder with DESCRIPTION).
# Assumes: devtools::check() passes, git is set up, remote is origin.


# ── Step 1: Place the pkgdown config files ──────────────────
#
# You should already have _pkgdown.yml at the project root.
# If you're replacing it with the new version, overwrite it now.
#
# The extra.css goes in pkgdown/extra.css (pkgdown picks it up
# automatically from that path — no config entry needed).
#
# From the R console:

dir.create("pkgdown", showWarnings = FALSE)
file.copy("path/to/extra.css", "pkgdown/extra.css", overwrite = TRUE)

# Or from the terminal:
#   mkdir -p pkgdown
#   cp extra.css pkgdown/extra.css


# ── Step 2: Add a logo (optional but recommended) ──────────
#
# The _pkgdown.yml references man/figures/logo_small.png.
# If you don't have a logo yet, create the directory and drop
# one in:
#
#   dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
#
# A 240×278 px PNG works well. If you skip this step, remove
# the logo: block from _pkgdown.yml to avoid a build warning.


# ── Step 3: Build the site ──────────────────────────────────
#
# This generates the full site into docs/.

pkgdown::build_site(dest_dir = "docs")

# If you only want to rebuild specific parts:
#   pkgdown::build_home()
#   pkgdown::build_reference()
#   pkgdown::build_articles()
#   pkgdown::build_news()
#
# Expect the vignette build to take ~10 min (Stan compilation).


# ── Step 4: Verify locally ─────────────────────────────────
#
# Open docs/index.html in a browser and check:
#   - Home page renders with description and README content
#   - Reference page shows all 7 sections with correct grouping
#   - Articles tab links to the bayesian-quickstart vignette
#   - News tab shows the 0.2.0 changelog
#   - Navbar logo and GitHub link work
#   - Search works (type a function name)


# ── Step 5: Commit and push ────────────────────────────────
#
# From the terminal (Git Bash or RStudio terminal):

# git add _pkgdown.yml pkgdown/ docs/
# git commit -m "Add pkgdown site"
# git push origin main


# ── Step 6: Enable GitHub Pages ─────────────────────────────
#
# Go to: https://github.com/immunoplex/curveRbayes/settings/pages
#
#   Source:  "Deploy from a branch"
#   Branch:  main
#   Folder:  /docs
#   → Save
#
# Wait 1-2 minutes, then visit:
#   https://immunoplex.github.io/curveRbayes/
#
# The site should match the curveRcore layout exactly.


# ── Step 7: Confirm cross-links ────────────────────────────
#
# Check that the ecosystem links work across all three sites:
#   https://immunoplex.github.io/curveRcore/
#   https://immunoplex.github.io/curveRfreq/
#   https://immunoplex.github.io/curveRbayes/


# ── Troubleshooting ────────────────────────────────────────
#
# | Problem                          | Fix                                          |
# |----------------------------------|----------------------------------------------|
# | "Topic not found" in reference   | Run devtools::document() then rebuild         |
# | Vignette missing from Articles   | Check VignetteBuilder: knitr in DESCRIPTION   |
# | Logo not showing                 | Verify man/figures/logo_small.png exists       |
# | 404 on GitHub Pages              | Settings → Pages: branch main, folder /docs   |
# | CSS not applied                  | Confirm pkgdown/extra.css exists, rebuild      |
# | MathJax not rendering            | The flatly theme + BS5 includes MathJax;       |
# |                                  | if equations break, add under template:        |
# |                                  |   math-rendering:                              |
# |                                  |     mathjax: true                              |
# | Search not working               | Requires docs/ served via HTTPS (GitHub Pages) |
