# MRDS duplicate matcher

A Shiny app that identifies duplicate detections between two tandem observers
on an aerial survey and builds the detection history (11 / 10 / 01) needed for
mark-recapture distance sampling. It is exported with
[shinylive](https://posit-dev.github.io/r-shinylive/) and runs entirely in the
browser as WebAssembly: there is no server, and uploaded data stay on the
user's computer.

## What it does

- Reads a CSV with one row per detection and lets you map your columns to
  observer, time, longitude, latitude, group size and, optionally, detection
  ID, transect, side, species, perpendicular distance and sighting angle.
- Scores candidate pairs with a Fellegi-Sunter model fitted by EM, the
  discrepancy mixture model of Hamilton et al. (2018), or fixed tolerances.
- Assigns duplicates by Hungarian one-to-one matching, Hamilton et al.'s
  stochastic cluster resolution, or a Metropolis-Hastings sample of matchings.
- Shows diagnostics and downloads the matched pairs, all draws, the candidate
  pairs and the detection history as CSV.

Abundance estimation is not part of the app.

## Files

| File | Contents |
|---|---|
| `app.R` | User interface and server |
| `R/01_mrds_duplicate_match.R` | Fellegi-Sunter matcher, assignment, detection history, diagnostics |
| `R/02_hamilton_matching.R` | Hamilton et al. scoring and cluster resolution |
| `R/03_app_helpers.R` | Example data, column mapping, time parsing, run wrapper |
| `tests/test_matching.R` | Checks; run with `Rscript tests/test_matching.R` |
| `.github/workflows/deploy-app.yaml` | Builds and deploys to GitHub Pages |

`R/01_*.R` is an unchanged copy of `scripts/mrds_duplicate_match.R` from the
analysis project, and `R/02_*.R` is sections 0-3 and 7 of
`scripts/hamilton_options.R`. Copy them across again when those scripts change.

## Run locally

```r
shiny::runApp()
```

Needs shiny, bslib, DT, dplyr, tidyr, purrr, tibble, geosphere and clue.

## Preview the WebAssembly build locally

```r
install.packages("shinylive")
shinylive::export(".", "site")
httpuv::runStaticServer("site")
```

## Deploy to GitHub Pages

1. Create a public repository on GitHub and push this folder to its `main`
   branch.
2. In the repository, open Settings > Pages and set Source to "GitHub Actions".
3. The workflow runs on every push to `main`. When it finishes the app is at
   `https://<user>.github.io/<repository>/`.

The first page load downloads R and the packages (tens of MB); later loads use
the browser cache.

## Limits

- Slower than desktop R and limited by browser memory. Jobs needing more than
  5 million within-transect A x B comparisons are refused; map a transect
  column to keep comparisons manageable.
- Coordinates must be decimal degrees (WGS84).

## Reference

Hamilton, O. N. P., Kincaid, S. E., Constantine, R., Kozmian-Ledward, L.,
Walker, C. G. & Fewster, R. M. (2018) Accounting for uncertainty in duplicate
identification and group size judgments in mark-recapture distance sampling.
*Methods in Ecology and Evolution* 9: 354-362.
