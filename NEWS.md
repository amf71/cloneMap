# cloneMap 1.0.1

## Installation

* `rgeos` is no longer required, and the manual installation steps
  previously described in the README (installing `geos` via Homebrew and
  `rgeos` from R-Forge) are no longer needed. All dependencies are now on
  CRAN and install automatically with
  `remotes::install_github("amf71/cloneMap")`.

  `rgeos` was never called directly by cloneMap. It was an indirect
  backend for `raster::rasterToPolygons(dissolve = TRUE)`, and since
  raster 3.5-29 (2022-08-14) that operation has been handled by `terra`
  instead, so the dependency had become stale. cloneMap now declares
  `raster (>= 3.5-29)` to make the requirement explicit.

* `qlcMatrix` is no longer a dependency. It was used in a single place,
  for `rowMax()` on a two-column matrix, which base `pmax()` does
  identically.

Plot output is unchanged by both of the above.

## Bug fixes

* Fixed misaligned indexing when selecting nuclei for daughter clones.
  The candidate nuclei were filtered to the most widely spaced fifth, but
  the corresponding vector of inter-nucleus distances was not filtered
  alongside them, so the distance passed to the clone growth step could
  belong to a different candidate than the one chosen.

* `cloneMap()` no longer aborts when a randomly chosen nucleus falls
  outside its parent clone. Clone positions are semi-random, so this
  happened intermittently (around 3% of runs for unrooted data plotted
  with `tissue_border = TRUE`). A fresh set of nuclei is now drawn
  instead, and an error is raised only if the problem persists across
  repeated attempts.

* Fixed a typo in the `cloneMap()` examples that referred to
  `tree_examkeple_poly` rather than `tree_example_poly`.

## Documentation

* Regenerated help files, which had drifted out of sync with the code:
  `inc_parents_ccf` was undocumented for both `cloneMap()` and
  `calc_clonal_diversity()`, `clone_names` was undocumented for
  `make.CCFs.tree.consistant()`, and the documented default for
  `smoothing.par` was stale.
