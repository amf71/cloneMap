# cloneMap 1.0.2

## Bug fixes

* Fixed column-major index decoding in the internal
  `matrix.index.to.coordinates()`. R matrices are stored column-major, but the
  helper decoded a matrix index as `x = index %% nrow` and
  `y = floor(index / ncol) + 1`. That form is correct for most indices but
  wrong for every index that is an exact multiple of `nrow`, i.e. the last row
  of each column: it returned row 0, which is not a valid index, and a column
  one too high. For a 100 x 100 raster that is 100 of the 10,000 positions.

  This was the main cause of the "chosen nucleus outside of parent" condition:
  a nucleus decoded to row 0 made the distance matrix be measured from outside
  the raster, so the `bad_nucleus` check fired. Correcting the decode removed it
  at the default resolution, where it had fired on roughly 3% of runs for
  polyclonal data. A second cause remained at coarse resolutions and is fixed
  separately below. The retry introduced in 1.0.1 is kept as a backstop.
  Because affected nucleus draws were placed wrongly before, clone positions
  shift for that small fraction of draws, so seeded output is not identical to
  1.0.1.

  The unused inverse `coordinates.to.matrix.index()` was row-major to match,
  and so did not invert the decode; it is now consistent and the two functions
  round-trip over all indices.

* `cloneMap()` no longer aborts when no candidate nucleus position is
  available. Daughter nuclei are drawn from an annulus at a set distance inside
  the parent clone, and for a small clone or a coarse `resolution.index` that
  annulus can contain no positions at all. The draw was made with
  `sample(1:n, ...)`, which for `n = 0` counts down to `c(1, 0)` and returns an
  index outside the candidate set, placing the nucleus outside the parent. The
  1.0.1 retry could not recover, because the candidate set was empty on every
  attempt: the failure was determined by the input and the seed rather than
  being transient. An empty annulus now falls back to the positions available
  within the parent, taking the most central ones when only one daughter is
  being placed. Output is unchanged wherever a candidate position existed
  before, verified `identical()` across 54 combinations of nine inputs and six
  seeds, and unrooted example data that previously failed at
  `resolution.index = 20` now completes on all of 120 seeds.

* Fixed `remove.clones.on.tree()` destroying the tree when called with
  `clones.to.remove` rather than `clones.to.keep`. In that case
  `clones.to.keep` is empty internally, and the guard that checks whether only
  the root clone remains tested it with `all()`, which is vacuously `TRUE` on an
  empty vector. The function therefore returned a single root -> root row and
  discarded the rest of the tree, so the `clones.to.remove` route was unusable
  for pruning. The `clones.to.keep` route was never affected, and because both
  internal calls inside `cloneMap()` pass `clones.to.keep`, no plot output was
  ever affected.

## Performance

* The clone growth model is about 16x faster with unchanged output. The row and
  column distance matrices in `make.distance.matrix()` are outer products and
  are now computed with `outer()` rather than by `cbind()`/`rbind()` over every
  row and column; candidate nucleus spacing is scored with a single
  `dist(method = "manhattan")` call on the candidate coordinates instead of a
  full `resolution.index^2` distance matrix per nucleus; and the growth cut off
  search counts positions with `findInterval()` on the sorted distances instead
  of rescanning the whole raster once per candidate cut off. Verified
  `identical()` over 54 combinations of 9 inputs and 6 seeds: 48.5s to 3.0s of
  compute, with byte-identical clone map objects.

* The clone continuity test is 20x to 600x faster. It previously built an
  eight-element character vector of neighbour ids for every position in a clone
  and grew the connected set by rescanning with `%in%` on each pass, which is
  quadratic in clone size. It now finds edges with shifted comparisons on a
  padded copy of the matrix and flood fills over integer matrix indices one
  front at a time. A 10,000 cell clone goes from 1.9s to under 5ms. Verdicts
  are identical on 20 test shapes, and the random draw that picks a starting
  edge is retained so the random number stream is unchanged.

* Plot extent setup no longer converts the whole raster to polygons, which for
  a 100 x 100 raster meant creating and invisibly drawing 10,000 single cell
  polygons when only the bounding box is used. A single rectangle from the
  raster extent gives an identical plot region (`par("usr")` and `par("plt")`
  match to the value) and saves about 0.35s per `cloneMap()` call. The raster
  construction is also hoisted out of the per-daughter plotting loop, where it
  was rebuilt identically for every clone.

## Documentation

* Added a `\value` section to the help page of every exported function, and
  expanded that of `cloneMap()` to describe the components of the returned
  clone map object and to state that the plotting path returns `NULL`
  invisibly. The documented but unexported `extract_daughters()` is now marked
  internal.

* Corrected the documentation of the `clone_names` argument of
  `make.CCFs.tree.consistant()`. It was described as restricting the output to
  the named clones. It does not: it is a lookup table mapping the internal clone
  names used during plotting to the original user-supplied ones, and is used
  only to report original names in warning messages. The returned values do not
  depend on it.

* Restructured the `cloneMap()` examples so that the example block fits within
  the check time budget. Two representative calls, one rooted and one unrooted,
  run live and the remaining variants are under `\donttest{}`. Note that
  `\donttest{}` does not exempt code from the time budget under `--as-cran`,
  so the live examples also pass `resolution.index = 40`, which renders every
  clone in all three example datasets without "clone missing from raster"
  warnings.

## Internal

* Renamed the internal `min.nuclei.distance()` to `min_nuclei_distance()`.
  `R CMD check` read the dotted name as an S3 method for `min()` and reported a
  generic/method consistency mismatch. The function is internal, so no user
  code is affected.

* Declared `methods`, `stats` and `utils` in `Imports`. `stats::quantile()` and
  `utils::txtProgressBar()` were already used without being declared, and
  `methods::as()` is used by the new plot extent helper.

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
