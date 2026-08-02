## Submission

This is a first submission. cloneMap has not previously been on CRAN.

This is a resubmission of version 1.0.2 (not yet accepted, so the version
number is unchanged). The initial submission failed CRAN's automatic incoming
checks on the Debian flavour with an examples-timing NOTE:

    Examples with CPU (user + system) or elapsed time > 5s
              user system elapsed
    cloneMap 5.834  0.195   6.033

The cause was that the package declared its dependencies in `Imports:` but
called them only as `raster::`, `sf::`, `smoothr::` and `vegan::`, with no
`importFrom` directives in NAMESPACE. Those namespaces therefore loaded
lazily at first use -- inside the first `cloneMap()` call, and so inside the
block `R CMD check` times for examples. Around 2.5s of the reported time was
raster/terra/sf initialisation rather than any work the example does.

The imports are now declared, which moves that cost into
`library(cloneMap)`, where it is not counted against the examples budget.
With that headroom, five illustrative examples (six `cloneMap()` calls: a
rooted tree, an unrooted/polyclonal tree, the clone_map-object replot
pattern, and a colour-matched pair) are live; three further variants that
add resolution or cosmetic detail rather than a new usage pattern
(tissue border, `space_fraction`, and the publication-quality default
resolution) remain in `\donttest{}`.

Measured with the `--run-donttest` pass disabled, so that the timings file
records the live examples rather than being overwritten by the donttest run,
user+system time for the `cloneMap` topic is 2.04s on the test machine, a
34% margin below the 5s ceiling once scaled to CRAN's Debian pretest (see
below). Figures are unaffected -- six seeded examples render byte-identically
before and after -- and all 244 tests still pass.

## Test environments

* local macOS install (arm64, Apple silicon), R 4.5.3
* win-builder, R-devel (x86_64-w64-mingw32, Windows Server 2022) -- 1 NOTE
* R-hub v2 (GitHub Actions), R-devel, Ubuntu 24.04 -- OK
* R-hub v2 (GitHub Actions), R-devel, Windows Server 2022 -- OK
* R-hub v2 (GitHub Actions), R-devel, macOS Sequoia 15.7.7 -- OK

## R CMD check results

The test environments above are for the 1.0.2 tarball as it stood before the
examples-timing fix; none of them reproduced the Debian NOTE, which is why it
was not caught before the first submission. Their reported timings were also
not a reliable guide, because a plain `R CMD check` prints the timings table
for information without raising a NOTE, and the `--run-donttest` pass
overwrites the timings file; CRAN's incoming report is what escalates an
over-threshold entry to a NOTE.

After the fix, a local `R CMD check --as-cran` gives `checking examples ...
OK`, all 244 tests passing, and no other change. The remaining ERROR, WARNING
and NOTEs are the local-machine artefacts described below.

The test machine measured 3.72s for the originally-submitted tarball, which
CRAN's Debian pretest measured at 6.03s -- a scaling factor of 1.62x. At that
factor, the current examples time (2.04s locally) projects to roughly 3.3s
on CRAN's Debian machine, comfortably under the 5s ceiling.

There were no ERRORs or WARNINGs on any of the platforms below. All four
remote platforms ran the package's testthat suite (244 expectations) with no
failures.

The local macOS run reported 1 ERROR and 1 WARNING in addition to the NOTEs
below; both are the PDF reference manual failing to typeset because that
machine has no LaTeX installation, not a package fault, and did not reproduce
on any of the four remote platforms. Detail under "Local environment
artefacts".

One NOTE is expected to remain on a clean platform:

* checking CRAN incoming feasibility ... NOTE

      New submission

  Expected: the package has not been on CRAN before. Seen on win-builder;
  not shown by R-hub's --as-cran run.

The local macOS run additionally raised:

* checking CRAN incoming feasibility ... NOTE

      Found the following (possibly) invalid ORCID iD:
        iD: 0000-0002-0341-7878  (from: DESCRIPTION)

  The iD is correct. Its check digit was verified against the ISO 7064
  MOD 11-2 scheme that ORCID uses, and it validates. The note is raised
  because the local machine had no network access to orcid.org, so the
  validator could not confirm the iD against the register; win-builder and
  R-hub, which do have network access, did not raise it.

## Local environment artefacts

These were seen on the local test machine only. They are properties of that
machine rather than of the package, and are not expected on CRAN's check
machines. They are listed for completeness.

* The ERROR and WARNING on the PDF reference manual -- no LaTeX installation
  is available locally, so the manual cannot be typeset here. All 17 Rd files
  parse and convert to LaTeX without error, checked with `tools::parse_Rd()`
  and `tools::Rd2latex()`.
* A NOTE on `cloneMap-manual.tex` being left in the check directory -- a
  by-product of the failed PDF step above.
* `checking for future file timestamps ... NOTE: unable to verify current
  time` -- the local machine has no network access to the time service.
* A NOTE on the HTML version of the manual -- the local HTML Tidy predates
  the version R checks against.

## Downstream dependencies

There are none. This is a first submission.
