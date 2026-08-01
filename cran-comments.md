## Submission

This is a first submission. cloneMap has not previously been on CRAN.

## Test environments

* local macOS install (arm64, Apple silicon), R 4.5.3
* win-builder, R-devel (x86_64-w64-mingw32, Windows Server 2022) -- 1 NOTE
* R-hub v2 (GitHub Actions), R-devel, Ubuntu 24.04 -- OK
* R-hub v2 (GitHub Actions), R-devel, Windows Server 2022 -- OK
* R-hub v2 (GitHub Actions), R-devel, macOS Sequoia 15.7.7 -- OK

## R CMD check results

There were no ERRORs or WARNINGs on any platform. All four remote platforms
above ran the package's testthat suite (244 expectations) with no failures.

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
