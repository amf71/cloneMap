## Submission

This is a first submission. cloneMap has not previously been on CRAN.

## Test environments

* local macOS install (arm64, Apple silicon), R 4.5.3
* win-builder, R devel -- NOT YET RUN
* win-builder, R release -- NOT YET RUN
* R-hub, Linux (Ubuntu, R release) -- NOT YET RUN

## R CMD check results

The local run reported 1 ERROR, 1 WARNING and 4 NOTEs. The ERROR and the
WARNING are both the PDF reference manual failing to typeset, because the
local machine has no LaTeX installation; they are not package faults and are
detailed under "Local environment artefacts" below. No ERROR or WARNING arose
from the package itself.

Two NOTEs are expected to remain on a clean platform:

* checking CRAN incoming feasibility ... NOTE

      New submission

  Expected: the package has not been on CRAN before.

* checking CRAN incoming feasibility ... NOTE

      Found the following (possibly) invalid ORCID iD:
        iD: 0000-0002-0341-7878  (from: DESCRIPTION)

  The iD is correct. Its check digit was verified against the ISO 7064
  MOD 11-2 scheme that ORCID uses, and it validates. The note is raised
  because the machine used for local testing had no network access to
  orcid.org, so the validator could not confirm the iD against the
  register.

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
