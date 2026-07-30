# Clone placement is semi-random: cloneMap() draws nuclei positions with
# sample() and grows clones outward from them. A bug in the nucleus check used
# to abort with "(BUG) chosen nucleus outside of parent" on a small fraction of
# random draws for unrooted data, so the guarantee worth testing is that many
# independent draws all complete.
#
# These tests deliberately do NOT compare pixels. Placement is semi-random and
# platform-dependent, so a snapshot of figure output would be meaningless on a
# different machine. What is asserted is that no draw errors and that the
# structural invariants of the result hold.
#
# A second, related cause was found and fixed while these tests were written:
# at low resolution.index the smallest clones cover only a handful of raster
# cells, so the annulus of candidate nucleus positions could come up empty.
# sample() was then called on an empty set via the 1:n idiom, which yields an
# out of range index and places a nucleus outside the parent, and no number of
# retries could recover because the candidate set was empty on every attempt.
# Empty candidate sets now fall back to the positions available inside the
# parent, so coarse rasters complete too; test-tree-helpers.R and the seed loops
# below cover this.

test_that( "unrooted data with a tissue border completes across many random draws", {

  n_seeds <- 30

  with_null_device({

    for( seed in seq_len( n_seeds ) ){

      set.seed( seed )

      expect_no_error(
        quietly( cloneMap( tree.mat = tree_example_poly,
                           CCF.data = CCF_example_poly,
                           tissue_border = TRUE,
                           resolution.index = 25 ) ),
        message = paste( "cloneMap() failed with seed", seed ) )

    }

  })

})

test_that( "coarse rasters complete, where an empty nucleus annulus used to abort", {

  # regression test for the empty candidate set fix. At resolution.index = 20
  # seed 10 aborted with "chosen nucleus outside of parent for clone 13", and
  # it was not retry-recoverable: the annulus of candidate positions was empty
  # on every attempt, so sample() drew an out of range index each time. This
  # was independent of tissue_border. Resolutions this coarse are not what you
  # would plot a figure at, but they must not error.

  with_null_device({

    for( border in c( TRUE, FALSE ) ){

      set.seed( 10 )

      expect_no_error(
        quietly( cloneMap( tree.mat = tree_example_poly,
                           CCF.data = CCF_example_poly,
                           tissue_border = border,
                           resolution.index = 20 ) ),
        message = paste( "seed 10 at resolution.index 20, tissue_border =", border ) )

    }

    # and a sweep of coarse resolutions, which stress the same path
    for( res in c( 12, 15, 20 ) ){

      for( seed in seq_len( 5 ) ){

        set.seed( seed )

        expect_no_error(
          quietly( cloneMap( tree.mat = tree_example_poly,
                             CCF.data = CCF_example_poly,
                             resolution.index = res ) ),
          message = paste( "resolution.index", res, "seed", seed ) )

      }

    }

  })

})

test_that( "unrooted data without a tissue border also completes across many draws", {

  with_null_device({

    for( seed in seq_len( 15 ) ){

      set.seed( seed )

      expect_no_error(
        quietly( cloneMap( tree.mat = tree_example_poly,
                           CCF.data = CCF_example_poly,
                           resolution.index = 25 ) ),
        message = paste( "cloneMap() failed with seed", seed ) )

    }

  })

})

test_that( "every random draw yields a structurally valid Clone map object", {

  # a semi-random placement must still respect the object contract and place
  # every clone the CCF table asked for
  for( seed in seq_len( 5 ) ){

    set.seed( seed )

    # the null device is needed even with plot.data = FALSE: the tissue border
    # is drawn while the tissue area is being established, before the
    # plot.data switch is consulted
    obj <- with_null_device( quietly( cloneMap( tree.mat = tree_example_poly,
                              CCF.data = CCF_example_poly,
                              tissue_border = TRUE,
                              output.Clone.map.obj = TRUE, plot.data = FALSE,
                              resolution.index = 25 ) ) )

    expect_s3_class( obj, "Clone map" )
    expect_equal( dim( obj$clone_matrix ), c( 25L, 25L ) )

    # The raster carries clone ids plus two sentinels: 1000 marks tissue inside
    # the plot area that no clone reached, and 0 marks positions outside the
    # tissue area altogether. Everything else must be a real internal clone id.
    clone_ids <- setdiff( unique( as.numeric( obj$clone_matrix ) ), c( 0, 1000 ) )
    expect_true( all( clone_ids %in% as.numeric( obj$names_match$internal ) ) )
    expect_false( any( is.na( obj$clone_matrix ) ) )

    # clones occupy some but not all of the tissue: these CCFs total well
    # under 100%, so unreached tissue (the 1000 sentinel) must remain
    tissue <- obj$clone_matrix != 0
    occupied_of_tissue <- sum( !obj$clone_matrix %in% c( 0, 1000 ) ) / sum( tissue )

    expect_gt( occupied_of_tissue, 0 )
    expect_lt( occupied_of_tissue, 1 )

  }

})

test_that( "space_fraction increases the white space left inside the tissue", {

  # space_fraction = 0.7 asks that 70% of the plot area be white space, so a
  # smaller share of the tissue should be covered by clones than under the
  # default, for the same random draw
  occupied_of_tissue <- function( obj ){
    sum( !obj$clone_matrix %in% c( 0, 1000 ) ) / sum( obj$clone_matrix != 0 )
  }

  for( seed in seq_len( 3 ) ){

    set.seed( seed )
    default <- with_null_device( quietly( cloneMap( tree.mat = tree_example_poly,
                                  CCF.data = CCF_example_poly,
                                  tissue_border = TRUE,
                                  output.Clone.map.obj = TRUE, plot.data = FALSE,
                                  resolution.index = 25 ) ) )

    set.seed( seed )
    sparse <- with_null_device( quietly( cloneMap( tree.mat = tree_example_poly,
                                 CCF.data = CCF_example_poly,
                                 tissue_border = TRUE, space_fraction = 0.7,
                                 output.Clone.map.obj = TRUE, plot.data = FALSE,
                                 resolution.index = 25 ) ) )

    expect_lt( occupied_of_tissue( sparse ), occupied_of_tissue( default ) )
    expect_gt( occupied_of_tissue( sparse ), 0 )

  }

})

test_that( "repeated draws actually differ, confirming placement is semi-random", {

  # if this ever stopped holding, the tests above would no longer be exercising
  # independent random draws at all
  set.seed( 1 )
  first <- quietly( cloneMap( tree_example_2, CCFs_example_2,
                              output.Clone.map.obj = TRUE, plot.data = FALSE,
                              resolution.index = 20 ) )$clone_matrix
  set.seed( 2 )
  second <- quietly( cloneMap( tree_example_2, CCFs_example_2,
                               output.Clone.map.obj = TRUE, plot.data = FALSE,
                               resolution.index = 20 ) )$clone_matrix

  expect_false( identical( first, second ) )

})

test_that( "the test suite leaves no plot files behind", {

  # CRAN checks for detritus in the package directory: every plotting call in
  # this suite is wrapped so that nothing is written
  expect_false( file.exists( "Rplots.pdf" ) )

})
