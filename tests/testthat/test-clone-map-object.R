# The "Clone map" object is cloneMap()'s reusable, reproducible output: it holds
# the rasterised clone positions so the same map can be replotted without
# re-running the semi-random placement.
#
# resolution.index is kept low throughout: these tests check structure, not
# picture quality, and placement cost grows quickly with resolution.

test_that( "output.Clone.map.obj returns a Clone map object with the documented components", {

  obj <- quietly( cloneMap( tree_example_2, CCFs_example_2,
                            output.Clone.map.obj = TRUE, plot.data = FALSE,
                            resolution.index = 20 ) )

  expect_s3_class( obj, "Clone map" )
  expect_type( obj, "list" )

  expect_true( all( c( "clone_matrix", "tree", "tree_internal", "names_match",
                       "CCFs", "shannon_diversity", "simpson_diversity" )
                    %in% names( obj ) ) )

  # the raster of clone positions, square at the requested resolution
  expect_true( is.matrix( obj$clone_matrix ) )
  expect_equal( dim( obj$clone_matrix ), c( 20L, 20L ) )

  # both views of the tree, and the mapping between them
  expect_true( is.matrix( obj$tree ) )
  expect_true( is.matrix( obj$tree_internal ) )
  expect_equal( dim( obj$tree ), dim( obj$tree_internal ) )
  expect_true( all( c( "orig", "internal" ) %in% names( obj$names_match ) ) )

  # the corrected CCF table, carrying the clone areas used for placement
  expect_s3_class( obj$CCFs, "data.frame" )
  expect_true( all( c( "clones", "CCF" ) %in% names( obj$CCFs ) ) )

  # precomputed diversity of the map
  for( d in c( obj$shannon_diversity, obj$simpson_diversity ) ){
    expect_length( d, 1 )
    expect_true( is.finite( d ) )
  }

})

test_that( "the object's tree is in original clone names and its raster holds those clones", {

  obj <- quietly( cloneMap( tree_example_2, CCFs_example_2,
                            output.Clone.map.obj = TRUE, plot.data = FALSE,
                            resolution.index = 20 ) )

  # $tree is expressed in the user's clone names, $tree_internal in the
  # internal ones used during placement
  expect_setequal( unique( as.character( obj$tree ) ),
                   as.character( unique( as.numeric( tree_example_2 ) ) ) )

  # every value in the raster is an internal clone id (or 0 for empty space)
  raster_values <- setdiff( unique( as.numeric( obj$clone_matrix ) ), 0 )
  expect_true( all( raster_values %in% as.numeric( obj$names_match$internal ) ) )

})

test_that( "plot.data = FALSE draws nothing and leaves no files behind", {

  files_before <- list.files()

  obj <- quietly( cloneMap( tree_example_2, CCFs_example_2,
                            output.Clone.map.obj = TRUE, plot.data = FALSE,
                            resolution.index = 20 ) )

  expect_s3_class( obj, "Clone map" )
  expect_setequal( list.files(), files_before )

  # no graphics device was opened, so nothing could have been drawn
  # (device 1 is the null device, i.e. no device is active)
  expect_equal( unname( grDevices::dev.cur() ), 1L )

})

test_that( "tissue_border draws while establishing the tissue area, even with plot.data = FALSE", {

  # Documented so the behaviour is not mistaken for a test bug: the border of
  # the tissue area is drawn as the plot area is being established, upstream of
  # the plot.data switch. Building an object from unrooted data with
  # tissue_border = TRUE therefore opens a device and, on a bare session, would
  # write Rplots.pdf into the working directory.
  files_before <- list.files()

  obj <- with_null_device( quietly(
    cloneMap( tree.mat = tree_example_poly, CCF.data = CCF_example_poly,
              tissue_border = TRUE, output.Clone.map.obj = TRUE,
              plot.data = FALSE, resolution.index = 25 ) ) )

  expect_s3_class( obj, "Clone map" )

  # with the drawing sent to the null device, still no detritus
  expect_setequal( list.files(), files_before )
  expect_false( file.exists( "Rplots.pdf" ) )

})

test_that( "a Clone map object can be passed back to reproduce its plot", {

  obj <- quietly( cloneMap( tree_example_2, CCFs_example_2,
                            output.Clone.map.obj = TRUE, plot.data = FALSE,
                            resolution.index = 20 ) )

  files_before <- list.files()

  expect_no_error(
    with_null_device( quietly( cloneMap( clone_map = obj ) ) ) )

  # replotting must not write anything into the package directory
  expect_setequal( list.files(), files_before )

})

test_that( "replotting the same object twice gives the same clone positions", {

  # this is the reproducibility guarantee of the object: placement already
  # happened, so the raster is fixed and replotting cannot move clones
  obj <- quietly( cloneMap( tree_example_2, CCFs_example_2,
                            output.Clone.map.obj = TRUE, plot.data = FALSE,
                            resolution.index = 20 ) )

  snapshot <- obj$clone_matrix

  with_null_device( quietly( cloneMap( clone_map = obj ) ) )
  with_null_device( quietly( cloneMap( clone_map = obj ) ) )

  expect_identical( obj$clone_matrix, snapshot )

})

test_that( "cloneMap needs either a tree plus CCFs or a Clone map object", {

  # suppressMessages keeps the progress signpost cloneMap() prints before
  # validating its input out of the test output
  expect_error( suppressMessages( cloneMap() ), "Please provide either" )

  # something that is not a Clone map object is refused
  expect_error( suppressMessages( cloneMap( clone_map = list( nonsense = TRUE ) ) ),
                "incorrect raster input" )

})

test_that( "a plain plotting call runs and leaves no files behind", {

  files_before <- list.files()

  expect_no_error(
    with_null_device( quietly(
      cloneMap( tree_example_1, CCFs_example_1, resolution.index = 20 ) ) ) )

  expect_setequal( list.files(), files_before )

})
