# make_clone_col_input() builds the named colour vector cloneMap() expects.

test_that( "one named hex colour is returned per clone", {

  clones <- as.character( 1:5 )
  cols <- make_clone_col_input( clones )

  expect_length( cols, length( clones ) )
  expect_identical( names( cols ), clones )
  expect_true( all( grepl( "^#[0-9A-Fa-f]{6}$", cols ) ) )
  expect_false( any( duplicated( cols ) ) )

})

test_that( "fewer than three clones still get the right number of colours", {

  # RColorBrewer cannot return fewer than 3 colours, so the function requests 3
  # and trims. The edge case is that the trim must match the clone count.
  two <- make_clone_col_input( c( "a", "b" ) )
  expect_length( two, 2 )
  expect_identical( names( two ), c( "a", "b" ) )
  expect_true( all( grepl( "^#[0-9A-Fa-f]{6}$", two ) ) )

  one <- make_clone_col_input( "solo" )
  expect_length( one, 1 )
  expect_identical( names( one ), "solo" )
  expect_true( grepl( "^#[0-9A-Fa-f]{6}$", one ) )

})

test_that( "palettes are extended past the brewer maximum without warning", {

  # brewer palettes top out at 12 colours, beyond which the function ramps
  many <- as.character( 1:20 )
  cols <- expect_silent( make_clone_col_input( many ) )

  expect_length( cols, 20 )
  expect_identical( names( cols ), many )
  expect_true( all( grepl( "^#[0-9A-Fa-f]{6}$", cols ) ) )

})

test_that( "clone names from a tree are preserved exactly", {

  clone_names <- unique( c( tree_example[, 1 ], tree_example[, 2 ] ) )
  cols <- make_clone_col_input( clone_names )

  expect_identical( names( cols ), as.character( clone_names ) )
  expect_length( cols, length( clone_names ) )

})

test_that( "an alternative brewer palette is honoured", {

  paired <- make_clone_col_input( c( "a", "b", "c" ) )
  dark <- make_clone_col_input( c( "a", "b", "c" ), brewer.palette = "Dark2" )

  expect_length( dark, 3 )
  expect_false( identical( unname( paired ), unname( dark ) ) )

})
