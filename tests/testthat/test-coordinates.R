# Regression tests for the column-major coordinate helpers.
#
# matrix.index.to.coordinates() and coordinates.to.matrix.index() are inverses,
# and R stores matrices column-major: index i sits at row ((i-1) %% nrow) + 1
# and column ((i-1) %/% nrow) + 1. A decoding bug that swapped rows for columns
# was fixed previously; a square test matrix would not catch a transposition,
# so every case below uses nrow != ncol.

test_that( "decoding then re-encoding every index is the identity", {

  nr <- 4L; nc <- 7L
  indices <- seq_len( nr * nc )

  coords <- cloneMap:::matrix.index.to.coordinates( indices, nrow = nr, ncol = nc )

  expect_s3_class( coords, "data.frame" )
  expect_identical( names( coords ), c( "x", "y" ) )
  expect_identical( nrow( coords ), length( indices ) )

  round_trip <- cloneMap:::coordinates.to.matrix.index( coords, nrow = nr, ncol = nc )

  expect_equal( round_trip, indices )

})

test_that( "decoding matches R's own column-major convention", {

  nr <- 4L; nc <- 7L
  indices <- seq_len( nr * nc )

  coords <- cloneMap:::matrix.index.to.coordinates( indices, nrow = nr, ncol = nc )

  # arrayInd() is R's reference implementation of the same decoding
  reference <- arrayInd( indices, .dim = c( nr, nc ) )

  expect_equal( coords$x, reference[, 1 ] )
  expect_equal( coords$y, reference[, 2 ] )

  # x is the ROW and y is the COLUMN: coordinates never exceed their extent
  expect_true( all( coords$x >= 1 & coords$x <= nr ) )
  expect_true( all( coords$y >= 1 & coords$y <= nc ) )

})

test_that( "encoding agrees with which() on a marked matrix", {

  nr <- 4L; nc <- 7L

  # mark one cell at a time and check the helper recovers the index which()
  # reports, for every cell of a non-square matrix.
  # Coordinates are passed as doubles: the scalar branch of the encoder tests
  # class(coordinates) == "numeric", which an integer vector does not satisfy.
  recovered <- vapply( seq_len( nr ), function( r ){
    vapply( seq_len( nc ), function( cc ){
      m <- matrix( FALSE, nr, nc )
      m[ r, cc ] <- TRUE
      encoded <- cloneMap:::coordinates.to.matrix.index(
        as.numeric( c( r, cc ) ), nrow = nr, ncol = nc )
      isTRUE( encoded == which( m ) )
    }, logical( 1 ) )
  }, logical( nc ) )

  expect_true( all( recovered ) )

})

test_that( "the data.frame branch handles many coordinates at once", {

  nr <- 4L; nc <- 7L

  # the encoder accepts either a single numeric coordinate pair or a whole
  # data.frame of them, and the two routes must agree
  coords <- data.frame( x = c( 1, 4, 2, 3 ), y = c( 1, 7, 5, 2 ) )

  vectorised <- cloneMap:::coordinates.to.matrix.index( coords, nrow = nr, ncol = nc )
  one_at_a_time <- vapply( seq_len( nrow( coords ) ), function( i ){
    cloneMap:::coordinates.to.matrix.index(
      as.numeric( c( coords$x[ i ], coords$y[ i ] ) ), nrow = nr, ncol = nc )
  }, numeric( 1 ) )

  expect_equal( vectorised, one_at_a_time )

  # and both agree with R's own column-major indexing
  expect_equal( vectorised, ( coords$y - 1 ) * nr + coords$x )

})

test_that( "a transposed decoding would be caught", {

  # index 5 of a 4 x 7 matrix is row 1, column 2. Under a row-major (or
  # transposed) reading it would decode to row 2 column 1 instead, so this
  # single assertion pins the orientation.
  coords <- cloneMap:::matrix.index.to.coordinates( 5, nrow = 4, ncol = 7 )

  expect_equal( coords$x, 1 )
  expect_equal( coords$y, 2 )

  # and the scalar (numeric vector) branch of the encoder round-trips it
  expect_equal( cloneMap:::coordinates.to.matrix.index( c( 1, 2 ), nrow = 4, ncol = 7 ), 5 )

})
