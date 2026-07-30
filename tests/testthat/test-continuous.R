# continuous.test() reports whether a clone occupies one connected patch or has
# broken into separate islands. Connectivity is 8-way (diagonal neighbours
# count), and only the clone's edge positions are walked.

test_that( "a single contiguous block is continuous", {

  clone <- matrix( FALSE, 10, 10 )
  clone[ 3:6, 3:6 ] <- TRUE

  expect_true( cloneMap:::continuous.test( clone ) )

})

test_that( "a clone split into two islands is not continuous", {

  clone <- matrix( FALSE, 10, 10 )
  clone[ 2:3, 2:3 ] <- TRUE
  clone[ 8:9, 8:9 ] <- TRUE

  expect_false( cloneMap:::continuous.test( clone ) )

})

test_that( "degenerate clone positions are treated as continuous", {

  # no positions at all: nothing to disconnect
  expect_true( cloneMap:::continuous.test( matrix( FALSE, 5, 5 ) ) )

  # a single occupied cell
  one <- matrix( FALSE, 10, 10 )
  one[ 5, 5 ] <- TRUE
  expect_true( cloneMap:::continuous.test( one ) )

})

test_that( "connectivity is 8-way, so diagonally touching blocks are continuous", {

  clone <- matrix( FALSE, 6, 6 )
  clone[ 2:3, 2:3 ] <- TRUE
  clone[ 4:5, 4:5 ] <- TRUE   # touches the first block only at a corner

  expect_true( cloneMap:::continuous.test( clone ) )

})

test_that( "the result does not depend on which edge position is drawn", {

  # the function picks its starting edge position with sample(); which one is
  # chosen cannot change the answer, so repeated calls must agree
  clone <- matrix( FALSE, 12, 12 )
  clone[ 2:4, 2:4 ] <- TRUE
  clone[ 9:11, 9:11 ] <- TRUE

  expect_true( all( !replicate( 10, cloneMap:::continuous.test( clone ) ) ) )

})
