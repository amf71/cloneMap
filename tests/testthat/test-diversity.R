# calc_clonal_diversity() converts CCFs into clonal prevalences (the fraction of
# cells carrying exactly that genotype) and passes them to vegan::diversity().

test_that( "shannon and simpson return a single finite number", {

  for( index in c( "shannon", "simpson" ) ){

    diversity <- suppressWarnings(
      calc_clonal_diversity( CCFs_example_2, tree_example, method = index ) )

    expect_type( diversity, "double" )
    expect_length( diversity, 1 )
    expect_true( is.finite( diversity ) )
    expect_gte( diversity, 0 )

  }

  # simpson is an index bounded above by 1; shannon is not
  simpson <- suppressWarnings(
    calc_clonal_diversity( CCFs_example_2, tree_example, method = "simpson" ) )
  expect_lte( simpson, 1 )

})

test_that( "a single clone tree has zero diversity", {

  # one clone occupying everything: no diversity of prevalences at all.
  # Both indices are 0 here (not NA and not NaN).
  solo_tree <- matrix( c( 1, 1 ), ncol = 2 )
  solo_CCFs <- data.frame( clones = 1, CCF = 1, stringsAsFactors = FALSE )

  shannon <- suppressWarnings(
    calc_clonal_diversity( solo_CCFs, solo_tree, method = "shannon" ) )
  simpson <- suppressWarnings(
    calc_clonal_diversity( solo_CCFs, solo_tree, method = "simpson" ) )

  expect_equal( shannon, 0 )
  expect_equal( simpson, 0 )

})

test_that( "more evenly spread clones are more diverse", {

  # one dominant clone with a small subclone vs an even split: the even split
  # must score higher on both indices
  tree <- matrix( c( 1, 2, 1, 3 ), ncol = 2, byrow = TRUE )

  skewed <- data.frame( clones = c( 1, 2, 3 ), CCF = c( 1, 0.02, 0.02 ),
                        stringsAsFactors = FALSE )
  even <- data.frame( clones = c( 1, 2, 3 ), CCF = c( 1, 0.45, 0.45 ),
                      stringsAsFactors = FALSE )

  for( index in c( "shannon", "simpson" ) ){
    expect_lt( suppressWarnings( calc_clonal_diversity( skewed, tree, method = index ) ),
               suppressWarnings( calc_clonal_diversity( even, tree, method = index ) ) )
  }

})

test_that( "the inc_parents_ccf switch is respected without error", {

  # the two settings resolve CCF/tree inconsistencies in opposite directions,
  # so both must produce a usable number
  for( inc in c( TRUE, FALSE ) ){

    diversity <- suppressWarnings(
      calc_clonal_diversity( CCFs_example_1, tree_example_1, inc_parents_ccf = inc ) )

    expect_length( diversity, 1 )
    expect_true( is.finite( diversity ) )

  }

})

test_that( "diversity of an unrooted tree is finite", {

  diversity <- suppressWarnings(
    calc_clonal_diversity( CCF_example_poly, tree_example_poly, method = "shannon" ) )

  expect_length( diversity, 1 )
  expect_true( is.finite( diversity ) )
  expect_gt( diversity, 0 )

})
