# Tests for the exported tree manipulation helpers.

test_that( "find_root returns the single root of a rooted tree", {

  root <- find_root( tree_example )

  expect_length( root, 1 )
  expect_equal( root, 1 )

  # the root is the only clone that never appears as a daughter
  expect_false( root %in% tree_example[, 2 ] )

})

test_that( "find_root returns every root of an unrooted (polyclonal) tree", {

  roots <- find_root( tree_example_poly )

  expect_gt( length( roots ), 1 )

  # clones with neither parent nor daughter are written clone -> clone and are
  # roots in their own right, alongside the true branch roots 3 and 7
  expect_setequal( roots, c( 1, 2, 3, 7, 9, 10, 11, 12, 13, 14 ) )

  # no root is ever some other clone's daughter
  daughters_of_others <- tree_example_poly[ tree_example_poly[, 1 ] != tree_example_poly[, 2 ], 2 ]
  expect_false( any( roots %in% daughters_of_others ) )

})

test_that( "logically.order.tree puts every parent before its daughters", {

  ordered <- quietly( logically.order.tree( tree_example ) )

  # same relationships, just reordered
  expect_equal( dim( ordered ), dim( tree_example ) )
  expect_setequal( apply( ordered, 1, paste, collapse = "->" ),
                   apply( tree_example, 1, paste, collapse = "->" ) )

  # the ordering invariant: a clone appears as a daughter only after the row in
  # which it appears as a parent, except roots which have no parent row
  expect_true( parents_precede_daughters( ordered ) )

  # the root leads the ordering
  expect_equal( ordered[ 1, 1 ], find_root( tree_example ) )

})

test_that( "logically.order.tree holds the invariant for an unrooted tree", {

  ordered <- quietly( logically.order.tree( tree_example_poly ) )

  expect_setequal( apply( ordered, 1, paste, collapse = "->" ),
                   apply( tree_example_poly, 1, paste, collapse = "->" ) )
  expect_true( parents_precede_daughters( ordered ) )

})

test_that( "logically.order.tree messages when handed a polyclonal tree", {

  # the message is how the function flags a tree with no common root
  expect_message( logically.order.tree( tree_example_poly ), "polyclonal" )

})

test_that( "remove.clones.on.tree reattaches daughters when removing an intermediate clone", {

  # clone 2 is an intermediate: 1 -> 2 -> 4, plus a sibling branch 1 -> 3.
  # Removing it must reattach 4 to its grandparent 1 rather than drop it.
  tr <- matrix( c( 1, 2,
                   2, 4,
                   1, 3 ), ncol = 2, byrow = TRUE )

  pruned <- quietly( remove.clones.on.tree( tr, clones.to.remove = 2 ) )

  edges <- apply( pruned, 1, paste, collapse = "->" )

  expect_setequal( edges, c( "1->3", "1->4" ) )

  # clone 2 is gone entirely, and no clone was lost along with it
  expect_false( 2 %in% as.numeric( pruned ) )
  expect_setequal( unique( as.numeric( pruned ) ), c( 1, 3, 4 ) )

})

test_that( "remove.clones.on.tree gives the same result via clones.to.keep", {

  tr <- matrix( c( 1, 2,
                   2, 4,
                   1, 3 ), ncol = 2, byrow = TRUE )

  by_remove <- quietly( remove.clones.on.tree( tr, clones.to.remove = 2 ) )
  by_keep   <- quietly( remove.clones.on.tree( tr, clones.to.keep = c( 1, 3, 4 ) ) )

  expect_setequal( apply( by_remove, 1, paste, collapse = "->" ),
                   apply( by_keep,   1, paste, collapse = "->" ) )

})

test_that( "remove.clones.on.tree does not discard the tree when clones.to.remove is used", {

  # Regression test. clones.to.keep is empty on the clones.to.remove route, and
  # all() on an empty vector is vacuously TRUE, so an unguarded root check
  # collapsed the whole tree to a bare root -> root matrix.
  tr <- matrix( c( 1, 2,
                   2, 4,
                   1, 3 ), ncol = 2, byrow = TRUE )

  pruned <- quietly( remove.clones.on.tree( tr, clones.to.remove = 2 ) )

  expect_equal( nrow( pruned ), 2 )
  expect_false( identical( unique( as.numeric( pruned ) ), 1 ) )

  # a larger tree keeps every clone except the removed one
  removed <- quietly( remove.clones.on.tree( tree_example, clones.to.remove = 3 ) )
  survivors <- unique( as.numeric( removed ) )

  expect_false( 3 %in% survivors )
  expect_setequal( survivors, setdiff( unique( as.numeric( tree_example ) ), 3 ) )

  # daughters of the removed clone are reattached to its parent, the root
  reattached <- removed[ removed[, 1 ] == 1, 2 ]
  expect_true( all( c( 5, 6, 7, 8 ) %in% reattached ) )

})

test_that( "remove.clones.on.tree collapses to root -> root when only the root is kept", {

  tr <- matrix( c( 1, 2,
                   2, 4,
                   1, 3 ), ncol = 2, byrow = TRUE )

  collapsed <- quietly( remove.clones.on.tree( tr, clones.to.keep = 1 ) )

  expect_equal( dim( collapsed ), c( 1L, 2L ) )
  expect_equal( as.numeric( collapsed ), c( 1, 1 ) )

})

test_that( "remove.clones.on.tree requires something to remove", {

  tr <- matrix( c( 1, 2, 2, 4, 1, 3 ), ncol = 2, byrow = TRUE )

  expect_error( remove.clones.on.tree( tr ), "not provided info" )

  # Asking to remove a clone that is not on the tree preserves every
  # relationship. Note the rows come back in trunk -> leaf order rather than
  # the input order, because the tree is reordered before pruning, and that
  # only the clones.to.keep route warns here: on the clones.to.remove route
  # the requested clone is absent from the tree but the removal list itself is
  # non-empty, so the no-op warning is not reached.
  unchanged <- quietly( remove.clones.on.tree( tr, clones.to.remove = 99 ) )

  expect_setequal( apply( unchanged, 1, paste, collapse = "->" ),
                   apply( tr, 1, paste, collapse = "->" ) )
  expect_equal( unchanged, quietly( logically.order.tree( tr ) ) )

  # via clones.to.keep, keeping everything is a no-op and does warn
  expect_warning( kept_all <- remove.clones.on.tree( tr, clones.to.keep = c( 1, 2, 3, 4 ) ),
                  "no clones specified to remove" )
  expect_setequal( apply( kept_all, 1, paste, collapse = "->" ),
                   apply( tr, 1, paste, collapse = "->" ) )

})

test_that( "extract_daughters returns all descendants at any depth", {

  # clone 3 has direct daughters 5, 6, 7, 8 and grandchildren 9 (via 5) and
  # 10 (via 6): all six are descendants, and 3 itself is excluded
  descendants <- cloneMap:::extract_daughters( tree_example, 3 )

  expect_setequal( descendants, c( 5, 6, 7, 8, 9, 10 ) )
  expect_false( 3 %in% descendants )

  # depth of one
  expect_setequal( cloneMap:::extract_daughters( tree_example, 5 ), 9 )

  # a leaf clone has no descendants
  expect_length( cloneMap:::extract_daughters( tree_example, 4 ), 0 )

  # from the root, every other clone on the tree is a descendant
  root <- find_root( tree_example )
  all_but_root <- setdiff( unique( as.numeric( tree_example ) ), root )
  expect_setequal( cloneMap:::extract_daughters( tree_example, root ), all_but_root )

  # several parents at once, still excluding the parents themselves
  pair <- cloneMap:::extract_daughters( tree_example, c( 5, 6 ) )
  expect_setequal( pair, c( 9, 10 ) )

})
