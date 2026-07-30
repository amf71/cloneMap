# make.CCFs.tree.consistant() corrects CCF tables in which a parent clone is
# smaller than the sum of its daughters, which is phylogenetically impossible.
# It warns by design whenever it makes such a correction.

# a deliberately inconsistent table: clone 2 has CCF 0.3 but its daughters
# 3 and 4 total 0.5
inconsistent_tree <- matrix( c( 1, 2,
                                2, 3,
                                2, 4 ), ncol = 2, byrow = TRUE )

inconsistent_CCFs <- data.frame( clones = c( 1, 2, 3, 4 ),
                                 CCF    = c( 1, 0.3, 0.25, 0.25 ),
                                 stringsAsFactors = FALSE )

test_that( "the test fixture really is inconsistent to begin with", {

  # guards the tests below: if this ever became consistent they would pass
  # trivially without exercising any correction
  expect_lt( min_parent_daughter_slack( inconsistent_tree, inconsistent_CCFs ), 0 )

})

test_that( "decreasing daughters removes every parent/daughter inconsistency", {

  corrected <- suppressWarnings(
    make.CCFs.tree.consistant( inconsistent_tree, inconsistent_CCFs ) )

  expect_s3_class( corrected, "data.frame" )
  expect_identical( names( corrected ), names( inconsistent_CCFs ) )
  expect_setequal( corrected$clones, inconsistent_CCFs$clones )

  # the contract: no parent is left below the sum of its daughters
  expect_gte( min_parent_daughter_slack( inconsistent_tree, corrected ), 0 )

  # daughters were scaled down proportionally, and their parent left alone
  expect_equal( corrected$CCF[ corrected$clones == 2 ], 0.3 )
  expect_equal( sum( corrected$CCF[ corrected$clones %in% c( 3, 4 ) ] ), 0.3 )
  expect_equal( corrected$CCF[ corrected$clones == 3 ],
                corrected$CCF[ corrected$clones == 4 ] )

  # the root of a rooted tree is rescaled back to a CCF of 1
  expect_equal( corrected$CCF[ corrected$clones == 1 ], 1 )

})

test_that( "increasing parents also removes every inconsistency", {

  corrected <- suppressWarnings(
    make.CCFs.tree.consistant( inconsistent_tree, inconsistent_CCFs,
                               increase.parents = TRUE ) )

  expect_gte( min_parent_daughter_slack( inconsistent_tree, corrected ), 0 )

  # the parent was raised to meet its daughters, who keep their CCFs
  expect_equal( corrected$CCF[ corrected$clones == 2 ], 0.5 )
  expect_equal( corrected$CCF[ corrected$clones == 3 ], 0.25 )
  expect_equal( corrected$CCF[ corrected$clones == 4 ], 0.25 )
  expect_equal( corrected$CCF[ corrected$clones == 1 ], 1 )

})

test_that( "a correction is reported as a warning", {

  expect_warning( make.CCFs.tree.consistant( inconsistent_tree, inconsistent_CCFs ),
                  "daughters have total CCF" )

})

test_that( "an already consistent table is left alone and warns nothing", {

  consistent <- data.frame( clones = c( 1, 2, 3, 4 ),
                            CCF    = c( 1, 0.6, 0.25, 0.25 ),
                            stringsAsFactors = FALSE )

  expect_gte( min_parent_daughter_slack( inconsistent_tree, consistent ), 0 )

  corrected <- expect_silent( make.CCFs.tree.consistant( inconsistent_tree, consistent ) )

  expect_equal( corrected$CCF, consistent$CCF )

})

test_that( "at least one correction direction must be requested", {

  expect_error(
    make.CCFs.tree.consistant( inconsistent_tree, inconsistent_CCFs,
                               decrease.daughters = FALSE, increase.parents = FALSE ),
    "cannot correct tree" )

})

test_that( "clone_names does not change the returned CCFs", {

  # clone_names is a lookup table (columns internal/orig) used only to report
  # original clone names in warning messages, so results must be unaffected
  name_lookup <- data.frame( internal = c( 1, 2, 3, 4 ),
                             orig     = c( "A", "B", "C", "D" ),
                             stringsAsFactors = FALSE )

  without <- suppressWarnings(
    make.CCFs.tree.consistant( inconsistent_tree, inconsistent_CCFs ) )
  with_names <- suppressWarnings(
    make.CCFs.tree.consistant( inconsistent_tree, inconsistent_CCFs,
                               clone_names = name_lookup ) )

  expect_identical( without, with_names )

})

test_that( "the example data passes through consistently", {

  corrected <- suppressWarnings(
    make.CCFs.tree.consistant( tree_example, CCFs_example_2 ) )

  ordered <- quietly( logically.order.tree( tree_example ) )
  in_CCFs <- quietly( remove.clones.on.tree( ordered, clones.to.keep = corrected$clones ) )

  expect_gte( min_parent_daughter_slack( in_CCFs, corrected ), 0 )

})
