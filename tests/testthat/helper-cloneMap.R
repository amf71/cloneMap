# Run plotting code on the null device so no Rplots.pdf (or any other file) is
# left behind in the package directory. The device is closed even if `code`
# throws, so a failing expectation cannot leak an open device into later tests.
with_null_device <- function( code ){
  grDevices::pdf( NULL )
  on.exit( grDevices::dev.off(), add = TRUE )
  force( code )
}

# cloneMap() is deliberately chatty (progress messages) and warns by design
# about CCF/tree inconsistencies. Silence both without hiding errors.
quietly <- function( code ) suppressWarnings( suppressMessages( force( code ) ) )

# For each parent -> daughter row, does the parent itself appear as a daughter
# in an earlier row? Root clones have no parent row, so they are exempt.
parents_precede_daughters <- function( tree ){
  roots <- as.character( find_root( tree ) )
  parents <- as.character( tree[, 1] )
  daughters <- as.character( tree[, 2] )
  all( vapply( seq_len( nrow( tree ) ), function( i ){
    if( parents[ i ] %in% roots ) return( TRUE )
    any( daughters[ seq_len( i - 1 ) ] == parents[ i ] )
  }, logical( 1 ) ) )
}

# Smallest CCF gap between any parent and the sum of its daughters. Negative
# means some parent is still smaller than its daughters, i.e. inconsistent.
min_parent_daughter_slack <- function( tree, CCF.data ){
  real_parents <- unique( tree[, 1 ][ tree[, 1 ] != tree[, 2 ] ] )
  if( length( real_parents ) == 0 ) return( Inf )
  min( vapply( real_parents, function( p ){
    daughters <- tree[ tree[, 1 ] == p, 2 ]
    CCF.data$CCF[ CCF.data$clones == p ] -
      sum( CCF.data$CCF[ CCF.data$clones %in% daughters ] )
  }, numeric( 1 ) ) )
}
