### code to prepare example input files for cloneMap() function and make example outputs ###

# example CCF tables, these could be from the same tumour #

CCFs_example_1 <- data.frame( clones = c( 1,   2,   3,   4 ), 
                              CCF    = c( 1, 0.4, 0.2, 0.1 ), 
                              stringsAsFactors = F)

CCFs_example_2 <- data.frame( clones = c( 1,   2,   3,   5,    6,    7,    8,   9,   10 ), 
                              CCF    = c( 1, 0.1, 0.7, 0.2, 0.25, 0.03, 0.06, 0.1, 0.05 ), 
                              stringsAsFactors = F )


# tree matricies are written with each relationship as a row, #
# the parent as the first column and child as the second      #

tree_example_1 <-   matrix( c(1, 2,
                              1, 3,
                              2, 4 ), ncol = 2, byrow = TRUE )

tree_example_2 <-   matrix( c(1, 2,
                              1, 3,
                              3, 5,
                              3, 6,
                              3, 7,
                              3, 8,
                              5, 9, 
                              6, 10 ), ncol = 2, byrow = TRUE )

tree_example <-   matrix( c(1, 2,
                            1, 3,
                            2, 4,
                            3, 5,
                            3, 6,
                            3, 7,
                            3, 8,
                            5, 9, 
                            6, 10 ), ncol = 2, byrow = TRUE )


# get all the clone names we want to plot #
clone.names <- unique( c( tree_example[,1], tree_example[,2] ) )


# can choose colour manually #
clone_colours_example <- c( "#1B9E77", "#A07125", "#B16548", "#8068AE", "#D03792",
                            "#A66753", "#7FA718", "#D9AA04", "#BF8B12", "#927132",
                            "#666666" ) 

# OR #

# can use RColorBrewer to get colours #
getPalette <- colorRampPalette( RColorBrewer::brewer.pal( 8, "Set1") )
clone_colours_example <- rev( getPalette( length( clone.names ) ) )


# names of colours vector are the clone names #
names(clone_colours_example) <- clone.names


usethis::use_data( CCFs_example_1, overwrite = TRUE )
usethis::use_data( CCFs_example_2, overwrite = TRUE )
usethis::use_data( tree_example_1, overwrite = TRUE )
usethis::use_data( tree_example_2, overwrite = TRUE )
usethis::use_data( tree_example, overwrite = TRUE )
usethis::use_data( clone_colours_example, overwrite = TRUE)


# install package again to update example data #

devtools::install()



## now plot the examples which we show in the README ##

library(cloneMap)


# first simple example #

pdf( "data-raw/example_outputs/example_1.pdf")

cloneMap( tree_example_1, CCFs_example_1, border.thickness = 3 )

invisible( dev.off() )


# second more complex example #

pdf( "data-raw/example_outputs/example_2.pdf")

cloneMap( tree_example_2, CCFs_example_2, border.thickness = 3 )

invisible( dev.off() )


# second example again using clone_map object #

clone_map_eg <- cloneMap( tree_example_2, CCFs_example_2, output.Clone.map.obj = TRUE, plot.data = FALSE )

pdf( "data-raw/example_outputs/example_3.pdf")

cloneMap( clone_map = clone_map_eg, border.thickness = 3 )

invisible( dev.off() )


# plot both egs as if they were from the same tumours to show how to specify same colours for clones #

layout <- matrix( c( 1, 2,
                     3, 4 ), ncol = 2, byrow = TRUE )

pdf( "data-raw/example_outputs/example_4.pdf", width = 12 )

layout( layout,
        heights = c(1, 5),
        widths = c(2, 2))

par( mar = c(1, 1, 1, 1), xpd = NA)

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0.5, 0.5, labels = "Tumour 1: Sample 1", cex = 4, font = 2)

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0.5, 0.5, labels = "Tumour 1: Sample 2", cex = 4, font = 2)

cloneMap( tree_example_1, CCFs_example_1, clone.cols = clone_colours_example, border.thickness = 3 )

cloneMap( tree_example_2, CCFs_example_2, clone.cols = clone_colours_example, border.thickness = 3 )

invisible( dev.off() )


### END ###

