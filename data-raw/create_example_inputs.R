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


# example tree and CCF table for normal tissue data - polyclonal no trunk
tree_example_poly <-   matrix( c(1, 1,
                                 2, 2,
                                 3, 4,
                                 3, 5,
                                 4, 6,
                                 7, 8,
                                 9, 9,
                                 10, 10,
                                 11, 11,
                                 12, 12,
                                 13, 13,
                                 14, 14), ncol = 2, byrow = TRUE )

CCF_example_poly <- data.frame( clones = c( 1,   2,    3,   4,    5,    6,    7,    8,    9,   10,     11,   12,    13,   14 ), 
                                CCF    = c( 0.03, 0.05, 0.2, 0.1, 0.02, 0.05,  0.1,  0.05, 0.1, 0.05, 0.02, 0.02, 0.05, 0.03 ), 
                                stringsAsFactors = F )



# get all the clone names we want to plot #
clone.names <- unique( c( tree_example[,1], tree_example[,2] ) )


# can choose colour manually #
clone_colours_example <- c( "#B15928", "#DDD399", "#9471B4", "#ED8F47", "#FDB762", 
                            "#E52829", "#B89B74", "#79C360", "#3F8EAA", "#A6CEE3" ) 

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
usethis::use_data( tree_example_poly, overwrite = TRUE )
usethis::use_data( CCF_example_poly, overwrite = TRUE)

# # clear env and install package again to update example data #
# 
# rm( ls() )
# 
# devtools::install()



## now plot the examples which we show in the README ##

library(cloneMap)


# first simple example #

png( "data-raw/example_outputs/example_1.png" )

cloneMap::cloneMap( tree_example_1, CCFs_example_1, border.thickness = 3 )

invisible( dev.off() )


# second more complex example #

png( "data-raw/example_outputs/example_2.png")

cloneMap::cloneMap( tree_example_2, CCFs_example_2, border.thickness = 3 )

invisible( dev.off() )


# second example again using clone_map object #

clone_map_eg <- cloneMap::cloneMap( tree_example_2, CCFs_example_2, output.Clone.map.obj = TRUE, plot.data = FALSE )

png( "data-raw/example_outputs/example_3.png" )

cloneMap::cloneMap( clone_map = clone_map_eg, border.thickness = 3 )

invisible( dev.off() )


# plot both egs as if they were from the same tumours to show how to specify same colours for clones #

layout <- matrix( c( 1, 2,
                     3, 4 ), ncol = 2, byrow = TRUE )

png( "data-raw/example_outputs/example_4.png", width = 800 )

layout( layout,
        heights = c(1, 5),
        widths = c(2, 2))

par( mar = c(1, 1, 1, 1), xpd = NA)

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0.5, 0.5, labels = "Tumour 1: Sample 1", cex = 3, font = 2)

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0.5, 0.5, labels = "Tumour 1: Sample 2", cex = 3, font = 2)

cloneMap::cloneMap( tree_example_1, CCFs_example_1, clone.cols = clone_colours_example, border.thickness = 3 )

cloneMap::cloneMap( tree_example_2, CCFs_example_2, clone.cols = clone_colours_example, border.thickness = 3 )

invisible( dev.off() )



# plot maps of polyclonal data similar to that found in normal tissues

png( "data-raw/example_outputs/example_polyclonal.png", width = 800 )

cloneMap::cloneMap( tree.mat = tree_example_poly, 
                    CCF.data = CCF_example_poly )

invisible( dev.off() )

# plot maps of polyclonal data similar to that found in normal tissues with border

png( "data-raw/example_outputs/example_polyclonal_border.png", width = 800 )

cloneMap::cloneMap( tree.mat = tree_example_poly, 
                    CCF.data = CCF_example_poly,
                    tissue_border = TRUE)

invisible( dev.off() )

# plot maps of polyclonal data similar to that found in normal tissues with border

png( "data-raw/example_outputs/example_polyclonal_spaced.png", width = 800 )

cloneMap::cloneMap( tree.mat = tree_example_poly, 
                    CCF.data = CCF_example_poly,
                    tissue_border = TRUE,
                    space_fraction = 0.7 )

invisible( dev.off() )


### END ###

