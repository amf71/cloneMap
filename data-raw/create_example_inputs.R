### code to prepare example input files for cloneMap() function ###

CCFs_example <- data.frame( clones = 1:11, 
			                   CCF = c(1, 0.4, 0.45, 0.05, 0.3, 0.15, 0.2, 0.1, 0.15, 0.07, 0.15) , stringsAsFactors = F)

CCFs_simple_example <- data.frame( clones = c(1, 2, 4, 5), 
                                   CCF = c(1, 0.4, 0.02, 0.2) , stringsAsFactors = F)

# tree matricies are written with each relaationship as a row, #
# the parent as the first column and child as the second      #

tree_example <- matrix( c(1, 2,
                          1, 3,
                          1, 4,
                          2, 5,
                          3, 7,
                          3, 8,
                          3, 9, 
                          5, 11,
                          5, 6, 
                          7, 10 ), ncol = 2, byrow = TRUE )


clones <- unique( as.numeric(c(tree_example[,1], tree_example[,2]) ) )
getPalette <- colorRampPalette( RColorBrewer::brewer.pal( 8, "Dark2") )
clone_colours_example <- getPalette( length( clones ) )

# OR #

clone_colours_example <- c( "#1B9E77", "#A07125", "#B16548", "#8068AE", "#D03792",
                            "#A66753", "#7FA718", "#D9AA04", "#BF8B12", "#927132",
                            "#666666" ) 

names(clone_colours_example) <- CCFs_example$clones


usethis::use_data( CCFs_example, overwrite = TRUE )
usethis::use_data( CCFs_simple_example, overwrite = TRUE )
usethis::use_data( tree_example, overwrite = TRUE )
usethis::use_data( clone_colours_example, overwrite = TRUE)

### END ###

