### code to prepare example input files for cloneMap() function ###

CCFs_example <- data.frame( clones = 1:10, 
			                   CCF = c(1, 0.4, 0.1, 0.05, 0.8, 0.3, 0.25, 0.1, 0.15, 0.17) , stringsAsFactors = F)

# tree matricies are written with each relaationship as a row, #
# the parent as the first column and child as th second        #

tree_example <- matrix( c(1, 5,
                          5, 2,
                          1, 3,
                          1, 4,
                          2, 7,
                          5, 6, 
                          6, 7,
                          6, 8,
                          2, 9,
                          2, 10 ), ncol = 2, byrow = TRUE )


clone_colours_example <- c( "#1B9E77", "#AE6D1C", "#A16864", "#9B58A5", "#D8367D",
                            "#749829", "#BBA90B", "#C9930D", "#97722D", "#666666" ) 

names(clone_colours_example) <- CCFs_example$clones


usethis::use_data( CCFs_example, overwrite = TRUE )
usethis::use_data( tree_example, overwrite = TRUE )
usethis::use_data( clone_colours_example, overwrite = TRUE)

### END ###

