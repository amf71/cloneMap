### code to prepare example input files for cloneMap() function ###

CCFs_example <- data.frame( clones = 1:10, 
			                   CCF = c(1, 0.4, 0.1, 0.05, 0.8, 0.3, 0.25, 0.1, 0.15, 0.17) , stringsAsFactors = F)

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


usethis::use_data(DATASET, overwrite = TRUE)

### END ###

