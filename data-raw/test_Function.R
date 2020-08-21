###### script to test that the function is still working properly after any changes ######

setwd("/Volumes/proj-tracerx-lung/tctProjects/frankella/R_packages/cloneMaps")

# source the latest code #
source('R/cloneMap.R')

# load the example data #

example_files <- system("ls data", intern = TRUE)

for( file in example_files ) load( paste0("data/", file) )

## run the examples and outputt to test folder ##

# first simple example #

png( "data-raw/test_outputs/example_1.png" )

cloneMap( tree_example_1, CCFs_example_1, border.thickness = 3, high_qualty_mode = T )

invisible( dev.off() )


# second more complex example #

png( "data-raw/test_outputs/example_2.png")

cloneMap( tree_example_2, CCFs_example_2, border.thickness = 3 )

invisible( dev.off() )


# second example again using clone_map object #

clone_map_eg <- cloneMap( tree_example_2, CCFs_example_2, output.Clone.map.obj = TRUE, plot.data = FALSE )

png( "data-raw/test_outputs/example_3.png" )

cloneMap( clone_map = clone_map_eg, border.thickness = 3 )

invisible( dev.off() )


# plot both egs as if they were from the same tumours to show how to specify same colours for clones #

layout <- matrix( c( 1, 2,
                     3, 4 ), ncol = 2, byrow = TRUE )

png( "data-raw/test_outputs/example_4.png", width = 800 )

layout( layout,
        heights = c(1, 5),
        widths = c(2, 2))

par( mar = c(1, 1, 1, 1), xpd = NA)

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0.5, 0.5, labels = "Tumour 1: Sample 1", cex = 3, font = 2)

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0.5, 0.5, labels = "Tumour 1: Sample 2", cex = 3, font = 2)

cloneMap( tree_example_1, CCFs_example_1, clone.cols = clone_colours_example, border.thickness = 3 )

cloneMap( tree_example_2, CCFs_example_2, clone.cols = clone_colours_example, border.thickness = 3 )

invisible( dev.off() )


### END ###



# for testing min function line by line

# tree.mat = NA; CCF.data = NA; clone_map = NA; output.Clone.map.obj = FALSE;
# plot.data = TRUE; high_qualty_mode = FALSE; track = NA; brewer.palette = "Paired";
# clone.cols = NA; border.colour = "grey20";  border.thickness = 1.5;
# resolution.index = 100;  smoothing.par = 10; repeat.limit = 4
 


