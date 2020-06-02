
# cloneMap


## Description

A function to map the distribution of somatic clones in sample or set of samples, accounting for the clone size(s), *ie* the Cancer Cell Fraction (CCF), and phylogenetic relationships between clones. Clone positions are semi-randomised in the plot while maintaining the two formally described factors.  


## Installation

You can use devtools::install_github() to install cloneMap from this repository:

`devtools::install_github("amf71/cloneMap")`


## Usage examples

load package:

```R
library(cloneMap)
```

### Example data

`tree_example_1`, `tree_example_2`, `tree_example`, `CCFs_example_1`, `CCFs_example_2`, and `clone_colours_example` are loaded into the R environment (hidden) upon package loading and also are defined below:


Example CCF tables, these could be from the same tumour

```R
CCFs_example_1 <- data.frame( clones = c( 1,   2,   3,   4 ), 
                              CCF    = c( 1, 0.4, 0.2, 0.1 ), 
                              stringsAsFactors = F)

CCFs_example_2 <- data.frame( clones = c( 1,   2,   3,   5,    6,    7,    8,   9,   10 ), 
                              CCF    = c( 1, 0.1, 0.7, 0.2, 0.25, 0.03, 0.06, 0.1, 0.05 ), 
                              stringsAsFactors = F )
```

Example tree matricies are written with each relationship as a row, the parent as the first column and child as the second    


```R
tree_example_1  <-  matrix( c(1, 2,
                              1, 3,
                              2, 4 ), ncol = 2, byrow = TRUE )

tree_example_2  <-  matrix( c(1, 2,
                              1, 3,
                              3, 5,
                              3, 6,
                              3, 7,
                              3, 8,
                              5, 9, 
                              6, 10 ), ncol = 2, byrow = TRUE )

tree_example   <-   matrix( c(1, 2,
                              1, 3,
                              2, 4,
                              3, 5,
                              3, 6,
                              3, 7,
                              3, 8,
                              5, 9, 
                              6, 10 ), ncol = 2, byrow = TRUE )
``` 


Example colours are define using a named vector of hexidecimal colours with clones as names


```R
clone.names <- unique( c( tree_example[,1], tree_example[,2] ) )

clone_colours_example <- c( "#B15928", "#DDD399", "#9471B4", "#ED8F47", "#FDB762", 
                            "#E52829", "#B89B74", "#79C360", "#3F8EAA", "#A6CEE3" ) 

names(clone_colours_example) <- clone.names
```

### Plot examples


Simple map:

```R
cloneMap( tree_example_1, CCFs_example_1 )
```

![example1](data-raw/example_outputs/example_1.png)

More complex map:

```R
cloneMap( tree_example_2, CCFs_example_2 )
```

![example2](data-raw/example_outputs/example_2.png)


Use a clone_map object to  plot cloneMaps reproducably and much faster:
 
```R
clone_map_eg <- cloneMap( tree_example_2, CCFs_example_2, output.Clone.map.obj = TRUE, plot.data = FALSE )
cloneMap( clone_map = clone_map_eg )
```

![example3](data-raw/example_outputs/example_3.png)


Specify the same clone colours accross several plots:
 
```R
cloneMap( tree_example, CCFs_example_1, clone.cols = clone_colours_example )
cloneMap( tree_example, CCFs_example_2, clone.cols = clone_colours_example )
```

![example4](data-raw/example_outputs/example_4.png)

