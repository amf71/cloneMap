
# cloneMap


## Description

A function to map the distribution of somatic clones in sample or set of samples, accounting for the clone size(s), *ie* the Cancer Cell Fraction (CCF), and phylogenetic relationships between clones. Clone positions are semi-randomised in the plot while maintains the two formally descrbied metrics.  


## Installation

You can use devtools::install_github() to install dndscv from this repository:

`devtools::install_github("amf71/cloneMap")`


## Usage examples

All objects below provided upon package loading

Simple map:

```R
cloneMap( tree_example, CCFs_simple_example )
```

More complex map:

```R
cloneMap( tree_example, CCFs_example )
```

Use a clone_map object to  plot cloneMaps reproducably and much faster
 
```R
clone_map_eg <- cloneMap( tree_example, CCFs_example, output.Clone.map.obj = TRUE, plot.data = FALSE )
cloneMap( clone_map = clone_map_eg )
```

specify the same clone colours accross several plots
 
```R
cloneMap( tree_example, CCFs_example, clone.cols = clone_colours_example )
cloneMap( tree_example, CCFs_simple_example, clone.cols = clone_colours_example )
```
