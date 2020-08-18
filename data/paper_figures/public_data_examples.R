####============================================================####
####                                                            ####
#### Script to extract public data examples and apply cloneMaps ####
####                                                            ####
####============================================================####

# Author : Alexander M Frankell
# Date : 10-Aug-2020


##=========================##
## load required libraries ##
##=========================##

library(readxl)
library(httr)
library(data.table)
library(cloneMap)


##======================================================================##
## download the example data from the web and format to expected input ##
##======================================================================##

# Data from McPherson et al 2016 Nature genetics Sorab Shah lab on Ovarian cancer #

data_loc <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fng.3573/MediaObjects/41588_2016_BFng3573_MOESM153_ESM.xlsx'

httr::GET(data_loc, httr::write_disk(tf <- tempfile(fileext = ".xlsx"))) # takes a few minutes
ovarian_ccfs <- data.table::as.data.table( readxl::read_excel( tf , sheet = "S9" ) )

# ovarian trees not present in supplementary but in figure 4 - copied from figures below:
# col 1 = parent and col 2 = child
ovarian_trees <- list( "1" = matrix( c("A", "B",
                                       "B", "C",
                                       "A", "D",
                                       "D", "E",
                                       "D", "F",
                                       "F", "G",
                                       "F", "H", 
                                       "H", "I"), byrow = T, ncol = 2 ),
                       
                       "2" = matrix( c("A", "B",
                                       "A", "C",
                                       "C", "D",
                                       "C", "E",
                                       "E", "F"), byrow = T, ncol = 2 ),
                       
                       "9" = matrix( c("A", "B",
                                       "A", "C"), byrow = T, ncol = 2 ),
                       
                       "3" = matrix( c("A", "B",
                                       "A", "C",
                                       "C", "D",
                                       "D", "E",
                                       "D", "F",
                                       "F", "G"), byrow = T, ncol = 2 ),
                       
                       "4" = matrix( c("A", "B",
                                       "A", "C",
                                       "C", "D",
                                       "C", "E",
                                       "E", "F",
                                       "E", "G",
                                       "E", "H"), byrow = T, ncol = 2 ),
                       
                       "7" = matrix( c("A", "B",
                                       "A", "C",
                                       "C", "D",
                                       "C", "E"), byrow = T, ncol = 2 ),
                       
                       "10" = matrix( c("A", "B",
                                        "A", "C",
                                        "C", "D",
                                        "C", "E",
                                        "E", "F"), byrow = T, ncol = 2 ) )


# 'prevalence' in ovarian data is not CCF but the CCF minus all daughter CCFs     #
# ie what % of cells actually contain only this genotype (and no other mutations) #
# convert back to CCFs using the tree 

ovarian_ccfs[, ccf := sapply( 1:.N, function( i ){
  # get tree for this case
  pat <- ovarian_ccfs[i, patient_id ] 
  tree <- ovarian_trees[ names(ovarian_trees) == pat ][[1]]
  
  # get all the daughters recursively for this clone
  clone <- ovarian_ccfs[i, clone_id ]
  
  daughters <- clone

  repeat{
    
    daughters.new <- unique( c( daughters, tree[ tree[,1] %in% daughters, 2] ) )
    if( all(daughters.new %in% daughters) ) break else daughters <- daughters.new
    
  }
  
  # add the CCF of all of these clones in this sample
  sample <- ovarian_ccfs[i, paper_id ]
  return( ovarian_ccfs[ patient_id == pat & paper_id == sample & clone_id %in% daughters, sum( prevalence ) ] )
  
})]



# Data from TRACERx 100 publication in NEJM : Jamal-Hanjani et al 2017 on lung cancer #

data_loc <- 'https://www.nejm.org/doi/suppl/10.1056/NEJMoa1616288/suppl_file/nejmoa1616288_appendix_2.xlsx'

httr::GET(data_loc, httr::write_disk(tf <- tempfile(fileext = ".xlsx")))
lung_ccfs <- data.table::as.data.table( readxl::read_excel( tf , sheet = "TableS3", skip = 19 ) )
lung_trees <- data.table::as.data.table( readxl::read_excel( tf , sheet = "TableS7", skip = 1 ) )
lung_clin <- data.table::as.data.table( readxl::read_excel( tf , sheet = "TableS2", skip = 1 ) )

## extract the CCFs per clone per sample from the lung data (currently per variant) ##
# remove cases / mutations without PyClone data
lung_ccfs <- lung_ccfs[ !PyClonePhyloCluster == "NA" ]

#split into a list with each element for as 1 patient
lung_ccfs <- lapply( lung_ccfs[, unique(SampleID) ], function( sample ) lung_ccfs[ SampleID %in% sample ] )

# split out the regions, melt into a region column, take an average ccf of all mutations for each clone and bind together 
lung_ccfs <- do.call( rbind, lapply( lung_ccfs, function( sample_table ){
  
  ccfs <- strsplit(sample_table$PyClonePhyloCCF, split = ";" ) 
  regions <- gsub( ":.*$", "" , ccfs[ sapply( ccfs, function(x) !all( x == "NA" ) ) ][[1]] )
  suppressWarnings( ccfs <- lapply( ccfs, function(x) as.numeric( gsub( "^.*:", "", x ) ) ) )
  ccfs <- data.table::as.data.table( do.call( rbind, ccfs ) )
  names(ccfs) <- regions
  
  out <- cbind( sample_table[, .(SampleID, PyClonePhyloCluster )], ccfs ) 
  region_cols <- names(out)[ 3:ncol(out) ]
  out <- data.table::melt( out, id.vars = c("SampleID", "PyClonePhyloCluster"),
                           measure.vars = region_cols,
                           variable.name = "sample",
                           value.name = "ccf" )
  out[, .(ccf = round( mean( ccf ), 3 ) ), by = .(SampleID, PyClonePhyloCluster, sample) ]
  
}))


## now extract the trees in the expected format ##

samples <- lung_trees$SampleID

lung_trees <- lapply( lung_trees$SampleID, function(sample){
  
  tree <- lung_trees[ SampleID == sample, PrimaryTreeStructure ]
  tree <- strsplit(tree, split = ";")[[1]]
  tree <- lapply( tree, function(x) strsplit(x, split = "->" )[[1]] )
  do.call( rbind, tree )
  
})
names(lung_trees) <- samples


### All ccf data at region level, also calculate ccfs over all samples in a patient ###
 
lung_ccfs[, patient_ccfs := mean( ccf ), by = .( SampleID, PyClonePhyloCluster )]
ovarian_ccfs[, patient_ccf := mean( ccf ), by = .(patient_id, clone_id)]

# some CCFs in ovarian data far lower than sensitivity, adjust to 0
ovarian_ccfs[ patient_ccf < 0.01, patient_ccf := 0]
ovarian_ccfs[ ccf < 0.01, ccf := 0]


##==================================##
## Plot clone maps at patient level ##
##==================================##

# get all the cloneMap objects for plotting #

ovarian_patient_maps <- lapply( unique(ovarian_ccfs$patient_id), function( pat ){
  
  # get ccfs with correct input names
  ccfs <- ovarian_ccfs[ patient_id == pat ][ ! duplicated( clone_id ), .(clone_id, patient_ccf) ]
  names(ccfs) <- c("clones", "CCF")
  
  #get tree
  tree <- ovarian_trees[ names(ovarian_trees) == pat ][[1]]
  
  return(
    
    cloneMap( tree.mat = tree, 
              CCF.data = ccfs, 
              output.Clone.map.obj = TRUE, 
              high_qualty_mode = TRUE, 
              plot.data = FALSE )
    
  )
  
})

names(ovarian_patient_maps) <- unique( ovarian_ccfs$patient_id )


ovarian_region_maps <- lapply( unique( paste( ovarian_ccfs$paper_id, ovarian_ccfs$patient_id) ) , function( sample_pat ){
  
  # get ccfs with correct input names
  ccfs <- ovarian_ccfs[ paste( paper_id, patient_id) == sample_pat ][ ! duplicated( clone_id ), .(clone_id, ccf) ]
  names(ccfs) <- c("clones", "CCF")
  
  #get tree
  tree <- ovarian_trees[ names(ovarian_trees) == strsplit(sample_pat, split = " ")[[1]][2] ][[1]]
  
  return(
    
    cloneMap::cloneMap( tree.mat = tree, 
                        CCF.data = ccfs, 
                        output.Clone.map.obj = TRUE, 
                        high_qualty_mode = TRUE, 
                        plot.data = FALSE )
    
  )
  
})

names(ovarian_region_maps) <- unique( paste( ovarian_ccfs$paper_id, ovarian_ccfs$patient_id) )

# lung


lung_patient_maps <- lapply( unique(lung_ccfs$SampleID), function( pat ){
  
  #a few tumour ahve no tree
  if( !any(names(lung_trees) == pat) ) return( NA )
  
  # get ccfs with correct input names
  ccfs <- lung_ccfs[ SampleID == pat ][ ! duplicated( PyClonePhyloCluster ), .(PyClonePhyloCluster, patient_ccfs) ]
  names(ccfs) <- c("clones", "CCF")
  
  #get tree
  tree <- lung_trees[ names(lung_trees) == pat ][[1]]
  
  return(
    
    cloneMap( tree.mat = tree, 
              CCF.data = ccfs, 
              output.Clone.map.obj = TRUE, 
              high_qualty_mode = TRUE, 
              plot.data = FALSE )
    
  )
  
})

names(lung_patient_maps) <- unique( lung_ccfs$SampleID )



#=============================================#
# save the save underlying data for the plots #
#=============================================#

# ccf data

save(lung_ccfs, file = 'data/paper_figures/data/TRACERX_lung_ccfs.rda')
save(ovarian_ccfs, file = 'data/paper_figures/data/Shah_ovarian_ccfs.rda')

# tree data

save(lung_trees, file = 'data/paper_figures/data/TRACERX_lung_trees.rda')
save(ovarian_trees, file = 'data/paper_figures/data/Shah_ovarian_trees.rda')


# clinical data for lung

save(lung_clin, file = 'data/paper_figures/data/TRACERX_lung_clinical_data.rda')

# maps

save(lung_patient_maps, file = 'data/paper_figures/data/TRACERX_lung_tumour_maps.rda')
save(ovarian_patient_maps, file = 'data/paper_figures/data/Shah_ovarian_tumour_maps.rda')
save(ovarian_region_maps, file = 'data/paper_figures/data/Shah_ovarian_region_maps.rda')

##================================##
## plot two clone maps to compare ##
##================================##

pdf( "data/paper_figures/output_plots/complex_eg_CRU0079.pdf" )

cloneMap::cloneMap( clone_map = lung_patient_maps[[ which(names(lung_patient_maps) == 'CRUK0079') ]] )

invisible( dev.off() )

pdf( "data/paper_figures/output_plots/complex_eg_CRU0094.pdf" )

cloneMap::cloneMap( clone_map = lung_patient_maps[[ which(names(lung_patient_maps) == 'CRUK0094') ]] )

invisible( dev.off() )

##========================================##
## function to plot groups of tumour maps ##
##========================================##

plot_grouped_tumour_maps <- function(tumour_group_data, group_name, sample_name, rasters.list, 
                                     sample_identifier, byRegion = FALSE, group_order = NA, track = TRUE, 
                                     no.cols = 20 ){
  
  groups <- tumour_group_data[ order(get(group_name)), .N, by = get(group_name)]
  setnames(groups, c("group","N"))
  
  if( !all(is.na(group_order)) ) groups <-  groups[ match( group_order, group), ]
  
  groups[, rounded_n := mceiling(N, no.cols) ][,
                                               nrow := rounded_n / no.cols ][,
                                                                             cummul_rounded_n := sapply( 1:nrow(groups), function(i)  sum(rounded_n[1:i])) ]
  
  ## now make the layout for the plot ##
  
  layout <- do.call(rbind, lapply( 1:nrow(groups), function(i){
    
    if(i == 1) cum_before <- 2 else cum_before <- groups[i-1,cummul_rounded_n]+ 1 + i
    
    cum_after <- groups[i,cummul_rounded_n] + i
    
    out <- matrix( cum_before:cum_after, nrow = groups[i,nrow], byrow = TRUE)
    
    out <- rbind( matrix(rep(cum_before-1, no.cols * 2), nrow = 2), out)
    
    return( out )
    
  }))
  
  ## now assign each tumour a piont in the layout ##
  
  positions <- data.table( pos = 1:max(layout),
                           to_plot = unlist(lapply(groups[,group], function( group ){
                             out <- group
                             out <- c(out, tumour_group_data[ get(group_name) == group, get(sample_name) ])
                             igroup <- which(groups$group == group)
                             out <- c(out, rep(NA, groups[igroup, rounded_n] - length(out) + 1))
                             return(out)
                           })) )
  
  
  ## now plot all the tumour maps ##
  
  
  #pdf(file = paste0(outputs.folder,"/",date,"_421_tumour_maps_by_Tumour_by_Stage.pdf"), width = 20, height = 30)
  
  layout(layout)
  
  for( i in 1:nrow(positions) ){
    
    is.sample <- grepl( sample_identifier, positions[i, to_plot] )
    
    par( mai = c(0, 0, 0, 0), xpd = NA)
    
    if( is.sample == TRUE ){
      
      tumour <- positions[i, to_plot]
      
      if(track == TRUE) cat( paste0( which(positions[ grepl(sample_identifier, to_plot), to_plot] == tumour ), " " ) )
      
      cloneMap::cloneMap( clone_map = rasters.list[[ which( names(rasters.list) == tumour ) ]] )
      
    } else {
      
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      
      if( !is.na( positions[i, to_plot] ))  text(x = 0.5, y = 0.5, positions[i, to_plot], cex =8, col = "black")
      
    }
    
  }
  
  #invisible( dev.off() )
  
}


standard.layout <- matrix( 1:9 , byrow = TRUE, ncol = 3)

# plot some ovarian data by region

patient <- " 3$"
regions <- names(ovarian_region_maps)[ grepl( patient, names(ovarian_region_maps) ) ]

clone_colours <- cloneMap::make_clone_col_input( unique( ovarian_ccfs[ patient_id == "3", clone_id ] ) )

pdf( "data/paper_figures/output_plots/Ovarian_3_regions.pdf" )

layout( standard.layout )

for(region in regions)  cloneMap::cloneMap( clone_map = ovarian_region_maps[[ which( names(ovarian_region_maps) == region ) ]], clone.cols = clone_colours )

invisible( dev.off() )

patient <- " 10$"
regions <- names(ovarian_region_maps)[ grepl( patient, names(ovarian_region_maps) ) ]

clone_colours <- make_clone_col_input( unique( ovarian_ccfs[ patient_id == "3", clone_id ] ) )

pdf( "data/paper_figures/output_plots/Ovarian_10_regions.pdf" )

layout( standard.layout )

for(region in regions)  cloneMap::cloneMap( clone_map = ovarian_region_maps[[ which( names(ovarian_region_maps) == region ) ]], clone.cols = clone_colours )

invisible( dev.off() )


# plot all ovarian tumours

pdf( "data/paper_figures/output_plots/Ovarian_tumours.pdf" )

layout( standard.layout )

for(tumour in names(ovarian_patient_maps) )  cloneMap::cloneMap( clone_map = ovarian_patient_maps[[ which( names(ovarian_patient_maps) == tumour ) ]] )

invisible( dev.off() )


# plot all lung tumours by stage

lung_clin <- as.data.table( lung_clin )

pdf( "data/paper_figures/output_plots/TRACERx100.pdf" )

plot_grouped_tumour_maps( tumour_group_data =  lung_clin, 
                          group_name = "Stage", 
                          sample_name = 'TRACERxID',
                          sample_identifier = "CRUK", 
                          rasters.list = lung_patient_maps, 
                          no.cols = 10 )

invisible( dev.off() )



