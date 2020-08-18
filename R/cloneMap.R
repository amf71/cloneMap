
#==================================================================#
# Function to produce plots showing clonal composotion of a tissue #
# requires phylogenetic tree data for each clone and the Cancer    #
# Cell Fraction (CCF) of each clone                                #
#==================================================================#

#####################
### Main Function ###
#####################


#' Construct clone maps
#'
#' Function to represent somatic clones present in a tissue accounting for clone
#' size, ie the Cancer Cell Fraction (CCF), and for phylogenetic relationships 
#' between clones
#'
#' @param tree.mat A phylogenetic tree matrix with two columns and one row per 
#' phylogenetic relationship (column1 = parent, column2 = child). All clones
#' must be contained in the tree or they will not be plotted (see example data).
#' 
#' @param CCF.data A table indicating CCFs for each clone in the tree (column
#'  names are "clones" and "CCF", see example data).
#' 
#' @param clone_map Instead of providing a tree and CCF table you may provide
#' a `clone_map` object which has previously been outputted by setting 
#' `output.Clone.map.obj` as TRUE. This contains data specifying tree structure
#' and the positions of clones. When inputting a tree and a CCF table, clone 
#' positions are semi-randomly generated but when inputting a `clone_map` postions
#' will be the identical each time as the postions are recorded in the input and 
#' once the `clone_map` has been generated the plotting is far quicker. 
#' 
#' 
#' @param output.Clone.map.obj If TRUE output a clone_map object containing
#' informtion on the positions of clones. This object can be saved and repeatly 
#' provided back to the `clone_map` argument of this function to reproduce precisely 
#' the same plot and once the `clone_map` has been generated the plotting is far 
#' quicker. 
#' 
#' @param plot.data Plot the cloneMap visualisation now. This may not be desirable if
#' `output.Clone.map.obj` is set to TRUE. 
#' 
#' @param high_qualty_mode This activates several mechanisms to ensure high quality plots
#' at the expense of longer running times, particularly so for some input data. These
#' mechanisms attempt to aviod plots where the same clone seperates into several islands,
#' and hence is not continuous with itself, as well as interweaving of clones. If TRUE Track
#' will also be TRUE unless Track is explicited set to FALSE.
#' 
#' @param track Provide extra feedback on progress of function. Default is FALSE in normal mode but TRUE in high quialty mode. 
#' 
#' @param brewer.palette A qualitative colour pallette from the `RColourBrewer` package to use for the
#' clones.  
#' 
#' @param clone.cols A named vector specifying hexidecimal colours for all clones in your 
#' input (see example). This allows you to keep consistent colours for the same clone 
#' accross several plots with different inputs, for instance if you wish to plot cloneMaps
#' for each of several samples from the same tumour. 
#' 
#' @param border.colour The colour of the clone borders.
#' 
#' @param border.thickness The thickness of the clone borders.
#' 
#' @param resolution.index This indicates the resolution of the model representing clonal
#' growth. More precisely it is the diameter of the rasterised matrix in which each piont
#' represents a position where a clone can reside. A higher resultion index will allow more
#' precise positioning of clones but also increase running time. 
#' 
#' @param smoothing.par Smoothing parameter for the clones. Edges of the clones are smoothened
#' for visual appeal. The higher this parameter the more smoothing will occur. 
#' 
#' @param repeat.limit When `high_qualty_mode` is TRUE the function will test each clone for 
#' continuity (ie that there are not multiple islands or patches of the same clone) and if
#' so it will plot this clone and all its sisters again to attempt to aviod this
#' as clone postions are semi-randomly generated. The `repeat.limit` indicates the maximum number of 
#' allowed repeats before the function will output a warning and continue with the plot.
#' These repeats increase run time significantly when high_qualty_mode` is TRUE. After two 
#' repeats the function will use various mechanisms to decrese the probability of non-contiinuous
#' clones even further but at the expense of even longer run times.
#'  
#' 
#' 
#' @return A `clone_map` object will be returned if `output.Clone.map.obj` is TRUE
#' containing information on clonal positions and tree structure. This output can 
#' be inputted using the `clone_map` argument to repeat the same plot many times.
#' 
#' @author 
#' 
#' Alexander M Frankell, Francis Crick institute, University College London, \email{alexander.frankell@@crick.ac.uk}
#' 
#' @examples 
#' # example objects provided in env after loading package #
#' 
#' cloneMap( tree_example, CCFs_example )
#' cloneMap( tree_example, CCFs_simple_example ) 
#'
#' 
#' # Use a clone_map object to  plot cloneMap reproducably #
#' 
#' clone_map_eg <- cloneMap( tree_example, CCFs_example, output.Clone.map.obj = TRUE, plot.data = FALSE )
#' cloneMap( clone_map = clone_map_eg )
#' 
#' 
#' # specify colours #
#' # plot two samples from the same tumour so the clone colours match #
#' 
#' cloneMap( tree_example, CCFs_simple_example, clone.cols = clone_colours_example )
#' cloneMap( tree_example, CCFs_example, clone.cols = clone_colours_example )
#' 
#' @export
cloneMap <- function( tree.mat = NA, CCF.data = NA, clone_map = NA, output.Clone.map.obj = FALSE,
                       plot.data = TRUE, high_qualty_mode = FALSE, track = NA, brewer.palette = "Paired",
                       clone.cols = NA, border.colour = "grey20",  border.thickness = 1.5,
                       resolution.index = 100,  smoothing.par = 10, repeat.limit = 4 ){
  
  # work out whther to track function in detail #
  if( high_qualty_mode & is.na( track )) track <- TRUE
  if( !high_qualty_mode & is.na( track )) track <- FALSE
  
  # how many core do you have access to for parrelellisation? #
  num_cores <- parallel::detectCores() / 2
  # only worth parrellelising if > 10 cores, otheriwise actually slows code! - need to test further #
  if( num_cores < 10 ) num_cores <- 1
  
  # signpost #
  cat( "formatting and cleaning input data...\n" )
  
  ## check we have some input ##
  
  clone_map_data_supplied <- !all( is.na(clone_map) )
  CCF_data_supplied <- !all( is.na(CCF.data) )
  tree_data_supplied <- !all( is.na(tree.mat) )
  
  if( !( clone_map_data_supplied | (CCF_data_supplied & tree_data_supplied) ) ) stop( "Please provide either a phylogenetic tree matrix and a CCF table or a clone_map object")
  
  
  # if tree supplied in raster input then extract from here #
  if( clone_map_data_supplied ){
    
    # check class is correct (ie is it expected output from this function) #
    if( ! class(clone_map) == "Clone map" ) stop( "incorrect raster input" )
    tree.mat <- clone_map$tree
    clone_names <- clone_map$names_match
    CCF.data <- clone_map$CCFs

  } else {
  
  # check that names are correct #
  correct.names <- all( names(CCF.data) == c("clones", "CCF") )
  if( correct.names == FALSE ) stop( "ensure that column names of CCF table are `clones` and `CCF`" )
  
  # check CCF.data is a table of some sort
  if( !any( class( CCF.data ) %in% c("matrix", "data.frame") ) ) stop( "please provide CCF data in tabular format as described" )
  # check CCF.data is got data in
  if( nrow(CCF.data) == 0 ) stop( "please provide data in CCF.data. No rows detected." )
  
  # Need to deal with numeric clone names for the matrix rasterisation #
  # create a conversion table
  orig_names <- unique( c(tree.mat[,1],  tree.mat[,2], CCF.data$clones ) )
  clone_names <- data.frame( orig = orig_names,
                             new = 1:length(orig_names))
  
  #now convert the clone names in the CCF table and tree
  CCF.data$clones <- clone_names[ match( CCF.data$clones, clone_names$orig ), "new" ]
  tree.mat[, 1] <- clone_names[ match( tree.mat[, 1], clone_names$orig ), "new" ]
  tree.mat[, 2] <- clone_names[ match( tree.mat[, 2], clone_names$orig ), "new" ]
  
  # make sure tree class is correct and it looks like a tree #
  tree.mat <- as.matrix( tree.mat )
  if( !ncol(tree.mat) == 2 ) stop( 'tree input should be a table with two columns with 
                                      each row indicating a parent (column 1) and corresponding 
                                      child (column 2' )
  
  }
  
  # get colours for plotting if these are not provided #
  clone_colours_supplied <- !all( is.na(clone.cols) )
  
  if( !clone_colours_supplied ){
    
    # order the tree so the trunk and earl clones are always the same colours across tumours #
    if( nrow(tree.mat) > 1 ) tree.mat <- logically.order.tree(tree.mat)
    
    clones <- unique( as.numeric(c(tree.mat[,1], tree.mat[,2]) ) )
    
    #subset for clones also in CCF table
    if( CCF_data_supplied ) clones <- clones[ clones %in% CCF.data$clones ]

   clone.cols <- make_clone_col_input( clones, brewer.palette )
    
  } else {
    
    # if clone colours provided then make clone names into internal numeric clone names
    clone.cols$clones <- clone_names[ match( clone.cols$clones, clone_names$orig ), "new" ]
    
  }
  
  ## if the raster data already supplied, skip the main steps and just plot the rasterised data ##

  if( !clone_map_data_supplied ){
    
    #####################################################
    ### Ensure all input data is clean and consistent ###
    #####################################################
    
    # ensure CCF.table is correct class #
    CCF.data <- as.data.frame( CCF.data )
    
    # check that names are correct #
    correct.names <- all( names(CCF.data) == c("clones", "CCF") )
    if( correct.names == FALSE ) stop( "ensure that column names of CCF table are `clones` and `CCF`" )
    
    # remove any clones with estimated CCF of 0 in all regions #
    if( any( CCF.data$CCF == 0 ) ) cat( "Some clones in CCF table have 0 CCF. These will not be plotted.")
    CCF.data <- CCF.data[ CCF.data$CCF > 0 ,]
    
    # limit CCF.table to clones on the tree #
    # if this is the case print those which have been removed #
    tree.clones <- unique( as.numeric(tree.mat) )
    
    if( !all( CCF.data$clones %in% tree.clones ) ) cat( paste0("        ", "clone(s) ", paste0(setdiff(CCF.data$clones, tree.clones), collapse = " "), " found in tree but not phylogenetic CCF table. These have been removed\n" ) ) 
    CCF.data <- CCF.data[ CCF.data$clones %in% tree.clones, ]
    
    # limit tree to clones in the CCF.table #
    # if this is the case print those which have been removed #
    if( !all( tree.clones %in% CCF.data$clones ) ){ 
      cat( paste0("        ", "clone(s) ", paste0(setdiff(tree.clones, CCF.data$clones), collapse = " "), " found in phylogenetic tree but not in CCF table. These have been removed\n" ) )
      tree.mat <- remove.clones.on.tree( tree.mat, clones.to.keep = CCF.data$clones ) 
    }
    
    # if tree has > one relationship, order the tree from Trunk -> branches -> leaves #
    if( nrow(tree.mat) > 1 ) tree.mat <- logically.order.tree( tree.mat ) 
    
    ### specify which clone is the root ###
    # as tree as has beeen ordered this been be the parent (ie column 1) in the first relationship #
    root <- tree.mat[1, 1]
    
    # clone CCF expected in percentages not fractions - if looks like fractions (clonal CCF < 2) then x 100 #
    if( CCF.data[ CCF.data$clones == root, "CCF" ] < 2 ) CCF.data$CCF <- CCF.data$CCF * 100
    
    ### In cases where CCF all daughters > CCF parent decrease the daughter CCFs proportionally to match parent ###
    # If this occurs a warning will be outputted with how much smaller the parent CCF is than its children #
    # given noise in CCF calculations mean we can accept some underestimate of parent CCF but if children # 
    # > 130% of parent this should be checked and tree/clones/CCFs may be incorrect                        #
    CCF.data <- make.CCFs.tree.consistant( tree.mat = tree.mat, CCF.data = CCF.data )   ## Function specified below and in FrankellA_functions.R script
    
    
    ######===============================================================######
    ######                                                               ######
    ######   Create rasterised data to indicate positions of each clone  ######
    ######                                                               ######
    ######===============================================================######
    
    # make matrix that specifies each position using the resolution index (argument) specifying number fo columns and rows #
    # the higher the specified resolution index the more precise the plots, but the longer this code will take to run #
    # populate the plot with 0s (indicates no presence of any clone) #
     
    clones_rasterised <- do.call( cbind, lapply( 1:resolution.index, function(i) rep( 0, resolution.index ) ) ) 
    
    ##############################
    ### Add the clonal cluster ###
    ##############################
    
    # clones are present in "patches" within the plot which simulates (v approximately) how clonal growth occurs #
    # Each clone will nucleate at a partially randomized postion within its parent clone and grow until it reaches its specified CCF #
    # sisters within the same parent will compete for space and grow around one another until both reach the appropriate CCF #
    # clones grow in this mathematically simple octagon shape unless they hit barrier (parent edge or sister) #
    
    # first get nucleus coordinate for clonal cluster (center of the plot) #
    nucleus <- c( resolution.index / 2, resolution.index / 2 )
    
    # now calculate distance from nucleus for all pixels #
    # this distance matrix is used to work out how the clone should grow #
    dist.mat <- make.distance.matrix( clones_rasterised, nucleus ) # function specified below
    
    # make clonal area the diameter of the whole plot #
    clonal.area <- sum( dist.mat < resolution.index / 2 )
    
    ## Calculate equivalent 'area' each subclone should consume depending on its CCF ##
    # CCF should be supplied as a percentage - converted to a fraction #
    CCF.data$area <- round( (CCF.data$CCF / 100) * clonal.area )
    
    # make sure clonal is 100% (sometimes ~99.5 how we calculate it) #
    CCF.data[CCF.data$clones == root, "area"] <- clonal.area
    
    # now determine cut off in distance matrix which results in desired amount of area for th clonal cluster #
    possible_cutoffs <- seq( 0.1, mround( max( dist.mat ), 0.1), 0.1)
    distance_cutoff <- min( possible_cutoffs[ sapply(possible_cutoffs, function(cut) sum(dist.mat<cut)) >= clonal.area ] )
    
    #save the blank version of the raster matrix #
    clones_rasterised_blank <- clones_rasterised
    
    # assign clonal clone to its positions
    clones_rasterised[ dist.mat < distance_cutoff ] <- root
    
    # assign all other positions in plot to Inf ( indicates no clones ) #
    dist.mat[ !dist.mat < distance_cutoff ] <- Inf
    
    #  plot.data == true if you want to plot the rasterised data as it is generate, if plot = FLASE   #
    #  you can set output.Rasterised.data == TRUE then instead of plotting the function will output   #
    #  the rasterised clone data to save a plot whenever you like. This is particularly useful as     #
    #  clones are seeded at random so by saving the rasterised data you can ensure the clones don't   #
    #  change positins ech ttime you plot them. You can input the rasterised data & instead of a CCF  #
    #  table with the Clone_map arguemnt                                                            #
    
    if( plot.data == TRUE ){
      
      # signpost #
      cat( "\nplotting outline...\n" )
      
      # ensure raster is class numeric not char #
      clones_rasterised_plot <- apply( clones_rasterised, 1, as.numeric )
      
      # make raster class #
      rasterPlot <- raster::raster( clones_rasterised_plot )
      
      # set up plot extent #
      raster::plot( raster::rasterToPolygons( rasterPlot ), col = NA, border = NA) 
      
      ## plot the clonal clone ##
      plot <- sf::st_as_sf( raster::rasterToPolygons( rasterPlot, function(x){x == root}, dissolve = TRUE) )
      
      # add smoothnig to make it more visually appealling #
      # note: smoothing will make the plotted area slightly underestmiate the true area but by only very little #
      # also has the advantage that you can see parent clones below at smooth edges and it becomes easier too   #
      # see tree relatinoships                                                                                  #
      plot.smooth <- smoothr::smooth(plot, method = "ksmooth", smoothness= smoothing.par)
      
      # specify border thickness & colour in arguments - default is 1.5 & grey #
      raster::plot( plot.smooth, col = clone.cols[ names( clone.cols ) == root ], border = border.colour, lwd = border.thickness, add = TRUE ) 
      
    }
    
    
    ##################################
    ### Add the subclonal clusters ###
    ##################################
    
    # get all parental clones and deal with the daughters of these in turn, starting with the daughters of the clonal cluster #
    # this will be in order of trunk -> branch -> leaf as tree has been ordered as such #
    # col 1 of tree always = parents #
    parents <- unique( tree.mat[, 1 ] )
    
    # check all parents have parents apart froom the root #
    check_parents_parents <- all( parents[ !parents %in% root ] %in% unique(tree.mat[,2]) )
    if( !check_parents_parents ) stop( "error in input tree, some branches are not connected to the rest of the tree")
    
    # when loop around parents, we need to record the distance matrices and cut offs used for each parental clone #
    # start by recording this for the clonal cluster #
    parent.dists <- list( dist.mat )
    names( parent.dists ) <- root
    parent_distance_cutoffs <- c( distance_cutoff )
    names( parent_distance_cutoffs ) <- root
    
    # signpost #
    cat( "\ngenerating subclone distributions...\n" )
    
    # loop around each parent #
    
    for( parent in parents ){
      
      if( track ) cat( paste0( "\n        ", "dealing with subclones of ", parent, "...\n" ) )
      
      # if only got the clonal cluster don't need to continue - already recorded and plotted the clonal cluster #
      if( all(tree.mat == root) ) break
      
      # extract daughters of this clone and how many there are #
      # daughters are always col 2 of tree #
      daughters <- tree.mat[ tree.mat[,1] == parent, 2]
      
      # order daughters by CCF - this helps to aviod non-continuous clones below # 
      daughters <- CCF.data[ order(CCF.data$CCF) & CCF.data$clones %in% daughters, "clones" ]
      
      # from the saved lists,  extract the parent's distance matrix and cut off for the boundry #
      parent.dist.mat <- parent.dists[[ which( names( parent.dists ) == parent ) ]]
      parent_distance_cutoff <- as.numeric(  parent_distance_cutoffs[ which( names( parent_distance_cutoffs ) == parent ) ] )
      
      ######################################################
      ### Grow all daughter subclones around one another ###
      ######################################################
      
      # occasionally clones are not continuous (ie they split into two patches) which is unlikely biologically #
      # if this occurs repeat the nucleation and growing of clones in the hope this doeesns't happen again     #
      # repeat occurs until clones are no longer continuous or reaches a speicified limit (default = 10) when  #
      # a warning is outputted and which the aannoying parent is and the last repeat is plotted                #
      
      # save raster so can start again with each repeat #
      clones_rasterised_parent <- clones_rasterised
      
      #object to record number of repeats in the below loop that have occured #
      repeati <- 0
      
      repeat{
        
        clones.finalised <- FALSE
        repeati <- repeati + 1
        clones_rasterised <- clones_rasterised_parent
        
        ########################################################
        ### determing position of nuclei for daughter clones ###
        ########################################################
        
        # choose nuclei positions near the centre of the parent but also not too near the edge #
        # if muliple clones they must nucleate a % of the parent diameter away from each other #
        # otherwise random #
        
        # if just one subclone nucleate subclone very near the centre of the parent #
        # if > 1 subclone nucleate subclones as far asway as possible from one another #
        # while stll maintainig a certain distance from the edge and centre of the parent clone #
        
        # signpost #
        if( track ) cat( paste0( "        ", "determining nuclei positions...\n" ) )
        
        # determine sum daughter CCF and paarent CCF #
        total.CCF.daughters <- sum( CCF.data$CCF[ CCF.data$clones %in% daughters ] )
        CCF.parent <- CCF.data$CCF[ CCF.data$clones == parent ]
        
        if( length(daughters) > 1 ){
          
          
          # normally allow clones to nucleate nither too fro away or too close to centre of parent
          nucleus.options <- parent.dist.mat > (parent_distance_cutoff * 0.40) & parent.dist.mat < (parent_distance_cutoff * 0.80) & clones_rasterised_parent == parent
          
          # space clones out as mch as possible if >2 and total CCF is high% 
          if( length(daughters) > 2 & total.CCF.daughters > CCF.parent * 0.5 ){
            nucleus.options <- parent.dist.mat > (parent_distance_cutoff * 0.7) & parent.dist.mat < (parent_distance_cutoff * 0.75) & clones_rasterised_parent == parent
          }
          
          # if already failed twice be less restrictive about positioning #
          if( repeati > 2 ){
            nucleus.options <- parent.dist.mat > (parent_distance_cutoff * 0.4) & parent.dist.mat < (parent_distance_cutoff * 0.8) & clones_rasterised_parent == parent
          } 
          
          # convert raster matrix TRUE positions to cordinates for possible nucleation #
          nucleus.options <- matrix.index.to.coordinates( matrix.index = which( nucleus.options ), nrow = nrow( clones_rasterised ), ncol = ncol( clones_rasterised ) ) ## Function specified below
          
          # randomly select 20 sets of n (n = number of clones) nuclei for clones - might be better which more but distance caluclatioons will take longer #
          nuclei_sample_number = 20
          
          # space clones out as mch as possible if >2 and total CCF is high #
          # try to space out even mor if already unsuccesful oon two attempts #
          # if there noly ne clone no need to have multiple nulcie really #
          if( length(daughters) > 2 & total.CCF.daughters > CCF.parent * 0.5 ) nuclei_sample_number = 50
          if( repeati > 2 ) nuclei_sample_number = 200
          if( length(daughters) == 1) nuclei_sample_number = 1
          
          nucleus.options.sel <- lapply( 1:nuclei_sample_number, function(x) nucleus.options[ sample( 1:nrow( nucleus.options ), length(daughters), replace = F), ] )
          
          # for each option work out the average distance between the sets of nuclei #
          # this can take a while so parrellelise #
          nucleus.options.min.dists <-  sapply( nucleus.options.sel, function( nucleus.option ){
            
            mindist <- sapply( 1:nrow(nucleus.option), function(i){
              
              dists <- make.distance.matrix( clones_rasterised_blank, nucleus = as.numeric( nucleus.option[i,  ] ) )
              
              # having tested, its faster to parrellelise here than in outer apply #
              dists <- unlist( parallel::mclapply( which( !1:nrow( nucleus.option ) == i ), function(i2) dists[ nucleus.option[ i2, "x" ], nucleus.option[ i2, "y" ] ], mc.cores = num_cores ))
              
              return( min( dists ) )
            })
            
            return( min( mindist ) )
          } )
          
          # randomly choose the set of nuclei which is in the upper fifth of distances apart #
          
          # space clones out as mch as possible if >2 and total CCF is high% 
          if(!( length(daughters) > 2 & total.CCF.daughters > CCF.parent * 0.5 )){
            
            nucleus.options.sel <- nucleus.options.sel[ nucleus.options.min.dists > quantile( nucleus.options.min.dists, 0.80 ) ]
            selected.i <- sample( 1:length( nucleus.options.sel ), 1 )
            nuclei <- nucleus.options.sel[[ selected.i ]]
            nuclei_min_distance <-  nucleus.options.min.dists[[ selected.i ]]
            
          } else {
            
            max_dist_i <- which( nucleus.options.min.dists == max( nucleus.options.min.dists, na.rm = T ))
            nuclei <- nucleus.options.sel[[ max_dist_i ]]
            nuclei_min_distance <-  nucleus.options.min.dists[[ max_dist_i ]]
            
          }
          
          nuclei <- lapply( 1:nrow(nuclei), function(x) as.numeric( nuclei[x, ] ) )
          names(nuclei) <- daughters # nuclei = the set of nuclei chosen for the daughters of this parent
          
        } else {
          
          # if only one daughter clone just place very near the centre (0.3 * distance to edge) of the parent #
          nucleus.options <- which( parent.dist.mat < (parent_distance_cutoff * 0.3) & clones_rasterised_parent == parent )
          nucleus.options <- nucleus.options[ sample( 1:length( nucleus.options ), 1 ) ]
          nuclei <- matrix.index.to.coordinates( nucleus.options, nrow = nrow( clones_rasterised ), ncol = ncol( clones_rasterised ) )
          nuclei <- lapply(1:nrow(nuclei), function(x) as.numeric(nuclei[x,]))
          names(nuclei) <- daughters
          
        }
        
        ##############################################################
        ### now allow clones to grow from nuclei around each other ###
        ##############################################################
        
        # extract new distance matricies for each of these new nuclei #
        nuclei.dists <- lapply( 1:length( nuclei ), function(i) make.distance.matrix( clones_rasterised, nucleus = nuclei[[i]] ) ) ## make.distance.matrix function specified below
        names( nuclei.dists ) <- daughters
        
        # just check that nuclei are defninitely within the parent (found this error a few times) #
        # if nucleus outside parent then will have been set to Inf therefore the matrix will lack a 0 #
        bad_nucleus <- sapply( sapply( 1:length(daughters), function(i) which( nuclei.dists[[i]] == 0 )), length) == 0
        if( any( bad_nucleus ) ){
          
          stop( paste0( "(BUG) chosen nucleus outside of parent for clone ", daughters[ bad_nucleus ], ". Try rerunning.") )
          
        }
        
        # grow clones regularly to allow clone to grow evenly relative to one another #
        
        # define targe area for each clone #
        clone.areas <- sapply(names(nuclei), function(clone) CCF.data[  CCF.data$clones == clone, "area" ])
        
        # define stage to grow the clones at the same rate - res/100 units of area at a time to the max of all daughter clones #
        
        # determnie how fast to grow - takes longer but more accurate in some nistances #
        rate.modifier <- 1
        
        # grow the clones more slowly if >2 and total CCF is high  #
        if( length(daughters) > 2 & total.CCF.daughters > CCF.parent * 0.5 & high_qualty_mode )  rate.modifier <- 4
        growth.rate <-  resolution.index / rate.modifier
        
        # caluclate max grow for biggest clone #
        max_growth <- mceiling( max(clone.areas), growth.rate)

        #don't strt grow from 0 to save time - grow in one step until the point you might encourter other clones #
        # if only one daughter clone then you don't need to grow them at all - go straight to max size #
        if(  length(daughters) == 1 ) min_growth <- max_growth else {
          
          # shuold be able to grow first to where none of them are touching to speed things up #
          # Use minimum distance between clones nucleo to calculate this #
          min_growth <-  mfloor( sum( nuclei.dists[[1]] < nuclei_min_distance * 0.45 ),  growth.rate ) ## function mfloor is below
          
        }
        
        # if min distance between clones is greater than max clone size growth then just grow in one step #
        if( min_growth > max_growth ) min_growth <- max_growth
        
        growth.stages <- seq( min_growth, max_growth, growth.rate ) 
        
        # if clones are growing make the last 5 rounds slower to avoid inter-weaving of clones #
        if(  min_growth < max_growth ){
          
          #normally do the lst 5 steps more slowly
          number_of_last <- 5
          
          # if < 5 steps do it for all steps #
          if( length(growth.stages) < 5 ) number_of_last <- length(growth.stages)
          
          # get the areas for the last expaanssions #
          last_expansions <- growth.stages[ ( length(growth.stages) - number_of_last ):length( growth.stages ) ]
          
          # remove them #
          growth.stages <- growth.stages[ !growth.stages %in% last_expansions ]
          
          # add the new more frequent grow stages to replace #
          growth.stages <- c( growth.stages, seq( min(last_expansions), max(last_expansions), resolution.index / 4 ))
          
        }
        
        # signpost #
        if( track ) cat( paste0( "        ", "growing clones to match CCF...\n" ) )
        
        # track this with a progressor bar #
        if( track ) pb <- txtProgressBar( min = 0, max = length(growth.stages), style = 3, width =  30 ) 
        
        for( area in growth.stages ){
          
          # for this round of growth determine area of expansion for each clone, limiting to area for this stage #
          areas <- rep( area, length(daughters) )
          names(areas) <- daughters
          area.too.large <- sapply(daughters, function(daughter) areas[ names(areas) == daughter ] > clone.areas[ names(clone.areas) == daughter ] ) 
          areas[ area.too.large ] <- clone.areas[ area.too.large ]
          
          # now determine distaance cut off which results in desired area #
          
          #determine where is availbelf  ro the clones to grow (area that is parne tbut not any of its sisters)
          is_parent_not_sister <- lapply( 1:length( daughters ), function(i) clones_rasterised_parent  == parent & !clones_rasterised %in% daughters[ !daughters %in% daughters[i] ] )
          
          avialible.space <- lapply( 1:length( daughters ), function(i) nuclei.dists[[i]][ is_parent_not_sister[[i]] ] )
          avialible.space <- sapply( 1:length( daughters ), function(i) max( avialible.space[[i]] ))
          distance_cutoffs <- sapply( 1:length(daughters), function(i){
            
            rounded.area <- mround( avialible.space[[i]], 0.1 ) 
            cut.options <- seq( resolution.index / 100 , rounded.area, resolution.index / 100 )
            clone_areas_for_cut_offs <- sapply(cut.options, function(cut) sum( nuclei.dists[[i]] < cut & is_parent_not_sister[[i]] ))
            less_then_target_area <- clone_areas_for_cut_offs <= areas[ names(areas) == daughters[i] ]
            
            # sometimes cannot expand anymore and no option is less_then_target_area hence max(c()) which return -Inf hence clone won't grow - appropriate action #
            cut <- suppressWarnings( max( cut.options[ less_then_target_area ] ) )
            
            return( cut )
          })
          names(distance_cutoffs) <- daughters
          
          # now record the areas on the rasterised plot with the clone names #
          for(i in 1:length(daughters)){
            clones_rasterised[ nuclei.dists[[ i ]] < as.numeric( distance_cutoffs[i] ) & is_parent_not_sister[[i]] ] <- daughters[i]
          }
          
          # update progress bar #
          if( track ) setTxtProgressBar(pb, which( growth.stages == area ) ) 
          
        }
        
        # test that all clones we are trying to plot are present in raster #
        clones.present <- sapply( daughters, function( clone ) any( clones_rasterised == clone ) )
        if( any( !clones.present ) ){
          warning( paste0( "Clone ", daughters[ !clones.present ], " is missing in the raster data, the CCF of this clone is too small for the given resolution of the plot\n this clone will not be shown\n" ) )
          daughters <- daughters[ clones.present ]
        }
        
        # test if all clones are continuous if >2 clones & high CCF - if not repeat #
        # if <= 2 clones its impossible they won't be continuous #
        if( length(daughters) > 1 & total.CCF.daughters > CCF.parent * 0.7 & high_qualty_mode ){
          
          # signpost #
          if( track ) cat( paste0( "\n        ", "checking that clones are continuoous...\n" ) )
          
          are_continuous <- sapply( daughters, function(clone) continuous.test( clone_position = clones_rasterised == clone ) ) ## continuous.test function specified below
          if( all( are_continuous ) ){
            
            clones.finalised <- TRUE
            
          } else  if( track ) cat( paste0( "        ", "clone(s) ", paste( daughters[ !are_continuous ], collapse = " "), " not continuous so repeat...\n" ) )
          
          # only allow certain number of repeats and then give up (default = 10) #
          non_cont <- paste(names( are_continuous )[ !are_continuous ], collapse = ", ")
          if( repeati == repeat.limit ){ cat( paste0( "reached repeat limit to achieve continuous daughter clones (", non_cont, ") in parent clone ", parent, "\n" ) ) ; clones.finalised <- TRUE }
          
        } else  clones.finalised <- TRUE
        
        if( clones.finalised ==  TRUE ){
        
        # recenter the clones depending on how they've grown  #
        nuclei.dists <- lapply( 1:length(daughters), function(i) recenter_distance_matrix( clone_position = clones_rasterised == daughters[i] ) )  # recenter_distance_matrix function below
        names(nuclei.dists) <- daughters
        # extract new cut offs - just use the max distance away from nucleus where each clone is #
        distance_cutoffs <- sapply( 1:length(daughters), function(i) max( nuclei.dists[[i]][ clones_rasterised == daughters[i] ] ) )
        names(distance_cutoffs) <- daughters
        
        ## Once we have the finial distributiions, for each of these clones record thier distrebution ##
        ## and distance cut offsso we can extract this information later if they are parents       ##
        
        for(clone in daughters){
          
          parent.dists <- c( parent.dists, list( nuclei.dists[[ which( names( nuclei.dists ) == clone ) ]] ) )
          names(parent.dists) <- c( names(parent.dists)[ 1:( length(parent.dists) -  1 ) ], clone )
          parent_distance_cutoffs <- c( parent_distance_cutoffs, distance_cutoffs[[ which( names( distance_cutoffs ) == clone ) ]] )
          names(parent_distance_cutoffs) <- c( names( parent_distance_cutoffs )[ 1:( length(parent_distance_cutoffs ) - 1 ) ], clone)
          
        }
        
        break
        
        }
        
      }
      
      ################################################################
      ### Now plot the new clones which we have assigned positions ###
      ################################################################
      
      if( plot.data == TRUE ){
        
        # signpost #
        if( track ) cat( paste0( "\n        ", "plotting subclones...\n" ) )
        
        # loop round daughters of this parent clone and plot ecah of them #
        
        for(clone in daughters){
          
          # ensure raster is numeric #
          clones_rasterised_plot <- apply( clones_rasterised, 1, as.numeric )
          
          # convert to class raster #
          rasterPlot <- raster::raster( clones_rasterised_plot )
          
          plot <- sf::st_as_sf( raster::rasterToPolygons( rasterPlot, function(x){ x == clone }, dissolve = TRUE) )
          plot.smooth <- smoothr::smooth(plot, method = "ksmooth", smoothness = smoothing.par) # smoothing par specified in arguemnts
          raster::plot( plot.smooth, col = clone.cols[ names(clone.cols) == clone], border = "grey20", lwd = border.thickness, add = TRUE ) # border thickness specified in arguemnts and col can be specified in arguments
          
        }
        
      }
      
    }
    
    # if specified in arguments output the raster matrix of clone positions #
    # this enables yu to repeatedly make the exact same plot or otherwise the clone positions will change each time #
    # also once you've made the raster the plotting alone is much faster #
    
    if( output.Clone.map.obj == TRUE ){
      
      clones_rasterised_plot <- apply( clones_rasterised, 1, as.numeric )
      
      clones_rasterised <- list( clone_matrix = clones_rasterised_plot, 
                                 tree = tree.mat, 
                                 names_match = clone_names,
                                 CCFs = CCF.data )
      
      class(clones_rasterised) <- "Clone map"
      
      return( clones_rasterised )
      
    }
    
    
  } else {  ## if clone_map_data_supplied
    
    ######===================================================================================######
    ######                                                                                   ######
    ######   if you have supplied a raster matrix of clone positions then use this to plot   ######
    ######                                                                                   ######
    ######===================================================================================######
    
    #  first object in Clone_map is the rsater matrix #
    clones_rasterised <- clone_map$clone_matrix
    
    # already extracted tree from raster data object #
    
    #extract all cloness to the plot from the tree (some f these may have been overtaken in raster plot) #
    #ensure they arre in correct order - tree will be rodered trrunk -> branches -> leaves
    clones <- unique( c(tree.mat[,1], tree.mat[, 2]) )
    
    # simulate what occurs when plotting & raster generated concurrently #
    # for each clone in the tree (trunk -> leaves) plots the area this occupies indluing all its daughters #
    
    clones_rasterised_blank <- clones_rasterised
    clones_rasterised_blank[] <- 0 
    clones_rasterised_blank <- apply( clones_rasterised_blank, 1, as.numeric )
    blank.plot <- raster::rasterToPolygons( raster::raster( clones_rasterised_blank ) )
    raster::plot( blank.plot, col = NA, border = NA ) # set up plot extent
    
    # signpost #
    cat( paste0( "\n        ", "plotting clones...\n" ) )
    
    for( clone in clones ) {
      
      # signpost #
      if( track ) cat( paste0(  "        ", "dealing with subclones of ", parent, "...\n" ) )
      
      # extract names of all daughters of this parent rercursively down the tree #
      
      clone.daughters <- clone
      
      # look down the tree until you run out of clones #
      
      repeat{ 
        clone.daughters.last <- clone.daughters
        clone.daughters <- unique( c(clone.daughters, tree.mat[ tree.mat[, 1] %in% clone.daughters, 2]) )
        if( all(clone.daughters %in% clone.daughters.last) ) break
      }
      
      # remove parent clone from dughter list #
      
      clone.daughters <- clone.daughters[!clone.daughters == clone]
      
      # make a new raster object where all the daughters are == parent clone #
      
      clones_rasterised_plot <- clones_rasterised
      clones_rasterised_plot[ clones_rasterised_plot %in% clone.daughters ] <- clone
      clones_rasterised_plot <- apply(clones_rasterised_plot, 1, as.numeric)
      
      # convert to class raster #
      rasterPlot <- raster::raster(clones_rasterised_plot)
      
      # plot the clone #
      plot <- sf::st_as_sf( raster::rasterToPolygons( rasterPlot, function(x){x == clone}, dissolve = TRUE))
      plot.smooth <- smoothr::smooth(plot, method = "ksmooth", smoothness= smoothing.par) # smoothing par speicified in arguemnts
      raster::plot(plot.smooth, col = clone.cols[ names(clone.cols) == clone ], border = "grey20", lwd = border.thickness, add = TRUE) # border thickness specified in arguemnts and col can be specified in arguments
      
    }
    
  }
  
  ###########
  ### END ###
  ###########
  
}

#' Function to create clone colour table with a RColourBrewer palette
#' 
#' 
#' @export
make_clone_col_input <- function( clones, brewer.palette = "Paired" ){
  
  if( length(clones) > 8 ){ 
    
    # suppress warning - gets max number of colours from pallette - none have > 12 #
    getPalette <- suppressWarnings( colorRampPalette( RColorBrewer::brewer.pal( 12, brewer.palette) ) ) # brewer.palette specified in arguments, default = "Paired"
    clone.cols <- getPalette( length( clones ) )
    
  } else {
    
    clone.cols <- RColorBrewer::brewer.pal(n = length(clones), name = brewer.palette)
    
  }
  
  names(clone.cols) <- clones
  
  return( clone.cols )
  
}


## function to generate a distnce martix, specifiying how far away each piont in the raster matrix    ##
## is from a central nucleus. So far managed two ways to calculaate this, one that specifies distance ##
## purely in two directions ( up and down ) hence casues a square-ish growth pattern and one that     ##
## accounts for a smaller than 2 distance for diagonal movement which creates an octagonal pattern    ##
## the latter is the default for now, ideally probably like to create a circular pattern but not      ##
## sure how to do this                                                                                ##  

make.distance.matrix <- function( matrix.outline = clones_rasterised, nucleus = nucleus, type = "octoagon" ){
  
  
  ### just count rows & cols from nucleus (makes 'square' shape ) ###
  
  if( type == "square" ){
    
    ## now calculate distance from nucleus for all pxls ##
    
    # dist with horizontal as primary #
    
    dist.mat.h <- do.call( cbind, lapply( 1:ncol( matrix.outline ), function( col ){
      
      # how many columns from nucleus? #
      base.dist <- abs( nucleus[2] - col ) 
      
      return( rep( base.dist, nrow( matrix.outline ) ) )
      
    } ))
    
    # dist with vertical as primary #
    
    dist.mat.v <- do.call( rbind, lapply( 1:nrow( matrix.outline ), function( row ){
      
      # how many rows nucleus? #
      base.dist <- abs( nucleus[1] - row )
      
      return( rep( base.dist, ncol( matrix.outline ) ) )
      
    } ))
    
    # now combine for shortest poss distance #
    dist.mat <- do.call( cbind, lapply( 1:ncol( matrix.outline ), function( col ){
      
      out <- cbind( dist.mat.h[, col ], dist.mat.v[, col ])
      return( as.numeric( rowSums( out ) ) )
      
    }))
    
  }
  
  ### if accounting for diagonal distance type == 'octogon' ###
  
  if( type == "octoagon" ){
    
    ## now calculate distance from nucleus for all pxls ##
    
    # dist with horizontal as primary #
    
    dist.mat.h <- do.call( cbind, lapply( 1:ncol( matrix.outline ), function( col ){
      
      # how many columns from nucleus? #
      base.dist <- abs( nucleus[2] - col )
      
      # addition allow us to account for diagnal movement with change in rows #
      horizontal.additions <- abs( 1:nrow( matrix.outline ) - nucleus[1] ) * 0.41
      
      return( base.dist + horizontal.additions )
      
    } ))
    
    
    # dist with vertical as primary #
    
    dist.mat.v <- do.call( rbind, lapply( 1:nrow( matrix.outline ), function( row ){
      
      # how many rows from nucleus? #
      base.dist <- abs( nucleus[1] - row )
      
      # addition allow us to account for diagnal movement with change in cols #
      vertical.additions <- abs( 1:ncol( matrix.outline ) - nucleus[2] ) * 0.41
      
      return( base.dist + vertical.additions )
      
    } ))
    
    # now combine for shortest poss distance #
    
    dist.mat <- do.call( cbind, lapply( 1:ncol( matrix.outline ), function( col ){
      
      out <- cbind( dist.mat.h[, col ], dist.mat.v[, col ])
      return( as.numeric( qlcMatrix::rowMax( out ) ) )
      
    }))
    
  } 
  
  # return output #
  return( dist.mat )
  
  #########
  ## END ##
  #########
  
}


#' function to order tree root -> branches -> leaves
#'
#'
#'
#'
#' @export
logically.order.tree <- function(tree){
  
  # if only 1 clone on tree then just output it as it is #
  if( all( tree == unique(tree)[1] ) ) return(tree)
  
  ### assign levels to each parent clone depending on how near the trunk ###
  
  # work out the root (the only clone tht's never a daughter ) #
  root <- unique( tree[,1] [! tree[,1] %in% tree[,2] ] )
  
  # check that there is one 'root' this means there is no true root output error #
  if( length(root) > 1 ) stop( "Tree is not rooted. Please provide one clone from which all others arise.")
  
  # check that no clone has two parents (impossible) #
  dups.i <- duplicated( tree[,2] )
  if( any( dups.i ) ) stop( paste0( "Clone(s)", tree[ dups.i, 2 ], "have multiple parents in tree. Not possible." ) )
  
  # make empty list for levels of each clone #
  levels <- rep(NA,nrow(tree))
  
  # aasign the root to level 1 #
  levels[tree[,1] %in% root ] <- 1
  
  # list aall the daughter's of root to work out the level #
  daughters <- tree[ tree[,1] %in% root,2]
  
  # go down the tree in levels and asssign clones correct levels #
  l <- 1
  repeat{
    l <- l + 1
    parents <- daughters
    levels[tree[,1] %in% parents ] <- l
    daughters <- tree[ tree[,1] %in% parents,2]
    if(all(!is.na(levels))) break
  }
  tree <- tree[order(levels),]
  
  if(class(tree)=="character" | class(tree)=="numeric") tree <- matrix(tree, ncol = 2, byrow = TRUE)
  
  return(tree)
  
}



#' function to remove clones from a phenogenetic tree matrix
#'
#' This maintains parent -> daughter relationships,  even if an intermediate
#' (branch) clone is being removed
#'
#'
#' @export
remove.clones.on.tree <- function(tree, clones.to.remove = NA, clones.to.keep = NA){
  
  # check that info have been provided on what clones to remove
  if(all(is.na(clones.to.remove)) & all(is.na(clones.to.keep))){
    stop("you have not provided info on which clones to remove")
  }
  
  # check classes are correct #
  clones.to.keep <- as.character( clones.to.keep )
  tree <- as.matrix( tree )
  
  # define list of all clones in the tree #
  all.clones <- unique(as.character(tree))
  
  # if clones to keep specific rather than clones         ##
  # to remove then work out from this which to be removed ##
  if(all(is.na(clones.to.remove))){
    clones.to.remove <- setdiff(all.clones,clones.to.keep)
  }
  
  # ensure class is correct and all clones specified to keep are in the tree matrix #
  clones.to.remove <- as.character(clones.to.remove)
  clones.to.keep <- clones.to.keep[clones.to.keep %in% all.clones]
  
  # if nothing to remove that is ctually ono the tree then just return the tree with a warning #
  if( length(clones.to.remove) == 0 ){
    warning( "no clones specified to remove that are on the tree\n")
    return(tree)
  }
  
  # order the tree trunk -> branches -> leaves #
  tree <- logically.order.tree(tree)
  
  # get the root clone #
  root <- tree[ 1, 1 ]
  
  # if only the root left then return root as parent and daughter to maintaian the same structure #
  if(all(clones.to.keep==root)){
    return(matrix(c(root,root),ncol = 2, byrow = TRUE))
  }
  
  # now loop round ecah clone to remove, get rid of all relationships its involved in and  #
  # then reassign aany daughter(s) to its parent                                           #
  for(clone in clones.to.remove){
    parent <- tree[tree[,2]==clone,1]
    if(any(tree[,1]==clone)){
      daughters <- tree[tree[,1]==clone,2]
      tree <- rbind(tree, matrix(c(rep(parent,length(daughters)),daughters),ncol = 2))
    }
    tree <- tree[!(tree[,1]==clone | tree[,2]==clone),]
    if(class(tree)=="character" | class(tree)=="numeric") tree <- matrix(tree,ncol = 2,byrow = TRUE)
  }
  
  # return pruned tree #
  return(tree)
  
  # END #
  
}


#' Function to correct CCFs when sum of daughters CCF > parent CCF 
#'
#'
#'
#'
#' @export
make.CCFs.tree.consistant <- function( tree.mat, CCF.data, warning.limit = 1 , parent.adjust = 1,
                                       decrease.daughters = TRUE, increase.parents = FALSE ){
  
  # one of decrease daughters or parents must be true
  if( all( !decrease.daughters & !increase.parents ) ) cat( "Please set either decrease.daughters or decrease.parents arguments to TRUE or cannot correct tree\n" )
  if( increase.parents == TRUE ) decrease.daughters <- FALSE
  
  # order tree trunk -> leaf #
  tree.mat <- logically.order.tree( tree.mat )
  
  # limit to clones in CCF table #
  tree.clones <- unique( as.numeric(tree.mat) )
  if( !all( tree.clones %in% CCF.data$clones) ) tree.mat <- remove.clones.on.tree( tree.mat, clones.to.keep = CCF.data$clones )
  
  # get root #
  root <- tree.mat[ 1, 1 ]
  
  # get ordered list of parents #
  parent.order <- unique(tree.mat[,1])
  if( decrease.daughters == FALSE ) parent.order <- rev( parent.order )
  
  # detect if fractions or percentages #
  is_frac <- CCF.data[ CCF.data$clones ==  root, "CCF" ] < 2
  if( is_frac ) clonal_CCF <- 1 else clonal_CCF <- 100
  
  # for each parent check if it has CCF < sum of daughters if so correct it #
  for(parent in parent.order){
    
    # get names of daughteer clones #
    daughters <- tree.mat[ tree.mat[, 1] == parent, 2 ]
    
    # get row indices for daughters and parent #
    parentrow <- which( CCF.data$clones == parent )
    daughterrows <- which( CCF.data$clones %in% daughters )
    
    # get CCF values #
    parent.CCF <- CCF.data[ parentrow, "CCF" ]
    daughter.total.CCF <- sum( CCF.data[ daughterrows, "CCF" ] )
    
    # if decrease parent == TRUE & daughter CCF > parent CCF then increase parent CCF to match daughters, then if #
    # parent CCF > 1 then adjustt all CCFs on the tree to allow parent CCF < 1 #
    # if decrease daughter == TRUE & daughter CCF > parent CCF then decrease daughter CCFs to match paarent, then if #
    # parent adjust allows soome wiggle room for niose in CCF data - default = 0  howeever #
    if( parent.adjust * parent.CCF < daughter.total.CCF ){
      if( daughter.total.CCF / parent.CCF > warning.limit ){
        if( decrease.daughters  )  type <- "Decreasing daughter CCFs proportionally" else type <- "Increasing parent CCF"
        cat( paste0("        ", "clone ", parent, "'s daughters have total CCF which is ", signif( (daughter.total.CCF * clonal_CCF) / parent.CCF, 3 ), "% its own CCF. ", type, "  so total CCF of daughters = parent\n") )
      }
      if( increase.parents  ) CCF.data[ parentrow, "CCF" ] <- daughter.total.CCF * parent.adjust
      if( decrease.daughters  ) CCF.data[ daughterrows, "CCF" ] <- sapply(daughterrows, function(rowi) (CCF.data[ rowi, "CCF" ] / daughter.total.CCF) * parent.CCF )
    }
  }
  
  # normalise to root #
  # if the clonal cluster CCF needed to be increased then adjust this back down too 1 and adjust all other clones by similar margin #
  clonalCCF.change <- CCF.data$CCF[1] / clonal_CCF
  if( clonalCCF.change > 1 ) cat( paste0("        Clonal CCF needed increasing by ", clonalCCF.change * 100, "% to accommodate daughters. Therefore decreasing all CCFs proportionally so clonal CCF == 1 again\n") )
  
  # ensure clonal cluster = 1 #
  CCF.data$CCF <- ( (CCF.data$CCF * clonal_CCF) / CCF.data[ CCF.data$clones == root, "CCF"] ) 
  
  return( CCF.data )
}






plot_grouped_cloneMaps <- function(tumour_group_data, group_name, rasters.list, byRegion = FALSE, group_order = NA, track = TRUE, no.cols = 20 ){
  
  groups <- tumour_clin[ order(get(group_name)), .N, by = get(group_name)]
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
                             out <- c(out, tumour_clin[ get(group_name) == group, tumour ])
                             igroup <- which(groups$group == group)
                             out <- c(out, rep(NA, groups[igroup, rounded_n] - length(out) + 1))
                             return(out)
                           })) )
  
  
  ## now plot all the tumour maps ##
  
  
  #pdf(file = paste0(outputs.folder,"/",date,"_421_tumour_maps_by_Tumour_by_Stage.pdf"), width = 20, height = 30)
  
  layout(layout)
  
  for( i in 1:nrow(positions) ){
    
    is.sample <- grepl( "LTX", positions[i, to_plot] )
    
    par( mai = c(0, 0, 0, 0), xpd = NA)
    
    if( is.sample == TRUE ){
      
      tumour <- positions[i, to_plot]
      
      if(track == TRUE) cat( paste0( which(positions[ grepl("LTX",to_plot), to_plot] == tumour ), " " ) )
      
      cloneMap::cloneMap( clone_map = rasters.list[[ which( names(trees) == tumour ) ]] )
      
    } else {
      
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      
      if( !is.na( positions[i, to_plot] ))  text(x = 0.5, y = 0.5, positions[i, to_plot], cex =8, col = "black")
      
    }
    
  }
  
  #invisible( dev.off() )
  
}



##########################################
### supprting functions (not exported) ###
##########################################

# function to extract (x = row, y = col) coordinates from TRUE in boolean matrix #

matrix.index.to.coordinates <- function( matrix.index, nrow, ncol ){
  
  return( data.frame( x = (matrix.index %% nrow), y = floor( matrix.index / ncol ) + 1, stringsAsFactors = F) )
  
}

# function to convert (x = row, y = col) coordinates into a boolean matrix #
# can input either single coordinate or data.frame of coordinates also     #
# need to input matrix outline (number of cols and rows)                   #

coordinates.to.matrix.index <- function(coordinates, nrow, ncol){
  
  if(class(coordinates)=="numeric"){
    return( ( ( coordinates[1] - 1 ) * ncol ) + coordinates[2] )
  }
  
  if(class(coordinates)=="data.frame"){
    return( ( ( coordinates$x - 1 ) * ncol ) + coordinates$y )
  }
  
}

# function to test whether a specific clone in a raster matrix of clone positions is entirely continuous #
# ie does the clone seperate into different 'islands'. This is porobably biologicallly implausible so    #
# should be avioded if possible                                                                          #

continuous.test <- function( clone_position ){

  # how many core do you haave access to for parrelellisation? #
  num_cores <- parallel::detectCores() / 2
  # only worth parrellelising if > 10 cores, otheriwise actually slows code! - need to test further #
  if( num_cores < 10 ) num_cores <- 1
  # get cordinataes of all position for this clone #
  coords <- matrix.index.to.coordinates( which( clone_position ), nrow = nrow( clone_position ), ncol = ncol( clone_position ) )  # matrix.index.to.coordinates specified above
  coords_ids <- paste(coords[, "x"] , coords[, "y"], sep = "_" )
  
  # make as list #
  coords.list <- lapply( 1:nrow( coords ), function(i) coords[ i ,])
  
  
  # get surronuding cords coordinates for each cooord #
  coords.surrounding.id.list <- parallel::mclapply( coords.list, function(coord) c( paste( coord$x + 1, coord$y, sep = "_" ),
                                                                          paste( coord$x - 1, coord$y, sep = "_" ),
                                                                          paste( coord$x, coord$y + 1, sep = "_" ),
                                                                          paste( coord$x, coord$y - 1, sep = "_" ),
                                                                          # diagonals # 
                                                                          paste( coord$x + 1, coord$y - 1, sep = "_" ),
                                                                          paste( coord$x + 1, coord$y + 1, sep = "_" ),
                                                                          paste( coord$x - 1, coord$y - 1, sep = "_" ),
                                                                          paste( coord$x - 1, coord$y + 1, sep = "_" ) ), mc.cores = num_cores )
  
  
  # choose a postion in clone at random #
  # get all its neighbours are in the clone, then all thier neighbours recursively until # of coords is no longer increasing #
  # then check if we've captured all clone coords, if not then there must must be non-continuous islands #
  
  # this also holds true for just the edges and this is faster so limit to these #
  
  # identify edge positions #
  # these are those coords surrounded by at least 1 coord which is not part of the clone #
  is_edge <- sapply( 1:length( coords.list ), function(i) any( !coords.surrounding.id.list[[i]] %in% coords_ids  ) )
  
  # now make a lists of all the edge coords and thier surrrounding coords #
  coords_ids_edge <- coords_ids[ which( is_edge ) ]
  coords_edge_surrounding_id_list <- coords.surrounding.id.list[ is_edge ]
  
  # chose an edge coord at random #
  edges_continuous <- coords_ids_edge[[ sample( 1:length( coords_ids_edge ), 1 ) ]]
  
  repeat{
    # add neightbouring edges #
    
    neighbours.index <- which( sapply(coords_edge_surrounding_id_list, function( surronding_coords ) any( surronding_coords %in% edges_continuous )) )
    
    edges_continuous_new <- unique( c( edges_continuous, coords_ids_edge[ neighbours.index ] ) )
    
    if( length(edges_continuous_new) == length(edges_continuous) ) break
    
    edges_continuous <- edges_continuous_new
    
  }
  
  # are all edge coords contained in these continuous edges from this nucleation? #
  is_continuous <- all( coords_ids_edge %in% edges_continuous ) 
  
  # if so must be continuous #
  
  return( is_continuous)
  
}

# function to recentre distance matrix  #

recenter_distance_matrix <- function( clone_position ){
  
  ## find approximately the new centre of the clone ##
  
  # determine in which columns and rows the clone is now present #
  cols_present <- which( apply( clone_position, 2, any ) )
  rows_present <- which( apply( clone_position, 1, any ) )
  
  # from this get approx centre cordinates #
  x <- min( rows_present ) + ( ( max( rows_present ) - min( rows_present ) ) / 2 )
  y <- min( cols_present ) + ( ( max(cols_present) - min(cols_present) ) / 2 )
  new_nucleus <- c( x, y )
  
  ## use this as nucleus for a new distrebution ##
  new.dist <- make.distance.matrix( clone_position, new_nucleus ) # make.distance.matrix function specified above
  
  ## make sure new nucleus is within clone
  new_nucleus_in_clone <- any( new.dist == 0 & clone_position )
  
  ## if not reassign the nucleus to the nearest part of the clone ##
  if( !new_nucleus_in_clone ){
    
    # matrix.index.to.coordinates ooupuuts data.frame just take the first option -       #
    # could be muplie nuclei if if several positions have minimum disttnce to new nulcei #
    nucleus.index <- which( new.dist == min( new.dist[ clone_position ]) )
    new_nucleus <- as.numeric( matrix.index.to.coordinates( matrix.index = nucleus.index, nrow =  nrow(clone_position),  ncol = ncol(clone_position))[ 1 ,]) # make.distance.matrix function specified above
    
    new.dist <- make.distance.matrix( clone_position, new_nucleus )
    
  }
  
  return(new.dist)
  
}

#' function to round to a specific base 

mround <- function( x, base ) base * round( x / base )

#' function to roun down to a specific base 

mfloor <- function( x, base ) base * floor( x / base )

#' function to round up to a specific base 

mceiling <- function( x, base ) base * ceiling( x / base )




#===========#
###  END  ###
#===========#
