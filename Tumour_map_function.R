#####

# Function to produce plots showing clonal composotion of a tissue #
# requires phylogenetic tree data for each clone and the Cancer    #
# Cell Fraction (CCF) of each clone                                #

Clonal_map <- function(tree.mat, CCF.data = NA, raster.data = NA, resolution.index = 100, 
                       brewer.palette = "Paired", smoothing.par = 10, clone.cols = NA, 
                       output.rastered.data = FALSE, plot.data = TRUE, border.thickness = 1.5,
                       border.colour = "grey20", repeat.limit = 10){
  
  
  ## load libraries required in this function ##
  suppressPackageStartupMessages( library(qlcMatrix) )
  suppressPackageStartupMessages( library(sf) )
  suppressPackageStartupMessages( library(smoothr) )
  suppressPackageStartupMessages( library(raster) )
  suppressPackageStartupMessages( library(RColorBrewer) )
  
  
  # ensure tree is correct class #
  tree.mat <- as.matrix( tree.mat )
  
  # get colours for plotting if these are not provided #
  if( all( is.na(clone.cols) ) ){
    
    # order the tree so the trunk and earl clones are always the same colours accross tumours #
    if( nrow(tree.mat > 1) ) tree.mat <- logically.order.tree(tree.mat)
    
    clones <- unique( as.numeric(tree.mat) )
    getPalette <- colorRampPalette( brewer.pal(9, brewer.palette) ) # brewer.palette specified in argumentss, default = "Paired"
    clone.cols <- getPalette( length( clones ) )
    names(clone.cols) <-clones
    
  }
  
  ## if the raster data already supplied, skip the main steps and just plot the rasterised data ##
  
  if( all( is.na(raster.data) ) ){
    
    #####################################################
    ### Ensure all input data is clean and consistent ###
    #####################################################
    
    # ensure CCF.table is correct class #
    CCF.data <- as.data.frame( CCF.data )
    
    # remove any clones with estimated CCF of 0 in all regions #
    CCF.data <- CCF.data[CCF.data$CCF>0,]
    
    # limit CCF.table to clones on the tree #
    # if this is the case print those which have been removed #
    tree.clones <- unique( as.numeric(tree.mat) )
    
    if( !all( CCF.data$clones %in% tree.clones ) ) print( paste0("clones ", paste0(setdiff(CCF.data$clones, tree.clones), collapse = " "), " found in CCF table but not phylogenetic tree. These have been removed." ) ) 
    CCF.data <- CCF.data[ CCF.data$clones %in% tree.clones, ]
    
    # limit tree to clones in the CCF.table #
    # if this is the case print those which have been removed #
    if( !all( CCF.data$clones %in% tree.clones ) ) print( paste0("clones ", paste0(setdiff(tree.clones, CCF.data$clones), collapse = " "), " found in phylogenetic tree but not in CCF table. These have been removed." ) )
    tree.mat <- remove.clones.on.tree(tree.mat,clones.to.keep = CCF.data$clones) ## Function specified below and in FrankellA_functions.R script
    
    # if tree has only one relationship, class becomes numeric/character after last step - correct this #
    if(class(tree.mat)=="character" | class(tree.mat)=="numeric"){
      
      tree.mat <- matrix(tree.mat,ncol = 2,byrow = TRUE)
      
    } else {
      
      # if tree has > one relationship, order the tree from Trunk -> branches -> leaves #
      if(nrow(tree.mat)>1){
        
        tree.mat <- logically.order.tree(tree.mat) ## Function specified below and in FrankellA_functions.R script
        
      }
    }
    
    ### specify which clone is the root ###
    
    #  as tree as has beeen ordered this been be the parent (ie column 1) in the first relationship #
    root <- tree.mat[1,1]
    
    ### In cases where CCF all daughters > CCF parent increase the parent CCF  to accomodate its daughters ###
    # If this occurs a warning will be outputted with how much smaller the parent CCF is than its childern #
    # given noise in CCF caluculations mean we can accept some underestimate of parent CCF but if children # 
    # > 130% of parent this should be checked and tree/clones/CCFs may be incorrect                        #
    CCF.data <- make.CCFs.tree.consistant(tree.mat = tree.mat, CCF.data = CCF.data,parent.adjust = 1)   ## Function specified below and in FrankellA_functions.R script
    
    
    ######===============================================================######
    ######                                                               ######
    ######   Create rasterised data to indicate postions of each clone   ######
    ######                                                               ######
    ######===============================================================######
    
    # make matrix that specifies each position using the resolution index (argument) specifying number fo columns and rows #
    # the higher the specified resultion index the more precise the plots, but the longer this code will take to run #
    # populate the plot with 0s (indicates no prescence of any clone) #
    
    clones_rasterised <- do.call(cbind, lapply(1:resolution.index,function(i) rep(0,resolution.index)))
    
    ##############################
    ### Add the clonal cluster ###
    ##############################
    
    # clones are present in "patches" within the plot which simulates (v approximately) how clonal growth occurs #
    # Each clone will nucleate at a partially randomised postion within its parent clone and grow until it reaches its specified CCF #
    # sisters within the same parent will compete for space and grow around one another until both reach the appropriate CCF #
    # clones grow in this mathematically simple octagon shape unless they hit barrier (parent edge or sister) #
    
    # first get nucleus coordinate for clonal cluster (center of the plot) #
    nucleus <- c( resolution.index / 2, resolution.index / 2 )
    
    # now calculate distance from nucleus for all pixels #
    # this distance matrix is used to work out how the clone should grow #
    dist.mat <- make.distance.matrix( martix.outline = clones_rasterised, nucleus = nucleus ) # function specified below
    
    # make clonal area the diameter of the whole plot #
    clonal.area <- sum( dist.mat < resolution.index / 2 )
    
    ## Calcualate equivilent 'area' each subclone should consume depending on its CCF ##
    # CCF should be supplied as a percentage - converted to a fraction #
    CCF.data$area <- round( (CCF.data$CCF / 100) * clonal.area )
    
    # make sure clonal is 100% (sometimes ~99.5 how we calculate it) #
    CCF.data[CCF.data$clones == root, "area"] <- clonal.area
    
    # now determine cut off in distance matrix which results in desired amount of area for th clonal cluster #
    possible_cutoffs <- seq(0.1,mround(max(dist.mat),0.1),0.1)
    distance_cutoff <- min( possible_cutoffs[ sapply(possible_cutoffs, function(cut) sum(dist.mat<cut)) >= clonal.area ] )
    
    #save the blank version of the raster matrix #
    clones_rasterised_blank <- clones_rasterised
    
    # assign clonal clone to its positions
    clones_rasterised[dist.mat < distance_cutoff] <- root
    
    # assign all other positions in plot to Inf ( indicates no clones ) #
    dist.mat[!dist.mat < distance_cutoff] <- Inf
    
    #  plot.data == true if you want to plot the rasterised data as it is generate, if plot = FLASE   #
    #  you can set output.Rasterised.data == TRUE then instead of plotting the function will output   #
    #  the rasterised clone data to save a plot whenever you like. This is particularly useful as     #
    #  clones are seeded at random so by saving the rasterised data you can ensure the clones don't   #
    #  change positins ech ttime you plot them. You can input the rasterised data & instead of a CCF  #
    #  table with the raster.data arguemnt                                                            #
    
    if( plot.data == TRUE ){
      
      # ensure raster is class numeric not char #
      clones_rasterised_plot <- apply(clones_rasterised,1,as.numeric)
      
      # set up plot extent #
      plot(rasterToPolygons(raster(clones_rasterised_plot)), col = NA, border = NA) 
      
      ## plot the clonal clone ##
      plot <- st_as_sf(rasterToPolygons(raster(clones_rasterised_plot), function(x){x == root}, dissolve = TRUE))
      
      # add smoothnig to make it more visually appealling #
      # note: smoothing will make the plotted area slightly underestmiate the true area but by only very little #
      # also has the advantage that you can see parent clones below at smooth edges and it becomes easier too   #
      # see tree relatinoships                                                                                  #
      plot.smooth <- smooth(plot, method = "ksmooth", smoothness= smoothing.par)
      
      # specify border thickness & colour in arguments - default is 1.5 & grey #
      plot(plot.smooth, col = clone.cols[names(clone.cols)==root], border = border.colour, lwd = border.thickness, add = TRUE) 
      
    }
    
    ##################################
    ### Add the subclonal clusters ###
    ##################################
    
    # get all parental clones and deal with the daughters of these in turn, starting with the daughters of the clonal cluster #
    # this will be in order of trunk -> branch -> leaf as tree has been ordered as such #
    # col 1 of tree always = parents #
    parents <- unique(tree.mat[,1])
    
    # when loop around parents, we need to record the distance matrices and cut offs used for each parental clone #
    # start by recording this for the clonal cluster #
    parent.dists <- list(dist.mat)
    names(parent.dists) <- root
    parent_distance_cutoffs <- c(distance_cutoff)
    names(parent_distance_cutoffs) <- root
    
    # loop around each parent #
    
    for(parent in parents){
      
      # for testing # print(parent)
      
      # if only got the clonal cluster don't need to continue - already recorded and plotted the clonal cluster #
      if( all(tree.mat == root) ) break
      
      # extract daughters of this clone and how many there are #
      # daughters are always col 2 of tree #
      daughters <- tree.mat[ tree.mat[,1] == parent, 2]
      levels <- length(daughters)
      
      # from the saved lists,  extract the parent's distance matrix and cut off for the boundry #
      parent.dist.mat <- parent.dists[[which(names(parent.dists) == parent)]]
      parent_distance_cutoff <- as.numeric(parent_distance_cutoffs[which(names(parent_distance_cutoffs) == parent)])
      
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
        
        repeati <- repeati + 1
        clones_rasterised <- clones_rasterised_parent
        
        ########################################################
        ### determing position of nuclei for daughter clones ###
        ########################################################
        
        # choose nuclei positions near the centre of the parent but also not too near the edge #
        # if muliple clones they must nucleate a % of the parent diameter away from each other #
        # otherwise random #
        
        # if just one subclone (levels = 1) nucleate subclone very near the centre of the parent #
        # if > 1 subclone (levels > 1) nucleate subclones as far asway as possible from one another #
        # while stll maintainig a certain distance from the edge and centre of the parent clone #
        
        
        if(levels > 1){
          
          
          # normally allow clones to nucleate nither too fro away or too close to centre of parent
          nucleus.options <- parent.dist.mat > (parent_distance_cutoff * 0.40) & parent.dist.mat < (parent_distance_cutoff * 0.80)
          
          # if struggling to find continuous solution may need to allow clones to nucleate nearer the edges #
          if( repeati > 5 ) nucleus.options <- parent.dist.mat > (parent_distance_cutoff * 0.30) & parent.dist.mat < parent_distance_cutoff 
          
          # convert raster matrix TRUE positions to cordinates for possible nucleation #
          nucleus.options <- matrix.index.to.coordinates(matrix.index = which(nucleus.options),nrow=nrow(clones_rasterised),ncol=ncol(clones_rasterised))
          
          # randomly select 20 sets of n (n = levels) nuclei for clones - might be better which more but distance caluclatioons take a while with more #
          nucleus.options.sel <- lapply(1:20, function(x) nucleus.options[sample(1:nrow(nucleus.options),levels,replace = F),])
          
          # for each option work out the average distance between the sets of nuclei
          nucleus.options.min.dists <- sapply( nucleus.options.sel, function( nucleus.option ){
            
            mindist <- sapply( 1:nrow(nucleus.option), function(i){
              
              dists <- make.distance.matrix( clones_rasterised.blank,nucleus = as.numeric( nucleus.option[i,] ) )
              dists <- sapply( which( !1:nrow( nucleus.option ) == i ), function(i2) dists[ nucleus.option[ i2, "x" ], nucleus.option[ i2, "y" ] ])
              
              return(min(dists))
            })
            
            return(min(mindist))
          })
          
          # randomly choose the set fo nuclei which is in the upper fifth of distances apart #
          nucleus.options.sel <- nucleus.options.sel[ nucleus.options.min.dists > quantile( nucleus.options.min.dists, 0.80 ) ]
          nuclei <- nucleus.options.sel[[ sample( 1:length(nucleus.options.sel ), 1 ) ]]
          nuclei <- lapply( 1:nrow(nuclei), function(x) as.numeric( nuclei[x, ] ) )
          names(nuclei) <- daughters # nuclei = the set of nuclei chosen for the daughters of this parent
          
        } else {
          
          # if only one adughter clone just place very near the centre (0.3 * distance to edge) of the parent #
          nucleus.options <- which( parent.dist.mat < (parent_distance_cutoff * 0.3) )
          nucleus.options <- nucleus.options[ sample( 1:length( nucleus.options ), 1 ) ]
          nuclei <- matrix.index.to.coordinates( nucleus.options, nrow = nrow( clones_rasterised ), ncol = ncol( clones_rasterised ) )
          nuclei <- lapply(1:nrow(nuclei), function(x) as.numeric(nuclei[x,]))
          names(nuclei) <- daughters
          
        }
        
        ##############################################################
        ### now allow clones to grow from nuclei around each other ###
        ##############################################################
        
        # extract new distance matricies for each of these new nuclei #
        nuclei.dists <- lapply( 1:length( nuclei ), function(i) make.distance.matrix( clones_rasterised, nucleus = nuclei[[i]], type = "octoagon" ) )
        
        # set regions which are off limits (outside of parental clone) to Inf #
        nuclei.dists <- lapply(nuclei.dists, function( nuclei.dist ){
          nuclei.dist[!clones_rasterised == parent] <- Inf
          return( nuclei.dist )
        })
        names( nuclei.dists ) <- daughters
        
        # just check that nuclei are defninitely within the parent (found this error a few times) #
        # if nucleus outside parent then will have been set to Inf therefore the matrix will lack a 0 #
        if( any( sapply( sapply( 1:length(daughters), function(i) which( nuclei.dists[[i]] == 0 )), length) == 0) ){
          stop("chosen nucleus outside of parent")
        }
        
        # grow clones regularly to allow clone to grow evenly relative to one another #
        
        # define targe area for each clone #
        clone.areas <- sapply(names(nuclei), function(clone) CCF.data[CCF.data$clones==clone,"area"])
        
        # define stage to grow the clones at the same rate - 100 nits fo area at a time to the max of all daughter clones #
        growth.stages <- seq(100, mceiling( max(clone.areas), 100), 100)
        
        for(area in growth.stages){
          
          # for this round of growth determine area of expansion for each clone, limiting to max area #
          areas <- rep( area, length(daughters) )
          names(areas) <- daughters
          area.too.large <- sapply(daughters, function(daughter) areas[ names(areas) == daughter ] > clone.areas[ names(clone.areas) == daughter ] ) 
          areas[ area.too.large ] <- clone.areas[ area.too.large ]
          
          # now determine cut off which results in desired number of pixels #
          avialible.space <- lapply( daughters, function(clone) as.numeric( nuclei.dists[[ which( names( nuclei.dists ) == clone ) ]] ) )
          avialible.space <- lapply(1:length( daughters ), function(i) avialible.space[[i]][ avialible.space[[i]] < Inf ])
          avialible.space <- sapply(1:length( daughters ), function(i) max( avialible.space[[i]] ))
          distance_cutoffs <- sapply( daughters, function(clone){
            rounded.area <- mround( avialible.space[[ which (names(nuclei) == clone) ]], 0.1 ) ## mround function specified below
            cut.options <- seq(0.1,rounded.area,0.1)
            cut.options <- cut.options[ sapply(cut.options, function(cut) sum( nuclei.dists[[ which( names( nuclei ) == clone ) ]] < cut)) <= areas[ names(areas) == clone ]]
            return( max(cut.options) )
          })
          
          # now record the areas on the rasaterised plot with the clone names #
          for(clone in daughters){
            clones_rasterised[ nuclei.dists[[ which( names( nuclei.dists ) == clone ) ]] < as.numeric( distance_cutoffs[ names( distance_cutoffs ) == clone ] ) ] <- clone
          }
          
          # set regions which are off limits to this clone (part of sister clone) to Inf #
          nuclei.dists <- lapply( 1:length( nuclei.dists ), function(i){
            sisters <- names( nuclei )[ which( !1:length(nuclei) == i ) ]
            nuclei.dist <- nuclei.dists[[i]]
            nuclei.dist[ clones_rasterised_plot %in% sisters ] <- Inf
            return( nuclei.dist )
          })
          names(nuclei.dists) <- daughters
          
        }
        
        for(clone in names(nuclei)){
          parent.dists <- c( parent.dists, list( nuclei.dists[[ which( names( nuclei.dists ) == clone ) ]] ) )
          names(parent.dists) <- c( names(parent.dists)[ 1:( length(parent.dists) -  1 ) ], clone )
          parent_distance_cutoffs <- c( parent_distance_cutoffs, distance_cutoffs[[ which( names( distance_cutoffs ) == clone ) ]] )
          names(parent_distance_cutoffs) <- c( names( parent_distance_cutoffs )[ 1:( length(parent_distance_cutoffs) - 1 ) ], clone)
        }
        
        # test if all clones are continuous - if not repeat #
        if(all(sapply(names(nuclei), function(clone) continuous.test(clone,clones_rasterised_plot)))) break
        
        # only allow certain number of repeats and then give up (default = 10) #
        if( repeati == repeat.limit ){ print( paste0( "reached repeat limit to achieve continuous daughter clones in parent clone", parent ) ) ; break }
        
        
      }
      
      ################################################################
      ### Now plot the new clones which we have assigned positions ###
      ################################################################
      
      if( plot.data == TRUE ){
        
        for(clone in names(nuclei)){
          clones_rasterised_plot <- apply(clones_rasterised_plot,1,as.numeric)
          plot <- st_as_sf(rasterToPolygons(raster(clones_rasterised_plot), function(x){x == clone}, dissolve = TRUE))
          plot.smooth <- smooth(plot, method = "ksmooth", smoothness= smoothing.par) # smoothing par speicified in arguemnts
          plot(plot.smooth, col = clone.cols[names(clone.cols)==clone], border = "grey20", lwd = border.thickness, add = TRUE) # border thickness specified in arguemnts and col can be specified in arguments
        }
        
      }
      
    }
    
    # if specified in arguments output the raster matrix of clone positions #
    # this enables yu to repeatedly make the exact same plot or otherwise the clone positions will change each time #
    
    if( output.rastered.data == TRUE ) return( clones_rasterised_plot )
    
    
  } else {  ## if !is.na(raster.data)
    
    ######===================================================================================######
    ######                                                                                   ######
    ######   if you have supplied a raster matrix of clone positions then use this to plot   ######
    ######                                                                                   ######
    ######===================================================================================######
    
    clones_rasterised_plot <- raster.data
    
    # remove any clones n the tree which are not in the rasterised plot data #
    
    tree.mat <- remove.clones.on.tree(tree.mat,clones.to.keep = unique( as.character(clones_rasterised_plot) ) )
    
    #restore tree to martix if only 1/2 clones left and if > 1 relationship then oordr the tree trunk -> branch -> leaf
    
    if( class( tree.mat ) == "character" | class( tree.mat ) == "numeric" ){
      tree.mat <- matrix( tree.mat, ncol = 2, byrow = TRUE )
    } else {
      if( nrow( tree.mat ) > 1 ){
        tree.mat <- logically.order.tree( tree.mat )
      }
    }
    
    # simulate what occurs when plotting & raster generated concurrently #
    # for each clone in the tree (trunk -> leaves) plots the area this occupies indluing all its daughters #
    
    clones_rasterised_plot <- clones_rasterised_plot
    clones_rasterised_plot[] <- 0 
    clones_rasterised_plot <- apply( clones_rasterised_plot, 1, as.numeric )
    blank.plot <- rasterToPolygons( raster( clones_rasterised_plot ) )
    plot( blank.plot, col = NA, border = NA ) # set up plot extent
    
    for( clone in unique( as.numeric(tree.mat) ) ) {
      
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
      
      clones_rasterised_plot <- clones_rasterised_plot
      clones_rasterised_plot[ clones_rasterised_plot %in% clone.daughters ] <- clone
      clones_rasterised_plot <- apply(clones_rasterised_plot,1,as.numeric)
      
      # plot the clone #
      
      plot <- st_as_sf(rasterToPolygons(raster(clones_rasterised_plot), function(x){x == clone}, dissolve = TRUE))
      plot.smooth <- smooth(plot, method = "ksmooth", smoothness= smoothing.par) # smoothing par speicified in arguemnts
      plot(plot.smooth, col = clone.cols[names(clone.cols)==clone], border = "grey20", lwd = border.thickness, add = TRUE) # border thickness specified in arguemnts and col can be specified in arguments
      
    }
    
  }
  
  ###########
  ### END ###
  ###########
  
}

## function to remove clones from a phenogenetic tree matrix while    ##
## maintaining parent -> daughter reltionships,  even if intermediate ##
## (branch) clone is being removed                                    ##

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
  if(length(clones.to.remove)==0){
    warning( "no clones specified to remove that are on the tree")
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

## function to order tree root -> branches -> leaves ##

logically.order.tree <- function(tree){
  
  # if only 1 clone on tree then just output it as it is #
  if(all(tree==unique(tree)[1])) return(tree)
  
  
  ### assign levels to each parent clone depending on how near the trunk ###
  
  # work out the root (the only clone tht's never a daughter ) #
  root <- tree[,1] [! tree[,1] %in% tree[,2] ]
  
  # make empty list for levels of ecah clone #
  levels <- rep(NA,nrow(tree))
  
  # aasign the root to level 1 #
  levels[tree[,1] %in% root ] <- 1
  
  # list aall the daughter's of root to work out the level #
  daughters <- tree[ tree[,1] %in% root,2]
  
  # go down the tree in levels and asssign clenns correct levels #
  l <- 1
  repeat{
    l <- l + 1
    parents <- daughters
    levels[tree[,1] %in% parents ] <- l
    daughters <- tree[ tree[,1] %in% parents,2]
    if(all(!is.na(levels))) break
  }
  tree <- tree[order(levels),]
  return(tree)
}

## function to generate a distnce martix, specifiying how far away each piont in the raster matric   ##
## is from a central nucleus. So far managed two ways to calculaate this, one that specifiesdistance ##
## purely in two directions ( up and down ) hence casues a square-ish growth pattern and one that    ##
## accounts for a smaller than 2 distaance for diagonal movement which creates an octagonal pattern  ##
## the latter is the default for now, ideally probably like to create a circular pattern but not     ##
## sure how too do this                                                                              ##  

make.distance.matrix <- function( martix.outline = clones_rasterised, nucleus = nucleus, type = "octoagon" ){
  
  
  ### just count rows & cols from nucleus (makes 'square' shape ) ###
  
  if( type == "square" ){
    
    ## now calculate distance from nucleus for all pxls ##
    
    # dist with horizontal as primary #
    
    dist.mat.h <- do.call( cbind, lapply( 1:ncol(clones_rasterised ), function( col ){
      
      # how many columns from nucleus? #
      base.dist <- abs( nucleus[2] - col ) 
       
      return( rep( base.dist, nrow( clones_rasterised ) ) )
      
    } ))
    
    # dist with vertical as primary #
    
    dist.mat.v <- do.call( rbind, lapply( 1:nrow( clones_rasterised ), function( row ){
      
      # how many rows nucleus? #
      base.dist <- abs( nucleus[1] - row )
      
      return( rep( base.dist, ncol( clones_rasterised ) ) )
      
    } ))
    
    # now combine for shortest poss distance #
    dist.mat <- do.call( cbind, lapply( 1:ncol( clones_rasterised ), function( col ){
      
      out <- cbind( dist.mat.h[, col ], dist.mat.v[, col ])
      return( as.numeric( rowSums( out ) ) )
      
    }))
    
  }
  
  ### if accounting for diagonal distance type == 'octogon' ###
  
  if( type == "octoagon" ){
    
    ## now calculate distance from nucleus for all pxls ##
    
    # dist with horizontal as primary #
    
    dist.mat.h <- do.call( cbind, lapply( 1:ncol( clones_rasterised ), function( col ){
      
      # how many columns from nucleus? #
      base.dist <- abs( nucleus[2] - col )
      
      # addition allow us to account for diagnal movement with change in rows #
      horizontal.additions <- abs( 1:nrow( clones_rasterised ) - nucleus[1] ) * 0.41
      
      return( base.dist + horizontal.additions )
      
    } ))
    
    
    # dist with vertical as primary #
    
    dist.mat.v <- do.call( rbind, lapply( 1:nrow( clones_rasterised ), function( row ){
      
      # how many rows from nucleus? #
      base.dist <- abs( nucleus[1] - row )
      
      # addition allow us to account for diagnal movement with change in cols #
      vertical.additions <- abs( 1:ncol( clones_rasterised ) - nucleus[2] ) * 0.41
      
      return( base.dist + vertical.additions )
      
    } ))
    
    # now combine for shortest poss distance #
    
    dist.mat <- do.call( cbind, lapply( 1:ncol( clones_rasterised ), function( col ){
      
      out <- cbind( dist.mat.h[, col ], dist.mat.v[, col ])
      return( as.numeric( rowMax( out ) ) )
      
    }))
    
  } 
  
  # return output #
  return( dist.mat )
  
  #########
  ## END ##
  #########
  
}

# function to extract (x = row, y = col) coordinates from TRUE in boolean matrix #

matrix.index.to.coordinates <- function(matrix.index,nrow,ncol){
  
  return( data.frame( x = (matrix.index %% nrow), y = floor( matrix.index / ncol ) + 1, stringsAsFactors = F) )

}

# function to convert (x = row, y = col) coordinates intoo a boolean matrix #
# can input either single coordinate or data.frame of coordinates           #
# also need to input matrix outline (number of cols and rows)               #

coordinates.to.matrix.index <- function(coordinates,nrow,ncol){
  
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

#at the moment this doesn't work very well - will only pick up really bad cases but this is fine for now
# disabled this for now - set P value theshold to 1.1
continuous.test <- function( clone, df ){
  
  # get cordinataes of all position for this clone #
  coords <- matrix.index.to.coordinates( which( df == clone), nrow = nrow( df ), ncol = ncol( df ) )
  
  # get the nucleus for this clone #
  nucleus <- as.numeric( coords[ sample( 1:nrow( coords ), 1 ) ,] )
  
  # get a distance matrix for this clone #
  dist.mat <- make.distance.matrix( df, nucleus = nucleus, type = "octagon")
  
  # specify Inf if not the clone #
  dist.mat[  !df == clone ] <- Inf
  
  
  p<-dip.test(dist.mat[df==clone])$p
  if(p > 1.1){
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# function to round to a specific base #
mround <- function( x, base ) base * round( x / base )

