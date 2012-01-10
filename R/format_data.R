#' format data in ape format into ouch format
#' @param tree a phylogenetic tree of class "phylo", ape format, or a tree
#' in ouch format with the data also in ouch format
#' @param traits a numeric with trait values, or a matrix or data frame of 
#' traits, rownames matching species or handed in
#' @param species_names in the order of entries in traits, if not given in rownames.  
#' @param regimes the column in traits containing regime labels
#' @return the ouch-formatted tree, traits, and regimes  
#' @seealso \code{\link{convert}} to toggle between formats, including regime paintings
#' @details Should become an internal function to handle data conversion to ape-type
#' @import geiger
#' @export 
format_data <- function(tree, traits, species_names = NULL, regimes = NULL ){
# Function checks that tree and trait match and convert them into a format used by wrightscape
# Function also will code tree by finding the common ancestor of all species with matching entry specified in the regimes list and assigning that codename as the regime of all descendents of that ancestor.  May not handle conflicts if corresponding to overlapping clades.  Alternatively, the regimes can be specified directly in ouch format.   
	if( is(tree, "character" ) ){ 
		tree <- read.nexus(tree) 
	} else if (is(tree, "ouchtree")) {
		# uses my ouch2ape tree conversion script
    # should handle the data conversion of ouch2ape as well! 
		outtree <- convert(tree, regimes=regimes)
    if(is(traits, "data.frame"))
      dat <- traits[tree@nodelabels!="",1]
    else if(is(traits, "numeric"))
      dat <-  traits[tree@nodelabels!=""]
    else if(is(traits, "list"))
      dat <-  traits[[1]][tree@nodelabels!=""]
    else 
      stop(paste("traits class", class(traits), "not recognized"))
    names(dat) <- tree@nodelabels[!tree@nodelabels==""]
    out <- list(tree=outtree, data=dat)
	} else if(is(tree, "phylo")){
    # Figure out if species names is already attached to traits 
    if(is.null(species_names) ){
      if( !is.null(names(traits)) & is(traits, "numeric") ){
        species_names <- names(traits)
      } else if( !is.null(rownames(traits)) ){ 
        species_names <- rownames(traits) 
      } else { 
        stop("Species names not found")
      }
    }
    ## attach species names to traits
    if(is(traits, "numeric") ){
      names(traits) <- species_names
    } else if (is(traits, "data.frame") | is(traits, "matrix")){
    rownames(traits) <- species_names
    } else { stop("traits format unrecognized") }

    # drop missing taxa using geiger's function
    matched <- treedata(tree, traits, species_names)
    # treedata returns a matrix which loses class type for trait data.frame
    # so we restore the the dataframe structure
    if( is(matched$data, "matrix") ){ 
      matched$data <- as.data.frame(matched$data, stringsAsFactors=FALSE)
      if(is(traits, "data.frame")){
      for(i in 1:length(traits)){
          tmp <- class(traits[[i]])
          # each datarow that isn't numeric should be a character string (for regime, etc)
          if(tmp == "factor") tmp <- "character"
          class(matched$data[[i]]) <- tmp
        }
      }
    } 
    traits <- matched$data
    species_names <- rownames(matched$data) 
    # Convert to OUCH format 
    tree <- ape2ouch(matched$phy) 
    # Makes sure data is reformatted to ouch format matching the tree
    if( is(traits, "data.frame") ){
      message("traits ouch-formatted as data.frame")
      if(length(traits)>1){
          dataIn <- traits[match(tree@nodelabels, species_names),]
      } else {
        #hack to get around size-1 data-frames 
        traits[[2]] = NA
        dataIn <- traits[match(tree@nodelabels, species_names),]
        dataIn[[2]] <- NULL
      }
    } else if(is(traits, "numeric") ){ 
  message("traits ouch-formatted as numeric")
        dataIn <- traits[match(tree@nodelabels, species_names)]
    } else { stop(paste("data of class", class(traits), "not recognized")) }
    rownames(dataIn) <- tree@nodes
    if(is.numeric(regimes)){	
      R <- compute_regimes(tree, traits, species_names, regimes)
      nr <- R$noregimes
    } else {
      regimes= as.factor(rep(" ", length=tree@nnodes))
      names(regimes) <- tree@nodes 
      R <- list(regimes=regimes) 
      nr <- R$regimes
    }
	out <- list(tree=tree, data=dataIn, regimes=R$regimes, noregimes=nr)
	} else { 
    stop("Problem with tree format") 
  }
  out 
} 

#' internal method for computing the regimes used by format_data()
#' @param tree ouch-formated phylogenetic tree
#' @param traits ouch-formated data set
#' @param species_names for all the nodes 
#' @param regimes list in ouch format
#' @return a list of the regimes as edge-labels for ape's phylo format 
#' @keywords internal
compute_regimes <- function(tree, traits, species_names, regimes){

	## attach species names to regimes
	if(is.null(regimes) ){
		if( is(traits, "data.frame") | is(traits, "matrix") ){
			n <- dim(traits)[1]
		} else if (is(traits, "numeric")){
			n <- length(traits)
		}
		regimes <- as.factor(character(n))
		names(regimes) <- species_names
	} else if(length(regimes) == 1) {
		regimes <- traits[,regimes]
	}

	# order regimes into ouch format
	regimes <- regimes[match(tree@nodelabels, species_names)]

	# add an ancestral regime instead of NA, which confuses the mrcaOUCH function
	regimes <- as.character(regimes)
	regimes[is.na(regimes)] <- "anc"
	regimes <- as.factor(regimes)

	# add a noregimes case for ou1 fits
	noregimes <- as.factor(character(length(regimes)))

	# Paint regimes by clade using functions from maticce
	clades <- which(levels(regimes) != "anc")
	clade_ancestors <- integer(length(levels(regimes))-1)
	clade_names <- character(length(clade_ancestors) )
	k <- 1
	for(i in clades ){
		ancestor <-	as.integer(
			mrcaOUCH( tree@nodelabels[ regimes == levels(regimes)[i] ], tree)
		)
		clade_ancestors[k] <- ancestor
		clade_names[k] <- levels(regimes)[i]
		k <- k+1
	}
	## paintBranches needs the regime names in order
	A <- data.frame(clade_ancestors, clade_names)
	A <- A[order(A$clade_ancestors), ]
	regimes <- paintBranches(A$clade_ancestors, tree, regimeTitles=as.character(A$clade_names))

	# both regimes need node numbers as their names for OUCH
	names(noregimes) <- names(regimes)
	list(regimes=regimes, noregimes=noregimes)
}






