# Copyright (c) 2020, ETH Zurich


#' Allows the user to define the ecological consequences for species within each site,
#' defining thus species survival and abundance
#'
#' @details The arguments of the function allows to apply abiotic and biotic ecological rules to species in each 
#' site. Based on those rules, the function updates the abundance of each species in each site. If the abundance 
#' is null, the species is absent or extinct. Ecology can account for local environmental conditions, the abundance of
#' species, and/or their traits.
#'
#' @param abundance a named vector of abundances with one abundance value per species
#' @param traits a named matrix containing the species traits, one row per species
#' @param local_environment the environmental values for the given site
#' @param config the config of the simulation
#'
#' @return an abundance vector with the new abundance values for every species.
#' An abundance value of 0 indicates species death, any other values indicates survival.
#' @export
apply_ecology <- function(abundance, traits, local_environment, config) {
  stop("this function documents the user function interface only, do not use it.")
}



#' Orchestrates for applying the ecology function to all sites
#'
#' @details The ecology is applied on a per site basis over all species occurring in each site.
#' Therefore this function iterates over all sites and collects the abundance and traits of any species occurring there.
#' It then calls the user supplied apply_ecology function to this collection and apply ecology to each site.
#'
#' @param config the general config of the simulation
#' @param data the general data list
#' @param vars the general variables list
#'
#' @return returns the standard val(config, data, vars) list
#' @noRd
loop_ecology <- function(config, data, vars) {
  # skip ecology function if config$exp$enable_eco_mec is FALSE
  if(config$gen3sis$general$verbose>=3){
    cat(paste("entering ecology module @ time", vars$ti, "\n"))
  }

  all_cells <- rownames(data$landscape$environment)
  all_species_presence <- do.call( cbind, lapply(data$all_species, FUN = function(sp) {all_cells %in% names(sp$abundance)}))
  rownames(all_species_presence) <- all_cells

  # take ids that have at least one species...
  #occupied_cells <- rownames(geo_sp_ti[rowSums(data$geo_sp_ti)>0, ,drop=FALSE])
  occupied_cells <- rownames(all_species_presence)[rowSums(all_species_presence)>0]

  for (cell in occupied_cells) { # strat loop over ids with at least one species...
    local_environment = data$landscape[["environment"]][cell, , drop=FALSE]

    coo_sp <- which(all_species_presence[cell,])
    #create coocuring species traits for idi
    traits <- matrix(nrow = length(coo_sp), ncol = length(config$gen3sis$general$trait_names))
    abundance <- numeric(length(coo_sp))

    #colnames(tr_sp) <- colnames(data$eco[[1]])[1:(length(config$exp$eco$trait_names)+1)]
    colnames(traits) <- config$gen3sis$general$trait_names

    i <- 1
    for (spi in coo_sp){ #loop over co-ocurring species @ idi
      # tr_sp_ti_idi[i,] <- data$eco[[spi]][idi,-(length(config$exp$eco$trait_names)+2)]
      traits[i,] <- data$all_species[[spi]][["traits"]][cell, config$gen3sis$general$trait_names]
      abundance[i] <- data$all_species[[spi]][["abundance"]][cell]
      i <- i+1
    }
    max_n_sp_idi <- config$gen3sis$general$max_number_of_coexisting_species
    if (length(coo_sp) > max_n_sp_idi) {
      vars$flag <- "max_number_coexisting_species"
      paste0("Maximum number of species per cell (i.e. max_n_sp_idi) reached. Specifically ",
                  length(coo_sp),"(>", max_n_sp_idi,") species @ t",vars$ti, " idi",cell )
      return(list(config = config, data = data, vars = vars))
    }

    rownames(traits) <- coo_sp
    names(abundance) <- coo_sp

    #species <- traits[, c("abd", config$gen3sis$general$trait_names), drop = FALSE]
  

    NEW_abd <- config$gen3sis$ecology$apply_ecology(abundance, traits, local_environment, config)

    # colnames(NEW_abd) <- coo_sp_ti_idi
    # TODO check if colnames(geo_sp_ti[,coo_sp_ti_idi]) should be used see line 622+-
    names(NEW_abd) <- coo_sp
    #abd_threshold <- config$exp$abundance_threshold
    shalldie <- NEW_abd == 0

    for (spi in coo_sp){
      data$all_species[[spi]][["abundance"]][cell] <- NEW_abd[[toString(spi)]]
    }

    die_sure <- as.integer(names(NEW_abd)[NEW_abd == 0])

    if (length(die_sure)>0) { #check if there are any species that should go extinct
      chars <- as.character(die_sure)
    } #end of check if any die_sure
  } #end loop over ids with at least one species...
  species_list <- list()
  for (species in data$all_species) {
    cells <- names(species[["abundance"]])[species[["abundance"]] != 0]
    updated_species <- limit_species_to_cells(species, cells)
    species_list <- append(species_list, list(updated_species))
  }
  data$all_species <- species_list

  if(config$gen3sis$general$verbose>=3){
    cat(paste("exiting ecology module @ time", vars$ti, "\n"))
  }
  return(list(config = config, data = data, vars = vars))
}
### Biotic interactions

#' niche model function to calculate niche position for species and their links 
#' based on Williams & Martinez 2000
#' @param C connectance of the food web
#' @param S number of trophic species in the food web
#' @return list with a matrics of the links betweeen each species and 
#' a vector with the niche position of each species

niche_model <- function(C,S) {
  
  # Niche position
  n <- sort(runif(S,0,1))
  
  # Range
  beta <- 1/(2*C)-1
  range <- rbeta(S,1,beta)*n
  
  #Centroid
  centroid <- runif(S,range/2,n)
  
  # Make the first niche a producer
  range[which.min(n)] <- 0
  
  # Evaluate the matrix of links
  L <- matrix(0,nr=S,nc=S)
  low <- centroid - range/2
  high <- centroid + range/2
  
  for(i in 1:S)
    for(j in 1:S)
      if(low[i] < n[j] && high[i] > n[j]) L[j,i] = 1
  
  return(list(L,n))
}


#' function to obtain an interaction matrix with alpha values based on the niche model
#' of William & Martinez
#' @param C connectance of the food web
#' @param S number of trophic species in the food web
#' @param name_file string with file name (relevant if export = T)
#' @param export boolean, whether or not interaction matrix should be exported, default is T
#' @param type which type of interaction matrix should be generated
#' type=1 totally random assignment of the alpha values
#' type=2 random attribution of the lower diagonal values and opposite values 
#' assigned to the upper diagonal
#' type=3 positive random values attribution for the lower diagonal and opposite
#' values assigned to the upper diagonal
#' type=4 Predation alpha matrix
#' type=5 Mutualism alpha matrix
#' type=6 Competition alpha matrix
#' @return interaction matrix with alpha values for the type indicated
#' @export

get_alpha_mat <- function (C, S, name_file = "alpha.txt", export=T, type=2){
  
  counter <- 0
  while (counter==0) {
    a <- niche_model(C,S)[[1]]
    if(sum(a)<S) {counter=0} # check if L matrix has more 1s than S
    else{
      counter=1
    }  # end of if else
  } # end of while
  
  rownames(a) <- colnames(a) <- seq(1,S,1)
  diag(a) <- 0  # make diagonal zero
  
  e <- cbind(stack(as.data.frame(a)),b=seq(1,S,1)) # convert to longformat
  c <- e[which(e[,1]==1),]                         # only keep 1s
  d <- do.call(rbind,(strsplit(paste(c[,2],c[,3],sep="-"),"-"))) # ji combination
  g <- do.call(rbind,(strsplit(paste(c[,3],c[,2],sep="-"),"-"))) # ij combination
  f <- rbind(d,g) # both combinations joint
  w <- matrix(ncol=S,nrow=S,0) # initialize new matrix
  
  if(type==1){
    for (i in 1:dim(f)[1]) {w[as.numeric(f[i,1]),as.numeric(f[i,2])]=1} # writes 1s where 1s have previously been
    w[which(w == 1)]<- runif(length(which(w ==1)),min=-0.7,max=0.7) # substitutes them with runif value
    diag(w) <- 0
    b<-w
  }
  
  if (type==2){
    for (i in 1:dim(f)[1]) {w[as.numeric(f[i,1]),as.numeric(f[i,2])] =1}
    w[which(w == 1)]<-runif(length(which(w ==1)),min=-0.9,max=0.9)
    diag(w) <- 0
    
    b <- upper.tri(w)*w # only fills upper triagle of the matrix, rest is set to 0
    f <- -t(b) # transpose so that only lower triangle is filled with negative values from upper triangle
    b[lower.tri(b)]<-f[lower.tri(f)] # lower triangle now mirrors upper triangle but with neg values
    
    # substitute first row and col with new runif draws
    F1=runif(S,-0.9,0.9)
    b[1,] = F1; b[,1] = F1*-1
    b[1,1] = 0
  }
  
  if (type==3){
    for (i in 1:dim(f)[1]) {w[as.numeric(f[i,1]),as.numeric(f[i,2])]=1}
    w[which(w==1)]<- runif(length(which(w==1)),min=0,max=0.9)
    diag(w) <- 0
    
    b <- upper.tri(w)*w
    f <- -t(b)
    b[lower.tri(b)]<-f[lower.tri(f)] # until here same as type 2
    b <- b*-1 # make positive values negative and vice versa
  }
  
  ### Predation
  if (type==4){
    for (i in 1:nrow(f)) {w[as.numeric(f[i,1]),as.numeric(f[i,2])]=1}
    w[which(w == 1)] <- runif(length(which(w ==1)),min=0.05,max=0.9)
    diag(w) <- 0
    
    b <- upper.tri(w)*w
    f <- -t(b)
    b[lower.tri(b)]<-f[lower.tri(f)]
    b <- b*-1 # until here same as type 3
    
    # generate vector where half of the species has a neg value btw. -0.9 and 0
    F1 <- rep(0,S)
    nb <-  runif(round(S/2),-0.9,0)
    F1[sample(1:S,length(nb))] <- nb
    # use vector for first col and row
    b[1,] <- F1; b[,1] <- F1*-1
    b[1,1] <- 0
    b <- round(b,2)
    
    vecta <- b[lower.tri(b)]
    vecta[vecta>0] <- vecta[vecta>0]-0.05
    b[lower.tri(b)] <- vecta
    
  }
  
  # mutualism
  if (type==5){
    for (i in 1:nrow(f)) {w[as.numeric(f[i,1]),as.numeric(f[i,2])]=1}
    w[which(w==1)] <- runif(length(which(w==1)),min=0,max=0.4)
    diag(w) <- 0
    
    b <- upper.tri(w)*w
    d <- b[upper.tri(b)]>0
    f <- runif(length(which(d==TRUE)),min=0,max=0.4)
    d[d==TRUE] <- f
    b[lower.tri(b)] <- d  # fill corresponding positive positions in upper tri. with
    # new runif draw (same interval) for lower triangle
    
    F1 <- F2 <- rep(0,S)
    nb <- runif(round(S/2),0,0.4)
    F1[sample(1:S,length(nb))] <- nb
    
    nb <-  runif(round(S/2),0,0.4)
    F2[F1>0] <- nb
    
    b[1,] <- F1; b[,1] <- F2
    b[1,1] <- 0
    b <- round(b,2)
  }
  
  # competition
  if (type==6){ # same as type 5 just with negative values
    for (i in 1:nrow(f)) {w[as.numeric(f[i,1]),as.numeric(f[i,2])]=1}
    w[which(w==1)] <- runif(length(which(w ==1)),min=-0.6,max=0)
    diag(w) <- 0
    
    b <- upper.tri(w)*w
    d <- b[upper.tri(b)]<0
    f <- runif(length(which(d==TRUE)),min=-0.6,max=0)
    d[d==TRUE] <- f
    b[lower.tri(b)]<- d
    
    F1 <- F2 <- rep(0,S)
    nb <- runif(round(S/2),-0.6,0)
    F1[sample(1:S,length(nb))] <- nb
    
    nb <- runif(round(S/2),-0.6,0)
    F2[F1<0] <- nb
    
    b[1,] <- F1; b[,1] <- F2
    b[1,1] <- 0
    b <- round(b,2)
    
  } # end of if
  
  if(export==T){
    write.table(round(b,2),file=name_file,sep=" ",dec=".",row.names = FALSE,col.names = FALSE)
  }
  return(b)
}
