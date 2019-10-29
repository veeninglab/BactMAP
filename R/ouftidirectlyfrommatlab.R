##7/19/2017

##import functions when importing right from Matlab file.

##Renske van Raaphorst

##part of shinyspots package

############################################################

############################################################

##extract cellList as a data frame with lists of data for the mesh, spots, etc.
extr_OuftiCellList <- function(matfile){
  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    inp <- readline("Package 'R.matlab' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){utils::install.packages("R.matlab")}else{stop("Canceled")}
  }
  matlist <- R.matlab::readMat(matfile)
  matcellList <- matlist$cellList[[1]]
  matcelnums <- matlist$cellList[[2]]
  cellList <- list()
  for(n in 1:length(matcellList)){
    cellListN <- list()
    if(length(matcellList[[n]][[1]])!=0){
      for(y in 1:length(matcellList[[n]][[1]])){
        cell <- matcelnums[[n]][[1]][[y]]
        celldatas <- as.data.frame(matcellList[[n]][[1]][[y]][[1]])
        celldatas <- as.data.frame(t(celldatas))
        celldatas$cell <- cell
        celldatas$frame <- n
        cellListN[[y]] <- celldatas
      }

    }
    cellList[[n]] <- do.call(plyr::rbind.fill, cellListN)
  }
  cellList <- do.call(plyr::rbind.fill, cellList)
  return(cellList)
}

#make seperate data frame of the mesh - combined with parts of the meshturn function
extr_OuftiMat <- function(cellList){
  for(n in 1:nrow(cellList)){
    meshline <- as.data.frame(cellList$mesh[n])
    if(nrow(meshline)>1){
      colnames(meshline) <- c("x0", "y0", "x1", "y1")
      cell <- cellList$cell[n]
      meshline$cell <- cell
      frame <- cellList$frame[n]
      meshline$frame <- frame
      meshline$num <- 0:(nrow(meshline)-1)
      xdist <- meshline$x0 - meshline$x1
      ydist <- meshline$y0 - meshline$y1
      widths <- sqrt(xdist^2 + ydist^2)
      max_width <- max(widths)
      meshline$max.width <- max_width
      for(z in 1:nrow(meshline)){
        if(z==1){
          meshline$length[z] <- 0
          meshline$steplength[z] <- 0
        }
        if(z > 1){
          xdistL0 <- meshline$x0[z]-meshline$x0[z-1]
          ydistL0 <- meshline$y0[z]-meshline$y0[z-1]
          distL0 <- sqrt(xdistL0^2 + ydistL0^2)
          angL0 <- atan(ydistL0/xdistL0)
          angd <- atan(ydist[z]/xdist[z])
          anglength <- pi-abs(angL0)-abs(angd)
          steplength <- sin(anglength)*distL0
          meshline$steplength[z] <- steplength
          meshline$length[z] <- sum(meshline$steplength[1:z])
        }
        if(z==max(meshline$num+1)){
          angK <- atan(xdist[z-1]/ydist[z-1])
          anglengthEnd <- (pi/2)-abs(angL0-abs(angK))
          steplengthK <- sin(anglengthEnd)*distL0
          meshline$steplength[z] <- steplengthK
          meshline$length[z] <- sum(meshline$steplength[1:z])

        }
      }
      meshline$max.length <- sum(meshline$steplength)
      if(n==1){MESH <- meshline}
      if(n>1){MESH <- rbind(MESH, meshline)}
    }
  }
  return(MESH)
}

extr_OuftiSpots <- function(cellList){
  spots <- cellList$spots
  u <- 1
  for(n in 1:length(spots)){
    l <- spots[[n]][1]
    d <- spots[[n]][2]
    x <- spots[[n]][3]
    y <- spots[[n]][4]
    position <- spots[[n]][5]
    adj_Rsquared <- spots[[n]][6]
    CI_xy <- spots[[n]][7]
    spotframe <- data.frame("l"=t(l[[1]]), "d"=t(d[[1]]), "x"=t(x[[1]]), "y"=t(y[[1]]), "position" = t(position[[1]]), "adj_Rsquared" = t(adj_Rsquared[[1]]), "CI_xy" =paste(CI_xy[[1]]))
    if(nrow(spotframe)>0){
      spotframe$frame <- cellList$frame[n]
      spotframe$cell <- cellList$cell[n]
      if(u!=1){
        spottotal <- rbind(spottotal, spotframe)
      }
      else{spottotal <- spotframe
      u <- u+1}
    }
  }
  return(spottotal)
}

#' @export
extr_Oufti <- function(matfile, mag="No_PixelCorrection", phylo=FALSE){
  CSV <- FALSE
  ##if CSV input; take cellList, matlab file & object file from there
  if(substr(matfile, nchar(matfile)-3, nchar(matfile))==".csv"|substr(matfile, nchar(matfile)-3, nchar(matfile))==".txt"){
    CSV <- TRUE
    outlist <- extr_OuftiCSV(matfile)
    Mesh <- outlist$mesh
    cellList <- outlist$cellList
    if("objectframe"%in%names(outlist)){
      OBJ <- outlist$objectframe
    }
  }
  ##original matlab file
  else{
    outlist <- list()
    message("Extracting original data from the matlab file... This step may take a while.")
    cellList <- extr_OuftiCellList(matfile)
    message("Finished extracting.")
    outlist$cellList <- cellList
    message("Taking the x/y coordinates of the mesh outlines of each cell... This step may also take a while.")
    Mesh <- extr_OuftiMat(cellList)
    message("Converting mesh file into standard BactMAP format...")
  }
  ##turn meshes to one x/y column
  if("signal1"%in%colnames(outlist$cellList)){
    if(length(unique(outlist$signal1))>1){
      outlist$cellList$mean.signal <- unlist(lapply(outlist$cellList$signal1, function(x) mean(x)))
      outlist$cellList$sd.signal <- unlist(lapply(outlist$cellList$signal1, function(x) stats::sd(x)))
    }
  }
  Mesh <- spotrXYMESH(Mesh)
  Mesh <- meshTurn(Mesh)
  ##pixel --> um
  outlist$pixel2um <- unlist(get(magnificationList, envir=magEnv)[mag])
  Mesh$max_um <- Mesh$max.length*outlist$pixel2um
  Mesh$maxwum <- Mesh$max.width*outlist$pixel2um
  Mesh$Xrot_micron <- Mesh$X_rot*outlist$pixel2um
  Mesh$Yrot_micron <- Mesh$Y_rot*outlist$pixel2um
  if("signal1"%in%colnames(outlist$cellList)){
    if(length(unique(outlist$signal1))>1){
      Mesh <- merge(Mesh, outlist$cellList[,c("cell", "frame", "mean.signal", "sd.signal")])
      meansignalList <- unique(Mesh[,c("cell", "frame", "max.width", "maxwum", "max.length", "max_um", "mean.signal", "sd.signal")])
      outlist$meansignalList <- meansignalList
    }
  }
  outlist$mesh <- Mesh
  if("objectframe"%in%names(outlist)){
    OM <- suppressWarnings(centrefun(OBJ))
    OM <- suppressWarnings(midobject(Mesh, OM, outlist$pixel2um))
    outlist$object_relative <- OM
  }
  ##then take the spots

  if(length(unique(cellList$spots))>1){
    if(CSV==FALSE){
      message("Taking the spot coordinates and information per cell...")
      spotframe <- extr_OuftiSpots(cellList)
      outlist$spotframe <- spotframe
    }
    message("Adding cell dimensions (length/width) to spot information...")
    if(missing(mag)){
      mag <- "No_PixelCorrection"
    }
    spot_mesh <- mergeframes(outlist$spotframe, Mesh, mag, ouf=TRUE)

    ##here are the spots put relative to mesh length/width
    outlist$spots_relative <- spot_mesh[!is.na(spot_mesh$cell),]
  }

  ##optional (default = OFF) add genealogy information as phylo objects.
  if(phylo == TRUE){
    if(length(cellList$descendants)>1){
      message("Getting phylogenies from ancestor/descendants information...")
      cellList$parent <- as.numeric(lapply(cellList$ancestors, function(x) utils::tail(as.numeric(x), n=1)))
      cellList$parent[is.na(cellList$parent)] <- 0
      cellList$child1 <- as.numeric(lapply(cellList$descendants, function(x) as.numeric(x)[1]))
      cellList$child2 <- as.numeric(lapply(cellList$descendants, function(x) as.numeric(x)[2]))
      cellList$death <- unlist(lapply(cellList$cell, function(x) max(cellList$frame[cellList$cell==x])))
      cellList$birthframe <- as.numeric(as.character(cellList$birthframe))
      cellList$edgelength <- cellList$death-cellList$birthframe
      cellList$root <- 0
      CL <- cellList[cellList$frame==cellList$death,]
      phylolist <- getphylolist_SupSeg(CL)
      if(ape::is.rooted.phylo(phylolist$generation_lists)!=T){
        u <-lapply(phylolist$generation_dataframes, function(x) is.data.frame(x))
        phylolist$generation_dataframes <- phylolist$generation_dataframes[c(1:length(u))[u==T]]
        phylolist$generation_lists <- phylolist$generation_lists[c(1:length(u))[u==T]]
      }
      if(length(cellList$spots)>1){
        message("Saving the relative spot localizations per phylogeny...")
        if(is.data.frame(phylolist$generation_dataframes)==FALSE){
          phlist_allcells <- lapply(phylolist$generation_dataframes, function(x) x[,c("node", "cell", "frame", "parent", "child1", "child2", "death", "birthframe", "edgelength", "root", "nodelabel")])
          MRTlist <- lapply(phlist_allcells, function(x) merge(x, spot_mesh[,c("cell","frame", "max.length", "max.width", "spot", "totalspot", "Lmid", "pole1", "pole2", "Dum")]), all=T)
          phylolist$spot_relative_list <- MRTlist
        }
        if(is.data.frame(phylolist$generation_dataframes)==TRUE){
          ph <- phylolist$generation_dataframes[,c("node", "cell", "frame", "parent", "child1", "child2", "death", "birthframe", "edgelength", "root", "nodelabel")]
          phylolist <- merge(ph, spot_mesh[,,c("cell","frame", "max.length", "max.width", "spot", "totalspot", "Lmid", "pole1", "pole2", "Dum") ], all=T)
        }
      }
      message("Saving cell outlines per phylogeny...")
      if(is.data.frame(phylolist$generation_dataframes)==FALSE){
        Mlist <- lapply(phylolist$generation_dataframes, function(x) merge(x[,c("node", "cell", "frame", "parent", "child1", "child2", "death", "birthframe", "edgelength", "root", "nodelabel")], Mesh))
      }
      if(is.data.frame(phylolist$generation_dataframes)==TRUE){
        Mlist <- merge(phylolist$generation_dataframes[,c("node", "cell", "frame", "parent", "child1", "child2", "death", "birthframe", "edgelength", "root", "nodelabel")], Mesh)
      }
      phylolist$meshdata <- Mlist
      outlist$timelapsedata <- phylolist
    }
    if(length(cellList$descendants)<1){
      warning("No genealogy information found. Please check your original data file.")
    }
  }
  message("Finished Oufti extraction.")
  return(outlist)
}

