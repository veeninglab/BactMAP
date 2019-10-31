##Supersegger

##3.08.2018
##Renske van Raaphorst


##Combination of functions to: * get mesh outlines from cell data * make it compatible with bactmap format


#function to get all cell files into R as a list per frame number
upperdir <- function(frameN, loc){
  cellnameList <- list.files(path=paste(loc, "\\xy", frameN, "\\cell", sep=""))
  matlist <- lapply(cellnameList, function(x) R.matlab::readMat(paste(loc, "\\xy", frameN, "\\cell\\", x, sep="")))
  return(matlist)
}

#getting a list of cell files from frame 0:frames-1
readallsegcells <- function(frames, loc, startframe){
  allcells <- lapply(c(startframe:(frames+startframe-1)), function(x) upperdir(x, loc))
  return(allcells)
}

#getting the cell mask out of the cell lists
celmask <- function(matf, n, b){
  ras <- matf$CellA[[n]][[1]][[4]]
  ras[ras==0] <- NA
  ras <- raster::raster(ras,
                        xmn=min(matf$CellA[[n]][[1]][[2]]),
                        xmx = max(matf$CellA[[n]][[1]][[2]]),
                        ymn = min(matf$CellA[[n]][[1]][[3]]),
                        ymx = max(matf$CellA[[n]][[1]][[3]])
  )
  rasp <- raster::rasterToPolygons(ras, dissolve=TRUE)
  rasp <- as.data.frame(rasp@polygons[[1]]@Polygons[[1]]@coords)
  rasp$frame <- (n-1+as.numeric(b))
  return(rasp)
}

getallmaskstime <- function(matf){
  if((matf$death-matf$birth)>0){
    rangeC <- c(1:length(matf$CellA))
  }
  if((matf$death==matf$birth)){
    rangeC <- 1
  }
  celmasklist <- lapply(rangeC, function(x) celmask(matf, x, matf$birth))
  return(celmasklist)
}

#getting a list of length x (amount of cells) containing the cell mask for each cell (y)
getallmasks <- function(matlist){
  matlength <- length(matlist)
  celmasklist <- lapply(1:matlength, function(x) lapply(matlist[[x]], function(y) getallmaskstime(y)))
  return(celmasklist)
}


#turn the dataframes
flipallcells <- function(celllists){
  cellflip <- lapply(1:length(celllists), function(x) lapply(celllists[[x]], function(z) lapply(z$CellA, function(y)  as.data.frame(t(as.data.frame(y[[1]]))))))
  return(cellflip)
}

#add frame number
addcellnumframe <- function(cellflip, x){
  cellflip$frame <- x
  return(cellflip)
}

#cell number
addcellnumcell <- function(cellflip, y){
  cellflip$cell <- y
  return(cellflip)
}

#final combination function

bindallcellsandmeshes <- function(cellflip, cellmask, timelapse=TRUE, cellListIn){
  cellflipout <- list()
  for(n in 1:length(cellflip)){
    for(u in 1:length(cellflip[[n]])){
      if(length(cellmask[[n]][[u]])==1){cellflip[[n]][[u]][[2]]<-NULL}
      for(z in 1:length(cellmask[[n]][[u]])){
        cellflip[[n]][[u]][[z]]$cell <- u
        cellmask[[n]][[u]][[z]]$cell <- u
        cellflip[[n]][[u]][[z]]$cellmask <- list(cellmask[[n]][[u]][[z]])
        cellmask[[n]][[u]][[z]]$num <- c(1:nrow(cellmask[[n]][[u]][[z]]))
        cellflip[[n]][[u]][[z]]$frame <- unique(cellmask[[n]][[u]][[z]]$frame)
      }
      cellflip[[n]][[u]] <- do.call('rbind', cellflip[[n]][[u]])
      cellmask[[n]][[u]] <- do.call('rbind', cellmask[[n]][[u]])
    }
    cellflipframe <- do.call('rbind', cellflip[[n]])
    cellmaskframe <- do.call('rbind', cellmask[[n]])
    if(timelapse==TRUE){

      if("frame"%in%colnames(cellflipframe)){
        cellflipframe$location <- n
      }
      if(!"frame"%in%colnames(cellflipframe)){
        cellflipframe$frame <- n
      }

      if("frame"%in%colnames(cellmaskframe)){
        cellmaskframe$location <- n
      }
      if(!"frame"%in%colnames(cellmaskframe)){
        cellmaskframe$frame <- n
      }
    }
    if(timelapse==FALSE){
      cellmaskframe$frame <- n
      cellflipframe$frame <- n
    }

    if(n==1){
      if(cellListIn==TRUE){
        cellflipout$cellList <- cellflipframe
      }
      cellflipout$mesh <- cellmaskframe
    }
    else{
      if(cellListIn==TRUE){
        cellflipout$cellList <- rbind(cellflipout$cellList, cellflipframe)
      }
      cellflipout$mesh <- rbind(cellflipout$mesh, cellmaskframe)
    }
  }
  return(cellflipout)
}


#####################################################
##extract function for exporting.

#' @export
extr_SuperSeggerCells <- function(loc, frames, mag, timelapse=FALSE, startframe=0, cellList=FALSE){
  if(paste(loc, "/xy", startframe, sep="")%in%list.dirs(loc)==FALSE){
    stop(paste("Cannot find SuperSegger output folder(s) starting with 'xy' in the directory '", loc, "'. Please make sure to set the correct path in variable 'loc'", sep=""))
  }
  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    inp <- readline("Package 'R.matlab' and 'raster' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(!requireNamespace("rgeos", quietly=TRUE)){
      message("Installing 'raster' dependency 'rgeos'..")
      utils::install.packages("rgeos")
    }
    if(inp=="y"|inp=="Y"){utils::install.packages(c("R.matlab", "raster"))}else{stop("Canceled")}
  }
  if (!requireNamespace("raster", quietly = TRUE)) {
    inp <- readline("Package 'raster' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){utils::install.packages("raster")}else{stop("Canceled")}
  }

  segcells <- readallsegcells(frames=frames, loc=loc, startframe = startframe)
  segmasks <- getallmasks(segcells)
  finalframe <- bindallcellsandmeshes(flipallcells(segcells), segmasks, timelapse, cellListIn=cellList)
  finalframe$mesh <- meshTurn(finalframe$mesh, "x", "y")
  finalframe$mesh$Y <- finalframe$mesh$Ymid + (finalframe$mesh$Ymid - finalframe$mesh$Y)
  finalframe$mesh$Y_rot <- -finalframe$mesh$Y_rot
  if(!missing(mag)){
    if(is.numeric(unlist(get(magnificationList,envir=magEnv)[mag]))==FALSE){
      stop("Magnification conversion factor not recognized. Please use addPixels2um('pixelName', pixelsize) to add your conversion factor")
    }
    finalframe$mesh$Yrot_micron <- finalframe$mesh$Y_rot * unlist(get(magnificationList,envir=magEnv)[mag])
    finalframe$mesh$Xrot_micron <- finalframe$mesh$X_rot * unlist(get(magnificationList,envir=magEnv)[mag])
    finalframe$mesh$max_um <- finalframe$mesh$max.length* unlist(get(magnificationList, envir=magEnv)[mag])
    finalframe$mesh$maxwum <- finalframe$mesh$max.width *  unlist(get(magnificationList, envir=magEnv)[mag])
    finalframe$mesh$area_um <- finalframe$mesh$area *  unlist(get(magnificationList, envir=magEnv)[mag])^2
    finalframe$pixel2um <- unlist(get(magnificationList,envir=magEnv)[mag])
  }
  if(missing(mag)){
    finalframe$pixel2um <- c("No_PixelCorrection" = 1)
  }
  return(finalframe)
}
