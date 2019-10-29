##collection of raw image plotting functions for BactMAP

#upload your image stack (one color)
#' @export
extr_OriginalStack <- function(picloc){
  if (!requireNamespace("tiff", quietly = TRUE)) {
    if(!requireNamespace("raster", quietly = TRUE)){
      inp <- readline("Package 'tiff' and 'raster' needed for this function to work. Press 'y' to install them, or any other key to cancel.")
      if(inp=="y"|inp=="Y"){utils::install.packages(c("tiff", "raster"))}else{stop("Canceled")}
    }else{
      inp <- readline("Package 'tiff' needed for this function to work. Press 'y' to install it, or any other key to cancel.")
      if(inp=="y"|inp=="Y"){utils::install.packages("tiff")}else{stop("Canceled")}}
  }
  if(!requireNamespace("raster", quietly = TRUE)){
    inp <- readline("Package 'raster' needed for this function to work. Press 'y' to install it, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){utils::install.packages("raster")}else{stop("Canceled")}
  }
  suppressWarnings(im <- tiff::readTIFF(picloc, all=T)) #if you want the best resolution, it needs to be a .tiff file
  im <- lapply(im, function(x) raster::raster(x))
  imdatframe <- lapply(im, function(x) as.data.frame(methods::as(x, "SpatialPixelsDataFrame"))) #get values
  imdatframe <- lapply(imdatframe, function(x) changecols(x))
  nx <- length(unique(imdatframe[[1]]$x))
  ny <- length(unique(imdatframe[[1]]$y))
  imdatframe <- lapply(imdatframe, function(x) changeres(x, nx, ny))
  return(imdatframe)
}


#' @export
extr_OriginalCells <- function(imdatframe, mesh, surroundings=FALSE, turnCell=TRUE){
  if("area"%in%colnames(mesh)){
    mesh <- mesh[mesh$area>2,]
  }

  allcellslist <- lapply(unique(mesh$frame), function(x) pipperframe(imdatframe, mesh, x, surroundings))
    allcellsframe <- do.call(rbind, allcellslist)
    allcellsframe$cell <- allcellsframe$pip
    allcellsframe$pip <- NULL
    if(turnCell==TRUE){
      outlist <- meshTurn(mesh, rawdatafile=allcellsframe)
    }else{outlist <- allcellsframe}
  return(outlist)
}


#' @export
plotCellsTime <- function(celdat,
                           updown = T,
                           movie = F,
                           viridisoption = "magma",
                           cellN,
                           minf,
                           maxf,
                           outlines=FALSE,
                           meshdata#,
                           #overlay=FALSE
                          ){

  if(movie==TRUE){
    if(!requireNamespace("gganimate", quietly = TRUE)){
      inp <- readline("Package 'gganimate' needed to make an animation. Press 'y' to install it, or any other key to cancel.")
      if(inp=="y"|inp=="Y"){utils::install.packages("gganimate")}else{stop("Canceled")}
    }
  }
  if(missing(minf)){
    minf <- min(celdat$frame)
  }
  if(missing(maxf)){
    maxf <- max(celdat$frame)
  }

  if(outlines==TRUE&!missing(meshdata)){
    meshdata <- meshdata[,c("frame", "cell", "X_rot", "Y_rot", "num")]
  }
  if(outlines==TRUE&missing(meshdata)){
    stop("Outlines is set to TRUE but no meshdata was found. Please specify your mesh dataframe.")
  }
  if(outlines==FALSE&!missing(meshdata)){
    warning("Meshdata was specified, but outlines are set to FALSE. No outlines will be drawn.")
  }
  #when no cell number is indicated:return a list of plots/movie objects
  if(missing(cellN)){
    plotout <- lapply( unique(celdat$cell), function(x) plotcellsframelist(celdat[celdat$cell==x,], maxframes=maxf, minframes=minf, updown, movie, viridisoption, outlines, meshdata#, overlay
                                                                           ) + ggplot2::ggtitle(x))
  }
  #and just plot one plot if cellN exists
  if(missing(cellN)!=T){
    if(length(cellN)==1){
      plotout <- plotcellsframelist(celdat[celdat$cell==cellN,], maxf, minf, updown, movie, viridisoption, outlines, meshdata#, overlay
                                    ) + ggplot2::ggtitle(cellN)
    }
    if(length(cellN)>1){
      plotout <- lapply(cellN, function(x) plotcellsframelist(celdat[celdat$cell==x,], maxf, minf,updown,movie,viridisoption, outlines, meshdata#, overlay
                                                              ) + ggplot2::ggtitle(x))
    }
  }
  return(plotout)
}

###################################|(*-*)/#############################################################
changecols <- function(x){
  colnames(x) <- c("values" , "x", "y")
  return(x)
}
changeres <- function(dat, nx, ny){
  dat$x <- dat$x*nx - 0.5
  dat$y <- max(dat$y)-dat$y
  dat$y <- dat$y*ny + 1
  return(dat)
}



pinping <- function(dat, mesh, x, surroundings=FALSE){
  message(paste("Cell", x))
  mesh <- mesh[mesh$cell==x,]

  if(surroundings==FALSE){
  minmeshx <- min(mesh$X)-2
  minmeshy <- min(mesh$Y)-2
  maxmeshx <- max(mesh$X)+2
  maxmeshy <- max(mesh$Y)+2
  dat <- dat[dat$x>minmeshx&dat$y>minmeshy&dat$x<maxmeshx&dat$y<maxmeshy,]
  p <- SDMTools::pnt.in.poly(dat[,c("x","y")], mesh[mesh$cell==x,][,c("X","Y")])
  #if pip == 1, the point is inside the polygon. if p==0, it is not.
  #I replace the pips which are 1 with the cell number x
  p$pip[p$pip!=0] <- x
  datje <- merge(dat, p[p$pip!=0,])
  }
  if(surroundings==TRUE){
  mesh <- unique(mesh[,c("Xmid","Ymid", "max.length", "max.width", "angle")])
  w <- 2
  ybox1 <- (0.5*mesh$max.length + w)*sin(pi-mesh$angle)-(0.5*mesh$max.width+w)*cos(pi-mesh$angle) + mesh$Ymid
  xbox1 <- (0.5*mesh$max.length + w)*cos(pi-mesh$angle)+(0.5*mesh$max.width+w)*sin(pi-mesh$angle)+ mesh$Xmid
  xbox2 <- (0.5*mesh$max.length+w)*cos(pi-mesh$angle)-(0.5*mesh$max.width+w)*cos(pi-mesh$angle) + mesh$Xmid
  ybox2 <- (0.5*mesh$max.length+w)*sin(pi-mesh$angle)+(0.5*mesh$max.width+w)*sin(pi-mesh$angle) + mesh$Ymid
  xbox4 <- mesh$Xmid-(xbox2-mesh$Xmid)
  ybox4<- mesh$Ymid-(ybox2-mesh$Ymid)
  xbox3 <- mesh$Xmid-(xbox1-mesh$Xmid)
  ybox3 <- mesh$Ymid-(ybox1-mesh$Ymid)
  box <- data.frame("x"=c(xbox1,xbox2,xbox3,xbox4), "y"=c(ybox1, ybox2, ybox3,ybox4))
  dat <- dat[dat$x>min(box$x)&dat$y>min(box$y)&dat$x<max(box$x)&dat$y<max(box$y),]
  p <- SDMTools::pnt.in.poly(dat[,c("x","y")], box)

  #add cell number (as pip so it matches the situation with the cell outlines)
  p$pip[p$pip!=0] <- x
  datje <- merge(dat, p[p$pip!=0,])
  }
  return(datje)
}

pipperframe <- function(dat, mesh, y, surroundings=FALSE){
  mesh <- mesh[mesh$frame==y,]
  dat <- dat[[y]]
  message(paste("Finding & saving the raw data per cell for frame", y))
  datjeslist <- lapply(unique(mesh$cell), function(x) pinping(dat, mesh, x, surroundings))
  datjesframe <- do.call(rbind, datjeslist)
  datjesframe$frame <- y
  return(datjesframe)
}

##################Plotting cells in a tower/row/movie per cell, per frame

plotcellsframelist <- function(TRframe, maxframes, minframes, updown=F, movie=F, viridisoption="magma", outlines=FALSE,meshdata #, overlay=FALSE
                               ){
  nframes <- length(TRframe$frame)
  if(missing(minframes)){
    minframes <- min(TRframe$frame)
  }
  if(missing(maxframes)){
    maxframes <- max(TRframe$frame)
  }
  if(nframes>0){

    minTRframe <- min(TRframe$frame)
    maxTRframe <- max(TRframe$frame)
    if(minTRframe>minframes){
      minframes<-minTRframe
    }
    if(maxframes>maxTRframe){
      maxframes <- maxTRframe
    }

    TRframe <- TRframe[TRframe$frame>=minframes&TRframe$frame<=maxframes,]

    if(!missing(meshdata)&outlines==TRUE){
      meshdata <- meshdata[meshdata$frame>=minframes&meshdata$frame<=maxframes,]
    }

    p <- ggplot2::ggplot(TRframe) +
      ggplot2::geom_polygon(ggplot2::aes_string(x='xt',y='yt',fill='values',group='pointN'),color=NA) +
      ggimage::theme_transparent() + ggplot2::coord_fixed() +
      ggplot2::scale_fill_viridis_c(option=viridisoption) +
      ggplot2::xlim(c(min(TRframe$xt, na.rm=TRUE), max(TRframe$xt, na.rm=TRUE))) +
      ggplot2::ylim(c(min(TRframe$yt, na.rm=TRUE), max(TRframe$yt, na.rm=TRUE))) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(NULL) +
      ggplot2::theme(strip.background=ggplot2::element_rect(fill="transparent", colour=NA),
            strip.text=ggplot2::element_blank(),
            legend.position="none",
            axis.text=ggplot2::element_blank(),
            axis.ticks=ggplot2::element_blank(),
            plot.margin=ggplot2::unit(c(0,0,0,0),"mm"),
            panel.spacing = ggplot2::unit(0, "mm"),
            panel.background=ggplot2::element_rect(fill="transparent", colour=NA),
            plot.background=ggplot2::element_rect(fill="transparent", colour=NA),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank())

    #if(overlay==TRUE){
      #p <- p + ggplot2::geom_path(data=meshdata, ggplot2::aes(x=X_rot,y=Y_rot), color="white")
    #}

    if(updown==T&movie==F){

      p <- p + ggplot2::facet_grid(frame~.)
    }

    if(updown==F&movie==F){
      p <- p + ggplot2::facet_grid(~frame)
    }

    if(movie==T){
      frame <-NULL
      p <- p + gganimate::transition_manual(frame)
    }
    return(p)
  }

  if(nframes==0){
    return(NA)
  }

}

###plot raw image

#' @export
plotRaw <- function(tiffdata,
                    meshdata,
                    spotdata,
                    frameN=1,
                    xrange,
                    yrange,
                    viridisoption="inferno",
                    meshcolor="white",
                    spotcolor="yellow",
                    valuerange,
                    legend=FALSE){

  if(missing(valuerange)!=T&missing(tiffdata)!=T){
    tiffdata[[frameN]] <- tiffdata[[frameN]][tiffdata[[frameN]]$value>valuerange[1]&tiffdata[[frameN]]$value<valuerange[2],]
  }
  if(legend==TRUE){
    lpos="right"
  }
  if(legend==FALSE){
    lpos="none"
  }
  if(missing(tiffdata)!=T){
  plotcells <- ggplot2::ggplot(tiffdata[[frameN]]) + #plot raw image
    ggplot2::geom_raster(ggplot2::aes_string(x='x',y='y',fill='values')) + #use geom_raster to remake image out of dataframe
    ggplot2::theme_classic() + #simple theme, no backgrounds
    ggplot2::scale_fill_viridis_c(option=viridisoption) + #well-working color scheme for gradient values
    ggplot2::theme(legend.position=lpos) # remove legend for easy viewing
  }
  if(missing(tiffdata)==T){
    plotcells <- ggplot2::ggplot() + ggplot2::theme_dark()
  }
  #add x and/or y range + fixed coordinates when indicated:
  if(missing(xrange)!=T&missing(yrange)!=T){
    plotcells <- plotcells + ggplot2::coord_fixed(xlim=xrange, ylim=yrange)  #sub-set of the image frame to zoom in
  }
  if(missing(xrange)!=T&missing(yrange)==T){
    plotcells <- plotcells + ggplot2::coord_fixed(xlim=xrange)
  }
  if(missing(xrange)==T&missing(yrange)!=T){
    plotcells <- plotcells + ggplot2::coord_fixed(ylim=yrange)
  }
  if(missing(xrange)==T&missing(yrange)==T){
    plotcells <- plotcells + ggplot2::coord_fixed()
  }


  #add mesh data when given:
  if(missing(meshdata)!=T){
    meshdata <- meshdata[order(meshdata$frame, meshdata$cell, meshdata$num),]
    plotcells <- plotcells + #plot made above
      ggplot2::geom_path(data=meshdata[meshdata$frame==frameN,], ggplot2::aes_string(x='X',y='Y', group='cell'), color=meshcolor) #add outline of cells, only frame one, white color
  }
  if(missing(spotdata)!=T){
    plotcells <- plotcells +
      ggplot2::geom_point(data=spotdata[spotdata$frame==frameN,], ggplot2::aes_string(x='x',y='y'), shape=1, color=spotcolor)# add yellow empty dots of our spot localizations on top
  }
  return(plotcells)
}


########plot time series (or by cell length) raw data kymographs
#' @export

prepForKymo <- function(turnedCells, dimension="length", bins = 25, sizeAV=FALSE){

  turnedCells <- unique(turnedCells[,c("values", "frame", "cell", "X_rot","Y_rot", "max.length", "max.width", "pointN")])
  originalCells <- list()
  for(n in unique(turnedCells$frame)){

    U <- turnedCells[turnedCells$frame==n,]
    if(dimension=="length"){
      U$group <- unlist(lapply(unique(U$cell), function(x) cut(U$X_rot[U$cell==x], breaks=bins, labels= c(1:bins))))
    }
    if(dimension=="width"){
      U$group <- unlist(lapply(unique(U$cell), function(x) cut(U$Y_rot[U$cell==x], breaks=bins, labels= c(1:bins))))
    }

    Umeans <- suppressWarnings(lapply(unique(U$cell), function(x) stats::aggregate(U[U$cell==x,], by=list(U$group[U$cell==x]), FUN=mean, na.rm=T)))

    Umeansall <- do.call('rbind', Umeans)


    if(sizeAV==TRUE){

      Umeansall$a <- Umeansall$frame + 0.5
      Umeansall$b <- Umeansall$frame - 0.5
      if(dimension=="length"){
        Umeansall$d <- unlist(lapply(unique(Umeansall$cell), function(x) Umeansall$X_rot[Umeansall$cell==x] + (max(Umeansall$X_rot[Umeansall$cell==x])/bins)*1.15))
        Umeansall$c <- unlist(lapply(unique(Umeansall$cell), function(x) Umeansall$X_rot[Umeansall$cell==x] - (max(Umeansall$X_rot[Umeansall$cell==x])/bins)*1.15))
      }
      if(dimension=="width"){
        Umeansall$d <- unlist(lapply(unique(Umeansall$cell), function(x) Umeansall$Y_rot[Umeansall$cell==x] + (max(Umeansall$Y_rot[Umeansall$cell==x])/bins)*1.15))
        Umeansall$c <- unlist(lapply(unique(Umeansall$cell), function(x) Umeansall$Y_rot[Umeansall$cell==x] - (max(Umeansall$Y_rot[Umeansall$cell==x])/bins)*1.15))
      }


      Umeansall <- reshape2::melt(Umeansall, measure.vars=c("a", "b"), value.name="x_coords")
      Umeansall$frameh <- Umeansall$variable
      Umeansall$variable <- NULL
      Umeansall <- reshape2::melt(Umeansall, measure.vars=c("d", "c"), value.name="y_coords")
      Umeansall$frameh <- as.character(Umeansall$frameh)
      Umeansall$frameh[Umeansall$frameh=="b"&Umeansall$variable == "c"] <- "c"
      Umeansall$frameh[Umeansall$frameh=="a"&Umeansall$variable == "c"] <- "d"
      Umeansall$variable <- NULL
      Umeansall <-  Umeansall[order( Umeansall$frame,  Umeansall$cell,  Umeansall$Group.1,  Umeansall$frameh),]

    }

    Umeansall$group <- Umeansall$Group.1
    Umeansall$Group.1 <- NULL

    originalCells[[n]] <- Umeansall

  }

  originalCells <- do.call("rbind", originalCells)

  originalCells <- originalCells[order(originalCells$max.length),]
  cellnumdatframe<- data.frame(cellnum.length = 1:length(unique(originalCells$max.length)), max.length = unique(originalCells$max.length))
  originalCells <- merge(originalCells, cellnumdatframe)

  originalCells <- originalCells[order(originalCells$max.width),]
  cellnumdatframe<- data.frame(cellnum.width = 1:length(unique(originalCells$max.width)), max.width=unique(originalCells$max.width))
  originalCells <- merge(originalCells, cellnumdatframe)

  return(originalCells)
}


##############plot#########
#' @export

bactKymo <- function(originalCells, timeD = FALSE, dimension = "length", bins=25, sizeAV=FALSE, cells="all", prep=TRUE, percDiv=FALSE, cutoff_demograph = 0.975, mag, legend=TRUE){

  measure <- "pixels"
  if(legend==TRUE){
    pos <- "right"
  }
  if(legend==FALSE){
    pos <- "none"
  }
  if(percDiv==TRUE){
    if("percentage_binned"%in%colnames(originalCells)!=TRUE){
      groupP <- perc_Division(unique(originalCells[,c("cell", "frame", "max.length")]), av=FALSE, plotgrowth=FALSE)$timelapse
      groupP <- unique(groupP[,c("cell","frame","percentage_binned")])
    }
    if("percentage_binned"%in%colnames(originalCells)==TRUE){
      groupP <- unique(originalCells[,c("cell", "frame", "percentage_binned")])
    }

  }

  if(prep==TRUE){
    originalCells <- prepForKymo(originalCells, dimension=dimension, bins=bins, sizeAV=sizeAV)
    if(percDiv==TRUE){
      originalCells <- merge(originalCells,groupP)
    }
  }

  if(percDiv==TRUE){
   originalCells$frame <- as.numeric(originalCells$percentage_binned)
   originalCells$cell <- 1

   if(sizeAV==FALSE){
     originalCells <- stats::aggregate(originalCells[,colnames(originalCells)[colnames(originalCells)!="group"&colnames(originalCells)!="frame"&colnames(originalCells)!="percentage_binned"]],
                         by=list(group=originalCells$group, frame=originalCells$frame),FUN=mean)
   }
   if(sizeAV==TRUE){
     originalCells <- stats::aggregate(originalCells[,colnames(originalCells)[colnames(originalCells)!="group"&colnames(originalCells)!="frameh"&colnames(originalCells)!="frame"&colnames(originalCells)!="percentage_binned"]],
                      by=list(group=originalCells$group, frameh=originalCells$frameh, frame=originalCells$frame), FUN=mean)

     originalCells$x_coords[originalCells$frameh=="a"|originalCells$frameh=="d"] <- originalCells$frame[originalCells$frameh=="a"|originalCells$frameh=="d"] + 0.5
     originalCells$x_coords[originalCells$frameh=="b"|originalCells$frameh=="c"] <- originalCells$frame[originalCells$frameh=="b"|originalCells$frameh=="c"] - 0.5
   }
  }

  if(timeD==FALSE&percDiv==F){
    #97.5 % cutoff for outlying superbright stuff

    if(cells!="all"){
      if(length(cells)>1){
        if(is.numeric(cells)!=TRUE){
          stop("'cells' is neither numeric nor 'all', please set 'cells' to a numeric vector or 'all'")
        }
        if(is.numeric(cells)==TRUE){
          originalCells <- originalCells[originalCells$cell%in%cells,]
        }
      }
      if(length(cells)<=1){
        stop("'cells' is of length 1, while there is no time dimension.
              \n if you want to plot a single cell over time, set 'timeD' to 'TRUE'
              \n if you want to plot all cells as a demograph, put 'cells' to 'all'.
              \n if you want to plot a specific group of cells as a demograph, put 'cells' to a vector identifying the cell numbers")
      }
    }
    originalCells <- originalCells[originalCells$values<stats::quantile(originalCells$values, cutoff_demograph),]

    if(dimension=="length"&sizeAV==FALSE){
      plot1 <- ggplot2::ggplot(originalCells, ggplot2::aes_string(x='cellnum.length',y='group', fill='values')) +
        ggplot2::geom_raster() +
        ggplot2::coord_fixed(ratio=20) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("n(th) cell ordered by cell length") +
        ggplot2::ylab("bin (by cell length)") +
        ggplot2::scale_fill_viridis_c(name="Fluorescence\nIntensity") +
        ggplot2::theme(legend.position=pos)
    }

    if(dimension=="length"&sizeAV==TRUE){
      originalCells$cellnum.length[originalCells$frameh=="a"|originalCells$frameh=="d"] <- originalCells$cellnum.length[originalCells$frameh=="a"|originalCells$frameh=="d"] + 0.5
      originalCells$cellnum.length[originalCells$frameh=="b"|originalCells$frameh=="c"] <- originalCells$cellnum.length[originalCells$frameh=="b"|originalCells$frameh=="c"] - 0.5

      if(!missing(mag)){
        originalCells$y_coords <- originalCells$y_coords * unlist(get(magnificationList, envir=magEnv)[mag])
        measure <- "micron"
      }
      originalCells$grouping <- paste(originalCells$X_rot, originalCells$cell, originalCells$frame, originalCells$pointN, sep="_")
      plot1 <- ggplot2::ggplot(originalCells, ggplot2::aes_string(x='cellnum.length',y='y_coords', group= 'grouping', fill='values')) +
        ggplot2::geom_polygon() +
        ggplot2::theme_minimal() +
        ggplot2::xlab("n(th) cell ordered by cell length") +
        ggplot2::ylab(paste("location on length axis (", measure, ")", sep="")) +
        ggplot2::scale_fill_viridis_c(name="Fluorescence\nIntensity") +
        ggplot2::theme(legend.position=pos)
    }
    if(dimension=="width"&sizeAV==FALSE){
      plot1 <- ggplot2::ggplot(originalCells, ggplot2::aes_string(x='cellnum.width',y='group', fill='values')) +
        ggplot2::geom_raster() +
        ggplot2::coord_fixed(ratio=20) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("n(th) cell ordered by cell width") +
        ggplot2::ylab("bin (by cell width)") +
        ggplot2::scale_fill_viridis_c(name="Fluorescence\nIntensity") +
        ggplot2::theme(legend.position=pos)
    }
    if(dimension=="width"&sizeAV==TRUE){
      originalCells$cellnum.width[originalCells$frameh=="a"|originalCells$frameh=="d"] <- originalCells$cellnum.width[originalCells$frameh=="a"|originalCells$frameh=="d"] + 0.5
      originalCells$cellnum.width[originalCells$frameh=="b"|originalCells$frameh=="c"] <- originalCells$cellnum.width[originalCells$frameh=="b"|originalCells$frameh=="c"] - 0.5

      if(!missing(mag)){
        originalCells$y_coords <- originalCells$y_coords * unlist(get(magnificationList, envir=magEnv)[mag])
        measure <- "micron"
      }
      originalCells$grouping <- paste(originalCells$X_rot, originalCells$cell, originalCells$frame, originalCells$pointN, sep="_")
      plot1 <- ggplot2::ggplot(originalCells, ggplot2::aes_string(x='cellnum.width',y='y_coords', group='grouping', fill='values')) +
        ggplot2::geom_polygon() +
        ggplot2::theme_minimal() +
        ggplot2::xlab("n(th) cell ordered by cell width") +
        ggplot2::ylab(paste("location on length axis (", measure, ")", sep="")) +
        ggplot2::scale_fill_viridis_c(name="Fluorescence\nIntensity") +
        ggplot2::theme(legend.position=pos)
    }
  }

  if(timeD==TRUE&sizeAV==FALSE&percDiv==FALSE){
    if(cells=="all"){
      plot1 <- lapply(unique(originalCells$cell), function(x) ggplot2::ggplot(originalCells[originalCells$cell==x,]) +
                        ggplot2::geom_raster(ggplot2::aes_string(x='frame', y='group', fill='values')) +
                        ggplot2::theme_minimal() +
                        ggplot2::xlab("Time (frames)") +
                        ggplot2::ylab(paste("bin (by cell ", dimension, ")", sep="")) +
                        ggplot2::scale_fill_viridis_c() +
                        ggplot2::ggtitle(paste("Cell", x, sep=" ")) +
                        ggplot2::theme(legend.position=pos)
      )
      names(plot1) <- paste("cell", unique(originalCells$cell), sep="")
    }
    if(is.numeric(cells)==TRUE&length(cells)>1){
      plot1 <- lapply(cells, function(x) ggplot2::ggplot(originalCells[originalCells$cell==x,]) +
                        ggplot2::geom_raster(ggplot2::aes_string(x='frame', y='group', fill='values')) +
                        ggplot2::theme_minimal() +
                        ggplot2::xlab("Time (frames)") +
                        ggplot2::ylab(paste("bin (by cell ", dimension, ")", sep="")) +
                        ggplot2::scale_fill_viridis_c() +
                        ggplot2::ggtitle(paste("Cell", x, sep=" ")) +
                        ggplot2::theme(legend.position=pos)
      )
      names(plot1) <- paste("cell", cells, sep="")
    }
    if(is.numeric(cells)==TRUE&length(cells)==1){
      plot1 <- ggplot2::ggplot(originalCells[originalCells$cell==cells,]) +
        ggplot2::geom_raster(ggplot2::aes_string(x='frame', y='group', fill='values')) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("Time (frames)") +
        ggplot2::ylab(paste("bin (by cell ", dimension, ")", sep="")) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::ggtitle(paste("Cell", cells, sep=" ")) +
        ggplot2::theme(legend.position=pos)
    }
    if(cells!="all" & is.numeric(cells)==FALSE){
      stop("'cells' must be either be a character string 'all', a number or a vector of numbers c(x, x,)")
    }
  }

  if(timeD==TRUE&sizeAV==TRUE&percDiv==FALSE){
    if(cells=="all"){

      if(!missing(mag)){
        originalCells$y_coords <- originalCells$y_coords * unlist(get(magnificationList, envir=magEnv)[mag])
        #originalCells$x_coords <- originalCells$x_coords * unlist(get(magnificationList, envir=magEnv)[mag])
        measure <- "micron"
      }
        originalCells$grouping <- paste(originalCells$X_rot, originalCells$pointN, originalCells$values, sep="_")
        plot1 <- lapply(unique(originalCells$cell[order(originalCells$cell)]), function(x) ggplot2::ggplot(originalCells[originalCells$cell==x,]) +
                          ggplot2::geom_polygon(ggplot2::aes_string(x='x_coords', y='y_coords', fill='values', group='grouping')) +
                          ggplot2::theme_minimal() +
                          ggplot2::xlab("Time (frames)") +
                          ggplot2::ylab(paste(dimension, " (in ", measure, ")", sep="")) +
                          ggplot2::scale_fill_viridis_c() +
                          ggplot2::ggtitle(paste("Cell", x, sep=" ")) +
                          ggplot2::theme(legend.position=pos)
        )
        names(plot1) <- paste("cell", unique(originalCells$cell[order(originalCells$cell)]), sep="")
        plot1 <- plot1[!is.na(plot1)]

    }
    if(is.numeric(cells)==TRUE&length(cells)>1){
      if(!missing(mag)){
        originalCells$y_coords <- originalCells$y_coords * unlist(get(magnificationList, envir=magEnv)[mag])
        originalCells$x_coords <- originalCells$x_coords * unlist(get(magnificationList, envir=magEnv)[mag])
        measure <- "micron"
      }
      originalCells$grouping <- paste(originalCells$X_rot, originalCells$pointN, originalCells$values, sep="_")
      plot1 <- lapply(cells, function(x) ggplot2::ggplot(originalCells[originalCells$cell==x,]) +
                        ggplot2::geom_polygon(ggplot2::aes_string(x='x_coords', y='y_coords', fill='values', group='grouping')) +
                        ggplot2::theme_minimal() +
                        ggplot2::xlab("Time (frames)") +
                        ggplot2::ylab(paste(dimension, " (in ", measure, ")", sep="")) +
                        ggplot2::scale_fill_viridis_c() +
                        ggplot2::ggtitle(paste("Cell", x, sep=" ")) +
                        ggplot2::theme(legend.position=pos)
      )
      names(plot1) <- paste("cell", cells, sep="")
    }
    if(is.numeric(cells)==TRUE&length(cells)==1){
      if(!missing(mag)){
        originalCells$y_coords <- originalCells$y_coords * unlist(get(magnificationList, envir=magEnv)[mag])
        originalCells$x_coords <- originalCells$x_coords * unlist(get(magnificationList, envir=magEnv)[mag])
        measure <- "micron"
      }
      originalCells$grouping <- paste(originalCells$X_rot, originalCells$pointN, originalCells$values, sep="_")
      plot1 <- ggplot2::ggplot(originalCells[originalCells$cell==cells,]) +
        ggplot2::geom_polygon(ggplot2::aes_string(x='x_coords', y='y_coords', fill='values', group='grouping')) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("Time (frames)") +
        ggplot2::ylab(paste(dimension, " (in ", measure, ")", sep="")) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::ggtitle(paste("Cell", cells, sep=" ")) +
        ggplot2::theme(legend.position=pos)
    }
    if(cells!="all" & is.numeric(cells)==FALSE){
      stop("'cells' must be either be a character string 'all', a number or a vector of numbers c(x, x,)")
    }
  }

  if(percDiv==TRUE&sizeAV==FALSE){
    plot1 <- ggplot2::ggplot(originalCells) +
      ggplot2::geom_raster(ggplot2::aes_string(x='frame', y='group', fill='values')) +
      ggplot2::theme_minimal() +
      ggplot2::xlab("Percentage of division") +
      ggplot2::ylab(paste("length (by cell ", dimension, ")", sep="")) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::theme(legend.position=pos)
  }


  if(percDiv==TRUE&sizeAV==TRUE){
    if(!missing(mag)){
      originalCells$y_coords <- originalCells$y_coords * unlist(get(magnificationList, envir=magEnv)[mag])
      #originalCells$x_coords <- originalCells$x_coords * unlist(get(magnificationList, envir=magEnv)[mag])
      measure <- "micron"
    }
    originalCells$grouping <- paste(originalCells$X_rot, originalCells$values, sep="_")
    originalCells$x_coords <- originalCells$x_coords*10

    plot1 <- ggplot2::ggplot(originalCells) +
      ggplot2::geom_polygon(ggplot2::aes_string(x='x_coords', y='y_coords', fill='values', group='grouping')) +
      ggplot2::theme_minimal() +
      ggplot2::xlab("Percentage of division") +
      ggplot2::ylab(paste("length by cell ", dimension, "( in", measure, ")", sep="")) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::theme(legend.position=pos)
  }

  return(plot1)

}

