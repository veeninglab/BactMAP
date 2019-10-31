##9-1-2015- last update 13-6-2016
##Renske van Raaphorst

#preparation of dataset containing spots for plotting.

#might be good to clean up first:

#before running the script, open datasets:
#"REP": output of peakfitter.
#"M": the excel output file of the meshes used for peakfitter
#aand start



#!! --> if using .txt extension: make sure the "length" cells are notated as numeric, with the same amount of decimal points
#       for both data frames before importing.
#!! --> also check your decimal seperator for the tab delimited .txt files.

###############################################################################################
#' @export
addPalette <- function(palList, palName){
  if(missing(palName)){
    palName <- readline("Give the name of your new Palette:\n")
  }
  if(missing(palList)){
    low <- readline("Give the color (hex) for the lowest value:\n")
    med <- readline("Give the color (hex) for the midpoint value:\n")
    hig <- readline("Give the color (hex) for the highest value:\n")
    palList <- list(low, med, hig)
  }
  palList <- list(palList)
  names(palList) <- palName
  colopts2 <- append(get(colopts, envir=colEnv), palList)
  assign(colopts, colopts2, envir=colEnv)
  showCurrentPalettes()

}

#' @export
showCurrentPalettes <- function(){
  colchoice = get(colopts, envir=colEnv)
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    if (!requireNamespace("mvtnorm", quietly = TRUE)) {
      inp <- readline("Package 'gridExtra' and 'mvtnorm' needed for this function to work. Press 'y' to install them, or any other key to cancel.")
      if(inp=="y"|inp=="Y"){utils::install.packages(c("gridExtra", "mvtnorm"))}else{stop("Canceled")}
    }
    else{inp <- readline("Package 'gridExtra' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){utils::install.packages("gridExtra")}else{stop("Canceled")}
    }
  }
  if(!requireNamespace("mvtnorm", quietly = TRUE)) {
    inp <- readline("Package 'mvtnorm' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){utils::install.packages("mvtnorm")}else{stop("Canceled")}
  }
  plotlist <- list()
  namelist <- names(colchoice)
  for(n in 1:length(colchoice)){
    p <- suppressWarnings(colorchoiceplot(colchoice[[n]], n, namelist[n]))
    plotlist <- append(plotlist,list(p))
  }
  return(suppressWarnings(do.call(gridExtra::grid.arrange,c(plotlist, ncol=2))))
}

#' @export
getPalette <- function(palName){
  if(missing(palName)){
    return(get(colopts, envir=colEnv))
  }else{
    collist <- get(colopts, envir=colEnv)
    return(collist[palName][[1]])
  }
}

#' @importFrom ggplot2 stat
#' @importFrom stats density
colorchoiceplot <- function(colchoice, nums, pname){

  z <- mvtnorm::rmvnorm(100, mean=c(3,5), sigma=matrix(c(1,0.5,0.5,2), nrow=2))
  z <- data.frame(z)
  return(suppressWarnings(ggplot2::ggplot(z, ggplot2::aes_string(x='X1', y='X2')) +
                            suppressWarnings(ggplot2::stat_density2d(ggplot2::aes(fill=stat(density)), geom="raster", contour=FALSE)) +
                            ggplot2::scale_fill_gradient2(low = colchoice[1], mid= colchoice[2], high = colchoice[3], midpoint=0.06) +
                            ggplot2::theme_minimal() + ggplot2::theme(legend.position="none") +
                            ggplot2::ggtitle(pname) + ggplot2::xlab("") +
                            ggplot2::ylab("") +
                            ggplot2::xlim(0,5) +
                            ggplot2::ylim(0,10) +
                            ggplot2::theme(axis.text=ggplot2::element_blank())))

}

#' @export
pal2Default <- function(){
  assign(colopts, cols, envir=colEnv)
  message("Colorpalettes back to default")
}

#wat basisplotfuncties
densityplot <- function(plot){
  return(
    plot + ggplot2::stat_density2d(ggplot2::aes(fill=stat(density)), geom="raster", contour = FALSE))

    }

LWplot <- function(plot, u="black", maxn){
  return(plot + #ggplot2::geom_rect(xmin=0, xmax=maxn, ymin=-1.5, ymax=1.5, fill=u)
           ggplot2::theme_minimal()  + ggplot2::scale_x_continuous(limits=c(0,maxn)))
}

heatmap <- function(pdens, mp, colchoice, viridis = F){
  if(viridis==F){
    return(pdens + ggplot2::scale_fill_gradient2(low = colchoice[1], mid= colchoice[2], high = colchoice[3], midpoint =mp, space = "Lab") + ggplot2::theme_minimal() + ggplot2::theme(legend.position="none", panel.background=ggplot2::element_rect(fill=colchoice[1])))
  }
  if(viridis==T){
    return(pdens + ggplot2::scale_fill_viridis_c(option=colchoice) + ggplot2::theme_minimal() + ggplot2::theme(legend.position="none", panel.background=ggplot2::element_rect(fill="black")))
  }
}

#function for goodlooking x/y coordinate plot:
#makes a plot sized as the max cell (width/length: xmax) of the quartile inside a plot
#which has a grey background as large as the largest cell in the dataset
#so all quartile plots will have the same dimensions.
#the title, y axis and x axis will also be drawn.
coplot <- function(pheat, xmax, ymax, u="black"){
  return(suppressMessages(pheat + ggplot2::xlab("Length (\u00B5m)") + ggplot2::ylab("Width (\u00B5m)") + ggplot2::coord_fixed(xlim = c(-xmax,xmax), ylim=c(-ymax,ymax)) + ggplot2::theme(panel.background = ggplot2::element_rect(fill = u), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())))
}

#plotlocation_histograms <- function(){}



######################### mid-points ################################################################################

superfun <- function(dat, bins,mag){
  dat <- dat[order(dat$frame, dat$cell, dat$num),]
  dat <- unique(dat)
  dat$av <- 0
  if("max.length"%in%colnames(dat)!=T){
    dat$max.length <- dat$length
  }
  dat$av <- dat$X_rot/dat$max.length*100
  dat <- dat[!is.na(dat$av),]

  cutpoints<-stats::quantile(dat$av,(0:bins)/bins)

  dat$binned1 <-cut(dat$av,cutpoints, include.lowest=TRUE, labels = 1:bins)

  xb <- split(dat$X_rot, dat$binned)
  ybp <- split(dat$Y_rot[dat$Y_rot>0], dat$binned[dat$Y_rot>0])
  ybm <- split(dat$Y_rot[dat$Y_rot<0], dat$binned[dat$Y_rot<0])
  xmeans <- lapply(xb, function(x)mean(x))
  ypmeans <- lapply(ybp, function(x)mean(x))
  ymmeans <- lapply(ybm, function(x)mean(x))
  Lb <- split(dat$max.length, dat$binned1)
  maxLmeans <- lapply(Lb, function(x)mean(x))
  meanframe <- data.frame(x=unlist(xmeans), yp=unlist(ypmeans), ym=unlist(ymmeans),max.length=unlist(maxLmeans))
  meanframe <- reshape2::melt(meanframe, c("x","max.length"), value.name="y")
  meanframe$variable <- NULL
  meanframe <- meanframe * mag
  meanframe <- meanframe[abs(meanframe$x)<=0.5*meanframe$max.length,]
  #add num as an order identifyer for drawing a path or a polygon as cell
  meanframe$num <- c(1:(nrow(meanframe)/2), nrow(meanframe):(nrow(meanframe)/2+1))
  meanframe <- rbind(meanframe, meanframe[1,])
  meanframe[nrow(meanframe),]$num <- nrow(meanframe)
  meanframe <- meanframe[order(meanframe$num),]
  return(meanframe)
}

#or two, by quartiles of the number of cells:
#' @export
createPlotList <- function(spotdata,  meshdata, groups =4 , colorpalette="GreenYellow", mag="No_PixelCorrection", AllPlot=F, Xm="X", Ym="Y", viridis=FALSE, showPlot=TRUE, getData=FALSE, getSummary=TRUE){
  if (!requireNamespace("MASS", quietly = TRUE)) {
    inp_P <- readline("Package 'MASS' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp_P=="y"|inp_P=="Y"){utils::install.packages("MASS")}else{stop("Canceled")}
  }
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    inp_P <- readline("Package 'gridExtra' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp_P=="y"|inp_P=="Y"){utils::install.packages("gridExtra")}else{stop("Canceled")}
  }
  if(AllPlot==T){
    if (!requireNamespace("grid", quietly = TRUE)) {
      inp_P <- readline("Package 'grid' needed for this function to work. Press 'y' to install, or any other key to cancel.")
      if(inp_P=="y"|inp_P=="Y"){utils::install.packages("grid")}else{stop("Canceled")}
    }
  }
  if(missing(mag)!=T&is.numeric(unlist(get(magnificationList,envir=magEnv)[mag]))==FALSE){
    stop("Magnification conversion factor not recognized. Please use addPixels2um('pixelName', pixelsize) to add your conversion factor")
  }
  if("Lmid"%in%colnames(spotdata)==T){
    MR <- spotdata
    MR$cellnum <- MR$num
  }
  if("l"%in%colnames(spotdata)==T&"Lmid"%in%colnames(spotdata)==F){
    MR <- LimDum(spotdata, pix2um = unlist(get(magnificationList, envir=magEnv)[mag]))
    MR <- spotMR(MR)
    #if(length(which(unique(meshdata$cell)%in%spotdata$cell!=T))==0){MR <- spotdata}else{
    MR <- merge(MR, unique(meshdata[,c("cell",
                                       "frame",
                                       "max.width",
                                       "max.length")]),
                all=TRUE)
    #}
    MR$totalspot[is.na(MR$totalspot)] <- 0
    MR <- MR[order(MR$max.length, MR$max.width),]
    MR <- MR[!is.na(MR$cell),]
    MR$num <- c(1:nrow(MR))
    MR$cellnum <- MR$num
  }
  if("Lmid"%in%colnames(spotdata)==F&"l"%in%colnames(spotdata)==F){
    message("Did not find cell data in spot dataset. Running spotsInBox to connect spots to meshes...")
    combineframes <- spotsInBox(spotdata, meshdata, Xm=Xm,Ym=Ym)
    meshdata <- combineframes$mesh
    spotdata <- combineframes$spots_relative

    MR <- LimDum(spotdata, pix2um = unlist(get(magnificationList, envir=magEnv)[mag]))
    MR <- spotMR(MR)
    #if(length(which(unique(meshdata$cell)%in%spotdata$cell!=T))==0){MR <- spotdata}else{
    MR <- merge(MR, unique(meshdata[,c("cell",
                                       "frame",
                                       "max.width",
                                       "max.length")]),
                all=TRUE)
    #}

    MR$totalspot[is.na(MR$totalspot)] <- 0
    MR <- MR[order(MR$max.length, MR$max.width),]
    MR <- MR[!is.na(MR$cell),]
    MR$num <- c(1:nrow(MR))
    MR$cellnum <- MR$num
  }
  if("obID"%in%colnames(spotdata)){
    if("max_um"%in%colnames(spotdata)!=T){
      MR <- unique(MR[,c("num", "frame", "cell", "max.length", "max.width", "Dum", "Lmid", "pole1", "pole2", "cellnum", "obnum", "obID")])
    }
    if("max_um"%in%colnames(spotdata)){
      MR <- unique(MR[,c("num", "frame", "cell", "max.length", "max.width", "Dum", "Lmid", "pole1", "pole2", "cellnum", "max_um", "maxwum", "obnum", "obID")])
    }
    MR <- MR[!is.na(MR$cell),]
  }

  if(viridis==FALSE){
    colc <- get(colopts, envir = colEnv)[colorpalette][[1]]
  }

  if(viridis==TRUE){
    colc <- colorpalette
  }

  if(groups>1){
    MR$q1 <- cut(MR$cellnum, breaks=groups, labels = 1:groups)
  }
  if(groups==1){
    MR$q1 <- 1
  }

  xmax <- 0.5*max(MR$max.length, na.rm=TRUE)
  ymax <- 0.5*max(MR$max.width, na.rm=TRUE)

  #Length and width corrected for average length per quartile for plotting the coordinate plots
  mframe <- suppressWarnings(stats::aggregate(MR, by=list(MR$q1), FUN=mean, na.rm=TRUE, na.action=stats::na.omit))
  meansL <- mframe$max.length
  meansW <- mframe$max.width

  #length correction
  MR$Lcor <- (MR$Lmid/MR$max.length*meansL[as.numeric(MR$q1)])
  #width correction
  MR$Dcor <- (MR$Dum/MR$max.width*meansW[as.numeric(MR$q1)])

  ###############################################################################################################################################
  ##plotting! -> coordinate plots
  allMRs <- split(MR, MR$q1)
  plist <- list()
  for(n in 1:groups){
    p <- ggplot2::ggplot(allMRs[[n]], ggplot2::aes_string(x='Lcor', y='Dcor'))
    p <- densityplot(p)
    plist <- append(plist, list(p))
  }

  #plotting! -> L and D ordered by cell length
  pL <- ggplot2::ggplot(MR, ggplot2::aes_string(x='num', y='Lmid'))
  if(viridis==TRUE){
    pL <- LWplot(pL, "black", max(MR$num, na.rm=T))
  }
  if(viridis==FALSE){
    pL <- LWplot(pL, colc[[1]], max(MR$num, na.rm=T))
  }
  #pLpoint <- pL + ggplot2::geom_point() + ggplot2::ggtitle("Spot location on length axis ordered by cell length") + ggplot2::xlab("Spot - ordered by cell length") + ggplot2::ylab("Y-position (\u03BCm)") + ggplot2::theme_bw()
  pLD <- densityplot(pL) +
    ggplot2::ggtitle("Spot location on length axis ordered by cell length") +
    ggplot2::xlab("Spot - ordered by cell length") +
    ggplot2::ylab("Y-position (\u03BCm)") +
    ggplot2::geom_line(data=MR, ggplot2::aes_string(x='num',y='pole1'),colour="white") +
    ggplot2::geom_line(data=MR, ggplot2::aes_string(x='num',y='pole2'),colour="white")

  pW <- ggplot2::ggplot(MR, ggplot2::aes_string(x='num', y='Dum'))
  if(viridis == TRUE){
    pW <- LWplot(pW, "black", max(MR$num,na.rm=T))
  }
  if(viridis==FALSE){
    pW <- LWplot(pW, colc[[1]], max(MR$num,na.rm=T))
  }
  #pWpoint <- pW + ggplot2::geom_point() + ggplot2::ggtitle("Spot location on width axis ordered by cell length") + ggplot2::xlab("Spot - ordered by cell length") + ggplot2::ylab("X-position (\u03BCm)") + ggplot2::theme_bw()
  pWD <- densityplot(pW) +
    ggplot2::ggtitle("Spot location on width axis ordered by cell length") +
    ggplot2::xlab("Spot - ordered by cell length") +
    ggplot2::ylab("X-position (\u03BCm)")

  #make heatmap using the half max densities:

  ##for all quartiles (or n-tiles)
  allMRLmids <- split(MR[,c("Lmid","Dum")][!is.na(MR$Lmid),], MR$q1[!is.na(MR$Lmid)])
  #get maximum half max density of all groups and apply for each plot to have same scaling.
  mp<- max(sapply(allMRLmids, function(x) stats::median(range(MASS::kde2d(x$Lmid,x$Dum)$z))))
  #create n-tile plots saved in a list.
  phlist <- list()
  for(n in 1:groups){
    p <- heatmap(plist[[n]], mp, colc, viridis)
    p <- coplot(p, xmax, ymax)
    phlist <- append(phlist, list(p))
  }

  u <- list()
  #do the same procedure for length/width heatmap plots
  mpL <- MASS::kde2d(MR$num[!is.na(MR$Lmid)], MR$Lmid[!is.na(MR$Lmid)])
  mpL1 <- stats::median(range(mpL$z))
  pLD <- heatmap(pLD, mpL1, colc,viridis)
  if(viridis==FALSE){
    pLD <- pLD + ggplot2::theme_minimal() +
      ggplot2::theme(panel.background=ggplot2::element_rect(fill=colc[[1]]), panel.grid = ggplot2::element_blank())
  }
  if(viridis==TRUE){
    pLD <- pLD + ggplot2::theme_minimal() +
      ggplot2::theme(panel.background=ggplot2::element_rect(fill="black"), panel.grid = ggplot2::element_blank())
  }
  u$lengthplot <- pLD

  mpW <- MASS::kde2d(MR$num[!is.na(MR$Dum)], MR$Dum[!is.na(MR$Dum)])
  mpW1 <- stats::median(range(mpW$z))
  pWD <- heatmap(pWD, mpW1, colc,viridis)
  if(viridis==TRUE){
    pWD <- pWD + ggplot2::theme_minimal() +
      ggplot2::theme(panel.background=ggplot2::element_rect(fill="black"), panel.grid = ggplot2::element_blank())
  }

  if(viridis==FALSE){
    pWD <- pWD + ggplot2::theme_minimal() +
      ggplot2::theme(panel.background=ggplot2::element_rect(fill=colc[[1]]), panel.grid = ggplot2::element_blank())
  }
  u$widthplot <- pWD

  if(missing(meshdata)){
    u$qplots <- phlist
  }

  #################add meshdata data###############################################################################################
  if(!missing(meshdata)){

    if("X_rot"%in%colnames(meshdata)!=T){
      message("Turning the cells... this may take a while.")

      meshdata <- meshTurn(meshdata, Xm, Ym)
      message("Finished turning the cells")
    }
    p2um <- as.numeric(get(magnificationList, envir=magEnv)[mag])
    if("max_um"%in%colnames(meshdata)!=T){
      meshdata$max_um <- meshdata$max.length*p2um
      meshdata$maxwum <- meshdata$max.width*p2um
    }


    meshdata <- merge(meshdata, MR[,c("cell", "frame", "q1")], all=T)
    MESHlist <- split(meshdata, meshdata$q1)
    message("Calculating mean cell outlines..")
    if(nrow(unique(meshdata[meshdata$max.length==max(meshdata$max.length),]))< 12){
      means <- lapply(MESHlist, function(x) superfun(x, 12, p2um))
    }
    if(nrow(unique(meshdata[meshdata$max.length==max(meshdata$max.length),]))>11){
      means <- lapply(MESHlist, function(x)superfun(x, round(min(x$max.length), digits=0), p2um))
    }
    u$mean_outlines <- means
    message("Finished calculating mean cell outlines")

    phmlist <- list()
    for(n in 1:groups){
      p <- phlist[[n]] +
        ggplot2::geom_path(data=means[[n]], ggplot2::aes_string(x='x',y='y'), colour="white") +
        ggplot2::coord_fixed(xlim=c(min(means[[groups]]$x)*1.2, max(means[[groups]]$x*1.2)), ylim=c(min(means[[groups]]$y)*1.2, max(means[[groups]]$y*1.2)))
      phmlist <- append(phmlist, list(p))
    }
    if(AllPlot==F){
      u$qplots <- gridExtra::arrangeGrob(grobs=phmlist, ncol=1)
      u$qplots_separate <- phmlist
    }
  }
  #and create the plots:
  if(AllPlot==T){
    phmalist <- list()
    for(n in 1:groups){
      if(missing(meshdata)!=TRUE){
        p1_all <- suppressMessages(allplot(phmlist[[n]], allMRs[[n]], max(means[[groups]]$x*1.2), max(means[[groups]]$y*1.2), empty))
      }
      if(missing(meshdata)==TRUE){
        p1_all <- suppressMessages(allplot(phlist[[n]], allMRs[[n]], mean(allMRs[[groups]]$max.length)*1.2, mean(allMRs[[groups]]$max.width)*1.2, empty))
      }
      phmalist <- append(phmalist,list(p1_all))
      u$qplots <- gridExtra::arrangeGrob(grobs=phmalist, ncol=1)
      u$qplots_separate <- phmalist
    }
  }

  #in case you want it for the whole thing instead of only quartiles:

  pall <- ggplot2::ggplot(MR, ggplot2::aes_string(x='Lcor', y='Dcor'))
  pall <- densityplot(pall)
  mppall <- mean(range(MASS::kde2d(MR$Lcor[!is.na(MR$Lcor)&!is.na(MR$Dcor)],MR$Dcor[!is.na(MR$Dcor)&!is.na(MR$Lcor)])$z))
  pall <- heatmap(pall, mppall, colc, viridis)


  if(nrow(unique(meshdata[meshdata$max.length==max(meshdata$max.length),]))< 12){
    meantotal <- superfun(meshdata, 12, p2um)
  }
  if(nrow(unique(meshdata[meshdata$max.length==max(meshdata$max.length),]))> 11){
    meantotal <- superfun(meshdata, 30, p2um)
  }
  pall <- coplot(pall,max(meantotal$x)*1.2, max(meantotal$y)*1.2)
  pall <- pall + ggplot2::geom_path(data=meantotal, ggplot2::aes_string(x='x',y='y'), colour="white")
  if(AllPlot==F){
    u$plottotal <- pall
  }
  if(AllPlot==T){
    pall_all <- suppressMessages(allplot(pall, MR, max(meantotal$x)*1.2, max(meantotal$y)*1.2, empty))
    u$plottotal <- pall_all

  }




  hislist <- list()
  #save all histograms (L coordinates) of the quartiles too:
  for(n in 1:groups){
    if(viridis==FALSE){
      p1his <- ggplot2::ggplot(allMRs[[n]], ggplot2::aes_string(x='Lcor')) +
        ggplot2::geom_density(fill=colc[[2]], color=colc[[2]]) +
        ggplot2::theme_minimal() +
        ggplot2::labs(x="Length(\u03BCm)") +
        ggplot2::coord_cartesian(xlim=c(min(MR$Lcor, na.rm=T),max(MR$Lcor, na.rm=T)))
    }
    if(viridis==TRUE){
      singopt <- viridissinglecols[colc]
      p1his <- ggplot2::ggplot(allMRs[[n]], ggplot2::aes_string(x='Lcor')) +
        ggplot2::geom_density(fill=singopt, color=singopt) +
        ggplot2::theme_minimal() +
        ggplot2::labs(x="Length(\u03BCm)") +
        ggplot2::coord_cartesian(xlim=c(min(MR$Lcor,na.rm=T), max(MR$Lcor,na.rm=T)))

    }
    hislist <- append(hislist, list(p1his))
  }

  u$histograms <- gridExtra::arrangeGrob(grobs=hislist, ncol=1)
  if(getData==TRUE){
    u$spotdata <- MR
    if(!missing(meshdata)){
      u$meshdata <- meshdata
    }
  }
  u$pixel2um <- p2um
  if(getSummary==TRUE){
    u$data_summary <- makePlotListSummary_2(MR, groups=groups)
  }
  message("Done plotting.")
  if(showPlot==TRUE){
    n <- 1
    while(n%in%c(1:9)){
      n <- readline("Press the corresponding number to view:\n \n 1. lengthplot \n 2. widthplot \n 3. mean_outlines \n 4. qplots \n 5. histograms \n 6. spotdata (if getData = TRUE) \n 7. meshdata (if getData = TRUE) \n 8. pixel2um \n 9. data_summary (if getSummary=TRUE) \n 10. exit")
      if(n==1){
        message("lengthplot: spot localization on length axis plotted over cells ordered by cell length - density plot")
        graphics::plot(u$lengthplot)
      }
      if(n==2){
        message("widthplot: spot localization on width axis plotted over cells ordered by cell length - density plot")
        graphics::plot(u$widthplot)
      }
      if(n==3){
        message("data frame with the mean x/y coordinates of the cell outlines per size group. summary displayed below:")
        print(summary(u$mean_outlines))
      }
      if(n==4){
        message(paste("cell projections grouped in ", groups, " groups", sep=""))
        graphics::plot(u$qplots)
      }
      if(n==5){
        message(paste("spot localization on length axis in ", groups, " groups", sep=""))
        graphics::plot(u$histograms)
      }
      if(n==6){
        if(getData==TRUE){
          message("spot dataframe including grouping (column 'q1'). summary:")
          print(summary(u$spotdata))
        }
        if(getData==FALSE){
          message("spot dataframe not included in data output.")
        }
      }
      if(n==7){
        if(getData==TRUE){
          message("mesh dataframe including grouping (column 'q1'). summary:")
          print(summary(u$meshdata))
        }
        if(getData==FALSE){
          message("spot dataframe not included in data output.")
        }
      }
      if(n==8){
        message("pixel to micron conversion factor:")
        print(message(u$pixel2um))
      }
      if(n==9){
        if(getSummary==TRUE){
          message("summary of data")
          print(u$data_summary)
        }
        if(getSummary==FALSE){
          message("No data summary created (getSummary=FALSE).")
        }
      }
    }


  }

  return(u)
}

###########################################double histograms!! yay!!!######################################3333
##allplot function combines histograms and density plot.

allplot <- function(plot, data, xmax, ymax, empty){

  #prepare seperate plots: histograms (hL, hD) and modified coordinate plots(remove legend )
  p1D <- plot + ggplot2::theme(legend.position = "none")
  p1hL <- ggplot2::ggplot(data, ggplot2::aes_string(x='Lcor')) +
    ggplot2::geom_histogram(fill="black") +
    ggplot2::coord_cartesian(xlim = c(-xmax, xmax)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
  p1hD <- ggplot2::ggplot(data, ggplot2::aes_string(x='Dcor')) +
    ggplot2::geom_histogram(fill="black") +
    ggplot2::coord_flip(xlim = c(-ymax, ymax)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())

  #align the plots properly before putting them together
  p1Dg <- ggplot2::ggplotGrob(p1D)
  p1hLg <- ggplot2::ggplotGrob(p1hL)
  p1hDg <- ggplot2::ggplotGrob(p1hD)

  maxWidth = grid::unit.pmax(p1Dg$widths[2:5], p1hLg$widths[2:5])
  p1Dg$widths[2:5] <- as.list(maxWidth)
  p1hLg$widths[2:5] <- as.list(maxWidth)

  #put the grids together using gridarrange
  return(gridExtra::arrangeGrob(p1hLg, empty, p1Dg, p1hD, ncol=2, nrow=2, widths=c(10*xmax, 2.5*xmax), heights=c(2.5*xmax, 10*ymax)))
}

##before using the function:
#create mockup plot to make space
empty <- ggplot2::ggplot(data.frame(u=1), ggplot2::aes_string('u','u')) +
  ggplot2::theme(
    plot.background = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  )

###################################################



makePlotListSummary_1 <- function(spotdata, timelapse=FALSE){

  totalcells <- nrow(unique(spotdata[,c("cell", "frame")]))
  localization_x <- summary(spotdata$Lmid)
  localization_y <- summary(spotdata$Dum)

  max_length <- summary(unique(spotdata[,c("frame","cell", "max_um")])$max_um)
  max_width <- summary(unique(spotdata[,c("cell", "frame","maxwum")])$maxwum)

  out<- list("Total_amount_of_cells" = totalcells,
             "Spot_distance_from_midcell_length_axis" = localization_x,
             "Spot_distance_from_midcell_width_axis" = localization_y,
             "Cell_length" = max_length,
             "Cell_width" = max_width
  )
  if("totalspot"%in%colnames(spotdata)){
    spotscell <- summary(as.factor(unique(spotdata[,c("cell", "frame","totalspot")])$totalspot))
    out$Distribution_spots_per_cell <- spotscell
  }
  if("obID"%in%colnames(spotdata)){
    obcell <- summary(as.factor(stats::aggregate(spotdata$obnum, by=list(spotdata$cell, spotdata$frame), FUN=max)$x))
    out$Distribution_objects_per_cell <- obcell
  }
  return(out)

}

makePlotListSummary_2 <- function(spotdata, groups){
  if(groups==1){
    return(makePlotListSummary_1(spotdata))
  }
  if(groups>1){
    Full_dataset <- makePlotListSummary_1(spotdata)
    grouplist <- lapply(c(1:groups), function(x) makePlotListSummary_1(spotdata[spotdata$q1==x,]))
    names(grouplist) <- paste("Group", c(1:groups), sep="_")
    grouplist$Full_dataset <- Full_dataset
    return(grouplist)
  }
}

###################################################

# color palette lists & options. possible to add one yourself manually

cols <-  list(OrangeHot=list("#000000", "#D55E00", "#F0E442"), GreenYellow = list("#000000", "#009E73", "#F0E442"), ColdBlue = list("#000000", "#0072B2", "#FFFFFF"), YellowHot = list("#000000", "#e69f00", "#F0E442"), RedHot = list("#000000", "#FF0000", "#FFFF00"), WhiteOrange = list("#FFFFFF", "#F0E442", "#D55E00"))
colopts <- "colopts"
colEnv <- new.env()
assign(colopts, cols, envir=colEnv)


singcollist <- list("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#000000", "#0072B2", "#D55E00", "#CC79A7")

#'@export
colorsCUD <- function(colornames){
  if(!missing(colornames)){
    outlist <- singcollist[1:length(colornames)]
    names(outlist) <- colornames
    return(unlist(outlist))
  }
  if(missing(colornames)){
    return(unlist(singcollist))
  }
}

xaxislist <- c("Cell Length (\u03BCm)", "Cell Width (\u03BCm)", "Cell Area (\u03BCm\u00B2)",
               "Spot location in length axis (\u03BCm)", "Spot location on width axis (\u03BCm)")


viridissinglecols <- list(magma = "#B63679FF", A = "#B63679FF", inferno = "#BB3754FF", B ="#BB3754FF", plasma = "#CC4678FF", C = "#CC4678FF", cividis = "#7C7B78FF", E = "#7C7B78FF", viridis = "#21908CFF", D = "#21908CFF")

