####28.08.2018

###Renske van Raaphorst

##BactMAP


#Merge datasets, both manually after already analyzing the dataset & already when uploading.
#Plot a multicolor or multi-condition plot. Either next to each other or overlay depending on the type.


#1. combine.analyses: function(listofdataframes, listofconditions, listofcolors)
#'@export
combineDataframes <- function(listofdataframes, listofconditions, listofchannels ){

  #make sure there's no misalignment in the amount of conditions & channels.
  if(!missing(listofchannels)){
      if(length(listofdataframes)!=length(listofchannels)){
        stop("The amount of input dataframes does not match the amount of channels.")
      }
    }
  if(!missing(listofconditions)){
      if(length(listofdataframes)!=length(listofconditions)){
        stop("The amount of input dataframes does not match the amount of conditions.")
      }
  }

  #input should be able to be output of other functions (e.g. plot output). in that case:  take first-order data frames out of the output & proceed.
  if(is.list(listofdataframes[[1]])){
    #take only the dataframes out of the listofdataframes.
    u <- lapply(listofdataframes, function(bp) bp[which(lapply(bp, function(x) is.data.frame(x))==TRUE)])
    #now apply the addcolumn to

  }

  if(is.data.frame(listofdataframes[[1]])){
    if(!missing(listofconditions)){
      listofdataframes <- lapply(c(1:length(listofdataframes)), function(x) addColumn(listofdataframes[[x]], listofconditions[x], "condition"))
      if(!missing(listofchannels)){
        names(listofdataframes) <- paste(listofconditions,listofchannels,sep="_")
      }
      if(missing(listofchannels)){
        names(listofdataframes) <- listofconditions
      }
    }
    if(!missing(listofchannels)){
      listofdataframes <- lapply(c(1:length(listofdataframes)), function(x) addColumn(listofdataframes[[x]], listofchannels[x], "channel"))
      if(missing(listofconditions)){
        names(listofdataframes) <- listofchannels
      }
    }
  }


  listofdataframes_trimmed <- returnCommonColumn(listofdataframes)
  finalframe <- do.call('rbind', listofdataframes_trimmed)


  return(list("finalframe" = finalframe, "originaldata" = listofdataframes))

  }

addColumn <- function(datframe, value, name){
  datframe[,name] <- value
  return(datframe)
}

returnCommonColumn <- function(datframelist){
  colnameslist <- lapply(datframelist, function(x) colnames(x))
  for(n in 2:length(datframelist)){
   colnameslist[[1]] <- colnameslist[[1]][colnameslist[[1]]%in%colnameslist[[n]]]
  }
  datframelist <- lapply(datframelist, function(x) x[,colnameslist[[1]]])
  return(datframelist)
}
#2. plotoverlay(combinedanalyses, by=c("condition", "color" or "both)), type= length/width/projection, quantiles=n)
#'@export
plotOverlay <- function(meshdata,
                        spotdata,
                        objectdata,
                        by="both",
                        type="all",
                        quantiles = 1,
                        quantiles_by = "max.length",
                        equal_groups=TRUE,
                        mag,
                        objectcolor=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                        spotcolor = c("#0072B2", "#D55E00", "#CC79A7", "000000"),
                        histogram_outline = NA,
                        his_scales = "free",
                        transparency  = 0.05,
                        meshcolor = "black",
                        dotsize=1){

  #check type input
  if(type!="all"&type!="projection"&type!="histogram"&type!="length"&type!="width"){type <- readline("Please give type of plot: 'histogram', 'length', 'width', 'projection', or 'all' and press enter to confirm.")
    if(type!="all"&type!="projection"&type!="histogram"&type!="length"&type!="width"){stop("Plot type not recognized")}
  }
  #check by input
  if(by!="both"&by!="channel"&by!="condition"){by <- readline("Please give the correct factor for seperating the plots: 'channel', 'condition', or 'both' and press enter to confirm.")
    if(by!="both"&by!="channel"&by!="condition"){stop("Factor 'by' not recognized.")}
  }

  if(by=="both"|by=="channel"){

    if(missing(objectdata)!=T){
      if("channel"%in%colnames(objectdata)!=T){
        stop("No 'channel' column found in object dataframe.")
      }
    }
    if(missing(spotdata)!=T){
      if("channel"%in%colnames(spotdata)!=T){
        stop("No 'channel' column found in spot dataframe.")
      }
    }
  }

  if(by=="both"|by=="condition"){
    if(missing(meshdata)!=T){
      if("condition"%in%colnames(meshdata)!=T)
      {stop("No 'condition' column found in mesh dataframe.")
      }
    }
    if(missing(objectdata)!=T){
      if("condition"%in%colnames(objectdata)!=T){
        stop("No 'condition' column found in object dataframe.")
      }
    }
    if(missing(spotdata)!=T){
      if("condition"%in%colnames(spotdata)!=T){
        stop("No 'condition' column found in spot dataframe.")
      }
    }
  }


  #check mesh for requirements (if they're there)
  if(missing(meshdata)!=T){
    if("Xrotum"%in%colnames(meshdata)){
      meshdata$Xrot_micron <- meshdata$Xrotum
      meshdata$Yrot_micron <- meshdata$Yrotum
      meshdata$Xrotum <- NULL
      meshdata$Yrotum <- NULL
    }
    if("Xrot_micron"%in%colnames(meshdata)!=T){
      if(missing(mag)){
        if("Xrotum"%in%colnames(meshdata)!=T){
          stop("No pixel-micron-corrected X/Y ggplot2::coordinates in mesh file 'Xrot_micron' and 'Yrot_micron' found and no magnification correction value ('mag') found. Add value 'mag' in function or add columns 'Xrot_micron' & 'Yrot_micron' to meshdata")}

      }
      meshdata$Xrot_micron <- meshdata$X_rot * unlist(get(magnificationList,envir=magEnv)[mag])
      meshdata$Yrot_micron <- meshdata$Y_rot * unlist(get(magnificationList,envir=magEnv)[mag])
    }
  }


  #check objects for requirements
  if(missing(objectdata)!=T){
    if("ob_out_x"%in%colnames(objectdata)!=T){
      if("l"%in%colnames(objectdata)!=T){stop("No object outlines correlated to cell outlines found ('ob_out_x'/'ob_out_y' or 'l'/'d'). Please correlate outlines to cells using 'spotsInBox' before entering to this function")}
      if(missing(mag)){stop("No magnification correction value ('mag') found. Add value 'mag' in function or add columns 'ob_out_y' & 'ob_out_x' to objectdata")}
      if("l"%in%colnames(objectdata)==T){ objectdata$ob_out_x <- objectdata$l  * unlist(get(magnificationList,envir=magEnv)[mag])
                                objectdata$ob_out_y <- objectdata$d  * unlist(get(magnificationList,envir=magEnv)[mag])
      }
    }
  }

  #and checking the spotdata requirements
  if(missing(spotdata)!=T){
    if("Lmid"%in%colnames(spotdata)!=T){
      if("l"%in%colnames(spotdata)!=T){stop("Spots are not related to any cells yet. Please use 'spotsInBox' to add cell information and add pixel-to-micron converter 'mag' in this function.")}
      if(missing(mag)){stop("Please add pixel-to-micron converter 'mag' to the function to calculate the spot ggplot2::coordinates in micron.")}
      spotdata$Lmid <- spotdata$l * unlist(get(magnificationList,envir=magEnv)[mag])
      spotdata$Dum <- spotdata$d * unlist(get(magnificationList,envir=magEnv)[mag])
      }
  }

  #group data if needed
    if(by=="both"|by=="condition"){
      collist <- c("frame", "cell", quantiles_by, "condition")
    }
    if(by=="channel"){
      collist <- c("frame", "cell", quantiles_by)
    }
    if(missing(meshdata)!=T){
      onlycells <- unique(meshdata[,collist])
    }
    if(missing(meshdata)==T&missing(spotdata)!=T){
      onlycells <- unique(spotdata[,collist])
    }
    if(missing(meshdata)==T&missing(spotdata)==T){
      if(missing(objectdata==T)){stop("No data input found.")}
      onlycells <- unique(objectdata[,collist])
    }

  onlycells <- onlycells[order(onlycells[,quantiles_by]),]
  onlycells$number <- c(1:nrow(onlycells))

  if(quantiles>1&equal_groups==TRUE){
    onlycells$quantiles <- cut(onlycells$number, breaks=quantiles, labels = c(1:quantiles))
  }
  if(quantiles>1&equal_groups==FALSE){
    onlycells$quantiles <- cut(onlycells[,quantiles_by], breaks=quantiles, labels = c(1:quantiles))
  }

  onlycells <- onlycells[,colnames(onlycells)[colnames(onlycells)!=quantiles_by]]
  #build plot

  if(type=="all"){
    plotout <- list()
  }
  if(type=="projection"|type=="all"){
    plot <- ggplot2::ggplot() + ggplot2::theme_minimal()
    if(missing(meshdata)!=T){
      meshdata <- merge(meshdata, onlycells)
      meshdata$cellframe <- paste(meshdata$cell, meshdata$frame, sep="_")
      plot <- plot + ggplot2::geom_polygon(data=meshdata, ggplot2::aes_string(x='Xrot_micron', y='Yrot_micron', group='cellframe'),
                                           fill=meshcolor, alpha=transparency, color=NA)
    }
    if(missing(objectdata)!=T){
      objectdata <- merge(objectdata, onlycells)
      objectdata <- objectdata[order(objectdata$frame, objectdata$cell, objectdata$obpath),]
      objectdata$frameOB <- paste(objectdata$frame, objectdata$obID, sep="_")
      if(by=="channel"|by=="both"){
        plot <- plot + ggplot2::geom_polygon(data=objectdata, ggplot2::aes_string(x='ob_out_x', y='ob_out_y', fill='channel', group='frameOB'),
                                             alpha=transparency, color=NA) + ggplot2::scale_fill_manual(values=objectcolor)
      }
      if(by=="condition"){
        plot <- plot + ggplot2::geom_polygon(data=objectdata, ggplot2::aes_string(x='ob_out_x', y='ob_out_y', group='frameOB'),
                                             alpha=transparency, color=NA, fill=objectcolor[1])
      }
    }
    if(missing(spotdata)!=T){
      spotdata <- merge(spotdata, onlycells)
      if(by=="channel"|by=="both"){
        plot <- plot + ggplot2::geom_point(data=spotdata, ggplot2::aes_string(x='Lmid', y='Dum', color='channel'), size=dotsize, alpha=transparency*10, shape=16) + ggplot2::scale_color_manual(values=spotcolor)
      }
      if(by=="condition"){
        plot <- plot + ggplot2::geom_point(data=spotdata, ggplot2::aes_string(x='Lmid', y='Dum'), color=spotcolor[1], size=dotsize, alpha=transparency*10, shape=16)
      }
    }

    plot <- plot + ggplot2::coord_fixed()

    if(quantiles>1){
      if(by=="condition"|by=="both"){
        plot <- plot + ggplot2::facet_grid(quantiles~condition)
      }
      if(by=="channel"){
        plot <- plot + ggplot2::facet_grid(quantiles~.)
      }
    }
    if(quantiles<=1){
      if(by=="condition"|by=="both"){
        plot <- plot + ggplot2::facet_grid(.~condition)
      }
    }
    if(type =="all"){
      plotout$projection <- plot
    }

  }


  if(type=="length"|type=="all"){
    plot <- ggplot2::ggplot() + ggplot2::theme_minimal()
    if(missing(objectdata)!=T){
      if(type!="all"){
        objectdata <- merge(objectdata, onlycells)
      }
      objectdata$frameOB <- paste(objectdata$frame, objectdata$obID, sep="_")
      if(by=="both"|by=="channel"){
        plot <- plot + ggplot2::geom_path(data=objectdata, ggplot2::aes_string(x='number', y='ob_out_x', color='channel', group='frameOB'), alpha=0.5) + ggplot2::scale_color_manual(values=objectcolor)

      }

      if(by=="condition"){
        plot <- plot + ggplot2::geom_path(data=objectdata, ggplot2::aes_string(x='number', y='ob_out_x'), color=objectcolor[1], alpha=0.5)
      }
    }
    if(missing(spotdata)!=T){
      if(type!="all"){
        spotdata <- merge(spotdata,onlycells)
      }
      if(by=="both"|by=="channel"){
        plot <- plot + ggplot2::geom_point(data=spotdata, ggplot2::aes_string(x='number', y='Lmid', color='channel'), alpha=0.5) + ggplot2::scale_color_manual(values=spotcolor)
      }
      if(by=="condition"){
        plot <- plot + ggplot2::geom_point(data=spotdata, ggplot2::aes_string(x='number', y='Lmid'), color=spotcolor[1], alpha=0.5)
      }
    }
    if(missing(spotdata)==T&missing(objectdata)==T&type=="length"){
      stop("Only mesh outlines found as input so nothing to correlate by length. Please select 'all', 'projection' or 'histogram' as plotting option")
    }
    if((missing(spotdata)!=T|missing(objectdata)!=T)){
      if(by=="condition"|by=="both"){
        plot <- plot + ggplot2::facet_grid(.~condition)
      }
      if(type=="all"){
        plotout$length <- plot
      }
    }
  }

  if(type=="width"|type=="all"){
    plot <- ggplot2::ggplot() + ggplot2::theme_minimal()
    if(missing(objectdata)!=T){
      if(type!="all"){
        objectdata <- merge(objectdata, onlycells)
      }
      objectdata$frameOB <- paste(objectdata$frame, objectdata$obID, sep="_")
      if(by=="both"|by=="channel"){
        plot <- plot + ggplot2::geom_path(data=objectdata, ggplot2::aes_string(x='number', y='ob_out_y', color='channel', group='frameOB'), alpha=0.5) +
                                            ggplot2::scale_color_manual(values=objectcolor)
      }
      if(by=="condition"){
        plot <- plot + ggplot2::geom_path(data=objectdata, ggplot2::aes_string(x='number', y='ob_out_y'), color=objectcolor[1], alpha=0.5)
      }
    }
    if(missing(spotdata)!=T){
      if(type!="all"){
        spotdata <- merge(spotdata, onlycells)
      }
      if(by=="both"|by=="channel"){
        plot <- plot + ggplot2::geom_point(data=spotdata, ggplot2::aes_string(x='number', y='Dum', color='channel'), alpha=0.5)
        if(missing(objectdata)==T){
        plot <- plot + ggplot2::scale_color_manual(values=spotcolor)
        }
      }
      if(by=="condition"){
        plot <- plot + ggplot2::geom_point(data=spotdata, ggplot2::aes_string(x='number', y='Dum'), color=spotcolor[1], alpha=0.5)
      }
    }
    if(missing(spotdata)==T&missing(objectdata)==T&type=="width"){
      stop("Only mesh outlines found as input so nothing to correlate by width. Please select 'all', 'projection' or 'histogram' as plotting option")
    }
    if((missing(spotdata)!=T|missing(objectdata)!=T)){
      if(by=="condition"|by=="both"){
        plot <- plot + ggplot2::facet_grid(.~condition)
      }
      if(type=="all"){
        plotout$width <- plot
      }
    }
  }

  if(type=="histogram"|type=="all"){
    histograms <- list()
    if(missing(meshdata)!=T){

      if(by=="channel"){
        meshdata <- unique(meshdata[,c("cell", "frame", quantiles_by)])
      }
      if(by=="condition"|by=="both"){
        meshdata <- unique(meshdata[,c("cell", "frame", quantiles_by, "condition")])
      }
      meshdata <- merge(meshdata, onlycells)
      meshdata$quant_by <- meshdata[,quantiles_by]
      his_mesh <- ggplot2::ggplot(meshdata) +
        ggplot2::geom_density(ggplot2::aes_string(x='quant_by'), fill=meshcolor, color=histogram_outline) + ggplot2::theme_minimal()
      if(by!="channel"){
        his_mesh <- his_mesh + ggplot2::facet_grid(.~condition, scales = his_scales)
      }


      histograms$mesh <- his_mesh
    }

    if(missing(objectdata)!=T){

      if(by=="channel"){
        objectdata <- unique(objectdata[,c("cell", "frame", "Lmid", "channel")])
      }
      if(by=="condition"){
        objectdata <- unique(objectdata[,c("cell", "frame", "Lmid", "condition")])
      }
      if(by=="both"){
        objectdata <- unique(objectdata[,c("cell", "frame", "Lmid", "condition", "channel")])
      }

      objectdata <- merge(objectdata, onlycells)


      if(by=="channel"|by=="both"){
        his_obj <- ggplot2::ggplot(objectdata) +
          ggplot2::geom_density(ggplot2::aes_string(x='Lmid', fill='channel'), color=histogram_outline, alpha=0.5) + ggplot2::theme_minimal()
      }
      if(by=="condition"){
        his_obj <- ggplot2::ggplot(objectdata) + ggplot2::geom_density(ggplot2::aes_string(x='Lmid'), fill="black", color=histogram_outline)
      }
      if(by!="channel"&quantiles>1){
        his_obj <- his_obj + ggplot2::facet_grid(quantiles~condition, scales = his_scales)
      }
      if(by=="channel"&quantiles>1){
        his_obj <- his_obj + ggplot2::facet_grid(quantiles~., scales = his_scales)
      }
      histograms$objects <- his_obj
    }

    if(missing(spotdata)!=T){

      if(by=="channel"){
        spotdata <- unique(spotdata[,c("cell", "frame", "Lmid", "channel")])
      }
      if(by=="condition"){
        spotdata <- unique(spotdata[,c("cell", "frame", "Lmid", "condition")])
      }
      if(by=="both"){
        spotdata <- unique(spotdata[,c("cell", "frame", "Lmid", "condition", "channel")])
      }
      spotdata<- merge(spotdata,onlycells)


      if(by=="channel"|by=="both"){
        his_spot <- ggplot2::ggplot(spotdata) +
          ggplot2::geom_density(ggplot2::aes_string(x='Lmid', fill='channel'), alpha=0.5, color=histogram_outline) +
          ggplot2::theme_minimal() +
          ggplot2::scale_fill_manual(values=spotcolor)
      }
      if(by=="condition"){
        his_spot <- ggplot2::ggplot(spotdata) +
          ggplot2::geom_density(ggplot2::aes_string(x='Lmid'), fill="black", color=histogram_outline)
      }
      if(by!="channel"&quantiles>1){
        his_spot <- his_spot + ggplot2::facet_grid(quantiles~condition, scales = his_scales)
      }
      if(by=="channel"&quantiles>1){
        his_spot <- his_spot + ggplot2::facet_grid(quantiles~. , scales = his_scales)
      }
      histograms$spots <- his_spot
    }

    if(type=="all"){
      plotout$histograms <- histograms
    }
    if(type=="histogram"){
      plot <- histograms
    }

  }

  if(type=="all"){return(plotout)}
  if(type!="all"){return(plot)}
}





#3. batch upload (mesh=c(filepath, "type"), spots="same" or c(filepath, "type"), object etc)
##for this: ask users to order their files according to condition/color/analysis type.
