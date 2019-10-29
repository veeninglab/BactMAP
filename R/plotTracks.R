#plot tracks
#12.18.18
#R van Raaphorst
#bactMAP

#'@export
plotTracks <- function(meshdata,
                       spotdata,
                       objectdata,
                       tracks=FALSE,
                       ignore_singles=FALSE,
                       movie=FALSE,
                       timepalette_lines="viridis",
                       timepalette_fill="magma",
                       dimension = c("length", "width"),
                       turn_cells = TRUE,
                       mag,
                       cell = "all",
                       transparency = 0.2
){
  outplot <- ggplot2::ggplot() + ggplot2::theme_minimal()
  if(cell[1]!="all"&is.numeric(cell)==FALSE){stop("Cell number not recognized. Give cell number or keep default (cell='all') to return all cells")}

  if(missing(meshdata)!=TRUE){
    if(cell[1]!="all"){
      meshdata <- meshdata[meshdata$cell%in%cell,]
    }
    if("length"%in%dimension&"width"%in%dimension){

      if(turn_cells==TRUE){
        if("X_rot"%in%colnames(meshdata)){
          if("Xrotum"%in%colnames(meshdata)){
            meshdata$Xrot_micron <- meshdata$Xrotum
            meshdata$Yrot_micron <- meshdata$Yrotum
            meshdata$Xrotum <- NULL
            meshdata$Yrotum <- NULL
          }
          if("Xrot_micron"%in%colnames(meshdata)!=T){
            if(missing(mag)){
              stop("Please specify the pixel-micron conversion value 'mag'.")
            }
            meshdata$Xrot_micron <- meshdata$X_rot*unlist(get(magnificationList,envir=magEnv)[mag])
            meshdata$Yrot_micron <- meshdata$Y_rot*unlist(get(magnificationList,envir=magEnv)[mag])
          }
        }
        outplot <- outplot +
          ggplot2::geom_polygon(data=meshdata,
                                ggplot2::aes_string(x='Xrot_micron', y='Yrot_micron', color='frame', group='frame'),
                                fill="black", alpha=transparency/5)
      }
      if(turn_cells==FALSE){
        outplot <- outplot +
          ggplot2::geom_polygon(data=meshdata,
                                ggplot2::aes_string(x='X', y='Y', color='frame', group='frame'),
                                fill="black", alpha=transparency/5)
      }
    }
  }

  if(missing(objectdata)!=TRUE){
    if(cell[1]!="all"){
      objectdata <- objectdata[objectdata$cell%in%cell,]
    }
    if(turn_cells==TRUE|("length"%in%dimension&!"width"%in%dimension)|(!"length"%in%dimension&"width"%in%dimension)){
      if("ob_out_x"%in%colnames(objectdata)!=TRUE){
        if("ob_x"%in%colnames(objectdata)!=TRUE){
          stop("Cannot find object cobjectdatardinate data (ob_out_x/ob_out_y for turned cells, ob_x/ob_y for raw cobjectdatardinates.")

          }
        if(missing(meshdata)){
            stop("Cannot find turned object data or mesh information to connect the objects to the cell localizations.
                     Replace object dataframe for a dataframe (object_relative) containing this information or add mesh data to convert dataframe.")
          }
        if(missing(mag)){
            stop("Need pixel to micron conversion factor 'mag' to correctly connect mesh to object data.")
          }

          objectdata <- suppressWarnings(centrefun(objectdata))
          objectdata <- suppressWarnings(midobject(meshdata, objectdata, get(magnificationList,envir=magEnv)[mag]))
      }
      objectdata$frameOB <- paste(objectdata$frame, objectdata$obID, sep="_")

      if("length"%in%dimension&"width"%in%dimension){

        outplot <- outplot +
          ggplot2::geom_polygon(data=objectdata,
                                ggplot2::aes_string(x='ob_out_x',y='ob_out_y', group='frameOB', fill='frame'),
                                color="black", alpha=transparency)
      }
      if("length"%in%dimension&!"width"%in%dimension){
        Om <- stats::aggregate(objectdata[,c("ob_out_x", "pole1", "pole2")], by=list("cell"=objectdata$cell, "frame"=objectdata$frame, "obID"=objectdata$obID), FUN=max)
        Omin <- stats::aggregate(objectdata[,c("ob_out_x", "pole1", "pole2")], by=list("cell"=objectdata$cell, "frame"=objectdata$frame, "obID"=objectdata$obID), FUN=min)
        Om$obxmax <- Om$ob_out_x
        Om$ob_out_x <- NULL
        Omin$obxmin <- Omin$ob_out_x
        Omin$ob_out_x <- NULL
        objectdata <- merge(Om, Omin)

        objectdata$framemin <- objectdata$frame - 0.1
        objectdata$framemax <- objectdata$frame + 0.1
        outplot <- outplot +
          ggplot2::geom_rect(data=objectdata,
                             ggplot2::aes_string(xmin = 'obxmin', xmax='obxmax', ymin='framemin', ymax = 'framemax', group='frameOB', fill='frame'))
        if(missing(spotdata)){
          objectdata <- objectdata[order(objectdata$cell, objectdata$frame),]
          outplot <- outplot +
            ggplot2::geom_path(data=objectdata, ggplot2::aes_string(x='pole1', y='frame')) +
            ggplot2::geom_path(data=objectdata, ggplot2::aes_string(x='pole2', y='frame'))
          outplot <- outplot + ggplot2::xlab("location on cell length axis (\u00b5m)") + ggplot2::ylab("Time (Frames)")
        }
      }
      if(!"length"%in%dimension&"width"%in%dimension){
        Om <- stats::aggregate(objectdata[,c("ob_out_y", "maxwum")], by=list("cell"=objectdata$cell, "frame"=objectdata$frame, "obID"=objectdata$obID), FUN=max)
        Omin <- stats::aggregate(objectdata[,c("ob_out_y", "maxwum")], by=list("cell"=objectdata$cell, "frame"=objectdata$frame, "obID"=objectdata$obID), FUN=min)
        Om$obymax <- Om$ob_out_y
        Om$ob_out_y <- NULL
        Omin$obymin <- Omin$ob_out_y
        Omin$ob_out_y <- NULL
        objectdata <- merge(Om, Omin)

        objectdata$framemin <- objectdata$frame - 0.1
        objectdata$framemax <- objectdata$frame + 0.1

        outplot <- outplot +
          ggplot2::geom_rect(data=objectdata,
                             ggplot2::aes_string(xmin = 'obymin', xmax='obymax', ymin='framemin', ymax = 'framemax', group='frameOB', fill='frame'))

        if(missing(spotdata)){
          objectdata$pole1 <- objectdata$maxwum/2
          objectdata$pole2 <- objectdata$pole1*-1
          objectdata <- objectdata[order(objectdata$cell, objectdata$frame),]

          outplot <- outplot +
            ggplot2::geom_path(data=objectdata, ggplot2::aes_string(x='pole1', y='frame')) +
            ggplot2::geom_path(data=objectdata, ggplot2::aes_string(x='pole2', y='frame'))

          outplot <- outplot +
            ggplot2::xlab("location on cell length axis (\u00b5m)") +
            ggplot2::ylab("Time (Frames)")}
      }
    }
    if(turn_cells!=TRUE&"length"%in%dimension&"width"%in%dimension){
      if("ob_x"%in%colnames(objectdata)!=TRUE){
        stop("Cannot find object coordinates 'ob_x'/'ob_y'.")
      }
      objectdata$frameOB <- paste(objectdata$frame, objectdata$obID, sep="_")
      outplot <- outplot +
        ggplot2::geom_polygon(data=objectdata, ggplot2::aes_string(x='ob_x', y='ob_y', group='frameOB', fill='frame'), color="black", alpha=transparency)
    }

    outplot <- outplot + ggplot2::scale_fill_viridis_c(option=timepalette_fill)
  }

  if(missing(spotdata)!=TRUE){
    if(cell[1]!="all"){
      spotdata <- spotdata[spotdata$cell%in%cell,]
    }
    if(tracks==TRUE){
      if(ignore_singles==TRUE){spotdata <- spotdata[spotdata$trajectory!=-1,]}}
    if(turn_cells==TRUE|("length"%in%dimension&!"width"%in%dimension)|(!"length"%in%dimension&"width"%in%dimension)){
      if("Lmid"%in%colnames(spotdata)!=T){
        if("l"%in%colnames(spotdata)!=T){
          if(missing(meshdata)!=T){
            A <- readline("Data doesn't include relative spot positions. Press 'Y'+enter to start spotsInBox() to relate spot positions to cells or any other key to stop the function.")
            if(A=="Y"|A=="y"){
              spotdata <- spotsInBox(spotdata, meshdata)$spots_relative
            }else{stop("Function stopped.")}
          }
        }
        if(missing(mag)){stop("Please specify the pixel-micron conversion value 'mag'.")}
        spotdata$Lmid <- spotdata$l * unlist(get(magnificationList,envir=magEnv)[mag])
        spotdata$Dum <- spotdata$d * unlist(get(magnificationList,envir=magEnv)[mag])
      }
      if("length"%in%dimension&"width"%in%dimension){
        outplot <- outplot + ggplot2::geom_point(data=spotdata, ggplot2::aes_string(x='Lmid', y='Dum', color='frame'), size=2)
        if(tracks==TRUE){
          outplot <- outplot + ggplot2::geom_path(data=spotdata[spotdata$trajectory!=-1,], ggplot2::aes_string(x='Lmid', y='Dum', group='trajectory', color='frame'), size=1)
        }
      }
      if("length"%in%dimension&!"width"%in%dimension){
        if("channel"%in%colnames(spotdata)){
          outplot <- outplot + ggplot2::geom_point(data=spotdata, ggplot2::aes_string(x='Lmid', y='frame', color='channel'), size=2)
        }
        if(!"channel"%in%colnames(spotdata)){
          outplot <- outplot + ggplot2::geom_point(data=spotdata, ggplot2::aes_string(x='Lmid', y='frame', color='frame'), size=2)
        }
        if(tracks==TRUE){
          if(!"channel"%in%colnames(spotdata)){
            outplot <- outplot + ggplot2::geom_path(data=spotdata[spotdata$trajectory!=-1,], ggplot2::aes_string(x='Lmid', y='frame', group='trajectory', color='frame'), size=1)
          }
          if("channel"%in%colnames(spotdata)){
            outplot <- outplot + ggplot2::geom_path(data=spotdata[spotdata$trajectory!=-1,], ggplot2::aes_string(x='Lmid', y='frame', group='trajectory', color='channel'), size=1)
          }
        }

        if(!"pole1"%in%colnames(spotdata)){
          if(!"maxum"%in%colnames(spotdata)){
            spotdata$maxum <- spotdata$max.length*unlist(get(magnificationList,envir=magEnv)[mag])
          }
          spotdata$pole1 <- 0.5*spotdata$maxum
          spotdata$pole2 <- -spotdata$pole1
          spotdata <- spotdata[order(spotdata$cell, spotdata$frame),]
          outplot <- outplot + ggplot2::geom_path(data=spotdata, ggplot2::aes_string(x='pole1', y='frame')) + ggplot2::geom_path(data=spotdata, ggplot2::aes_string(x='pole2', y='frame'))
        }
      }

      if(!"length"%in%dimension&"width"%in%dimension){
        if("channel"%in%colnames(spotdata)){
          outplot <- outplot + ggplot2::geom_point(data=spotdata, ggplot2::aes_string(x='Dum', y='frame', color='channel'), size=2)
        }
        if(!"channel"%in%colnames(spotdata)){
          outplot <- outplot + ggplot2::geom_point(data=spotdata, ggplot2::aes_string(x='Dum', y='frame', color='frame'), size=2)
        }
        if(tracks==TRUE){
          if(!"channel"%in%colnames(spotdata)){
            outplot <- outplot + ggplot2::geom_path(data=spotdata[spotdata$trajectory!=-1,], ggplot2::aes_string(x='Dum', y='frame', group='trajectory', color='frame'), size=1)
          }
          if("channel"%in%colnames(spotdata)){
            outplot <- outplot + ggplot2::geom_path(data=spotdata[spotdata$trajectory!=-1,], ggplot2::aes_string(x='Dum', y='frame', group='trajectory', color='channel'), size=1)

          }
        }
        if(!"pole1"%in%colnames(spotdata)){
          if(!"maxwum"%in%colnames(spotdata)){
            spotdata$maxwum <- spotdata$max.width*unlist(get(magnificationList,envir=magEnv)[mag])
          }
          spotdata$pole1 <- 0.5*spotdata$maxwum
          spotdata$pole2 <- -spotdata$pole1
          spotdata <- spotdata[order(spotdata$cell, spotdata$frame),]
          outplot <- outplot + ggplot2::geom_path(data=spotdata, ggplot2::aes_string(x='pole1', y='frame')) + ggplot2::geom_path(data=spotdata, ggplot2::aes_string(x='pole2', y='frame'))
        }

      }
    }
      if(turn_cells==FALSE&"length"%in%dimension&"width"%in%dimension){
        outplot <- outplot + ggplot2::geom_point(data=spotdata, ggplot2::aes_string(x='x',y='y', color='frame'), size=2)
        if(tracks==TRUE){
          outplot <- outplot + ggplot2::geom_path(data=spotdata[spotdata$trajectory!=-1,], ggplot2::aes_string(x='x',y='y',group='trajectory', color='frame'),size=1)
        }

      }
    }

    if(missing(meshdata)!=T|missing(spotdata)!=T){
      if(("length"%in%dimension&!"width"%in%dimension)|(!"length"%in%dimension&"width"%in%dimension)){
        if(!"channel"%in%colnames(spotdata)){
          outplot <- outplot + ggplot2::scale_color_viridis_c(option=timepalette_lines)
        }
        if("channel"%in%colnames(spotdata)){
          outplot <- outplot + ggplot2::scale_color_manual(values=colorsCUD())
        }
      }
      if("length"%in%dimension&"width"%in%dimension){
        outplot <- outplot + ggplot2::scale_color_viridis_c(option = timepalette_lines)
      }
    }
    if(turn_cells==TRUE&"length"%in%dimension&"width"%in%dimension){
      if("channel"%in%colnames(spotdata)&"condition"%in%colnames(spotdata)){
        outplot <- outplot + ggplot2::coord_fixed() + ggplot2::facet_grid(channel~condition~cell)
      }
      if("channel"%in%colnames(spotdata)&!"condition"%in%colnames(spotdata)){
        outplot <- outplot + ggplot2::coord_fixed() + ggplot2::facet_grid(channel~cell)
      }
      if(!"channel"%in%colnames(spotdata)&"condition"%in%colnames(spotdata)){
        outplot <- outplot + ggplot2::coord_fixed() + ggplot2::facet_grid(condition~cell)
      }
      if(!"channel"%in%colnames(spotdata)&!"condition"%in%colnames(spotdata)){
        outplot <- outplot + ggplot2::coord_fixed() + ggplot2::facet_grid(~cell)
      }
      outplot <- outplot + ggplot2::xlab(paste("x (", micron(), ")")) + ggplot2::ylab(paste("y (", micron(), ")"))
    }
    if(turn_cells==FALSE|("length"%in%dimension&!"width"%in%dimension)|(!"length"%in%dimension&"width"%in%dimension)){
      if("condition"%in%colnames(spotdata)){
        outplot <- outplot + ggplot2::facet_grid(condition~cell, scales="free")
      }
      if(!"condition"%in%colnames(spotdata)){
        outplot <- outplot + ggplot2::facet_grid(~cell, scales="free")
      }
      if("length"%in%dimension&"width"%in%dimension){
        outplot <- outplot + ggplot2::xlab("x (pixels)") + ggplot2::ylab("y (pixels)")
      }
      if("length"%in%dimension&!"width"%in%dimension){
        outplot <- outplot + ggplot2::xlab(paste("distance from mid-cell on length axis (", micron(), ")"))
      }
      if(!"length"%in%dimension&"width"%in%dimension){
        outplot <- outplot + ggplot2::xlab(paste("distance from mid-cell on width axis (", micron(), ")"))
      }
    }
    return(outplot)
  }


