##24.7.2018
##Renske van Raaphorst

##function to plot object projections (now only possible with oufti output)

#' @export
plotObjects <- function(obdat, meshdat, groups=1, cellcolor="black", objectcolor="green", transparency=0.1, getdata=FALSE){
  if(getdata==TRUE){outlist <- list()}

  if(groups>1){
    obcellsonly <- unique(obdat[,c("cell","max.length", "frame")])
    obcellsonly <- obcellsonly[order(obcellsonly$max.length),]
    obcellsonly$lnum <- c(1:nrow(obcellsonly))
    obcellsonly$q1 <- cut(obcellsonly$lnum, breaks=groups, labels=c(1:groups))
    obdat <- merge(obcellsonly, obdat, all=T)
  }
  if(groups==1){
    obdat$q1 <- 1
  }

  obdat <- obdat[order(obdat$frame, obdat$cell, obdat$obnum, obdat$obpath),]
  obdat$frameOB <- paste(obdat$frame, obdat$obID, sep="_")
  if(getdata==TRUE){outlist$object_data <- obdat}
  if(missing(meshdat)==TRUE){
    outplot <- ggplot2::ggplot(obdat, ggplot2::aes_string(x='ob_out_x', y='ob_out_y', group='frameOB')) +
      ggplot2::geom_polygon(alpha=transparency, fill=objectcolor, color=NA) +
      ggplot2::facet_grid(q1~.) +
      ggplot2::coord_fixed() +
      ggplot2::theme_classic() +
      ggplot2::xlab("X coordinates (micron)") +
      ggplot2::ylab("Y coordinates (micron)")

    if(getdata==TRUE){outlist$objectplot <- outplot}
  }
  if(missing(meshdat)!=TRUE){
    meshdat <- merge(meshdat, obdat[,c("cell", "frame", "q1")])
    meshdat <- meshdat[order(meshdat$frame, meshdat$cell, meshdat$num),]
    meshdat$cellframe <- paste(meshdat$cell, meshdat$frame, sep="_")
    obdat <- obdat[order(obdat$frame, obdat$cell, obdat$obnum, obdat$obpath),]
    outplot <- ggplot2::ggplot(meshdat) +
      ggplot2::geom_polygon(ggplot2::aes_string(x='Xrot_micron', y='Yrot_micron', group='cellframe'), alpha=transparency, fill=cellcolor, color=NA) +
      ggplot2::geom_polygon(data=obdat, ggplot2::aes_string(x='ob_out_x', y='ob_out_y', group='frameOB'), alpha=transparency, fill=objectcolor, color=NA) +
      ggplot2::facet_grid(q1~.) +
      ggplot2::coord_fixed() +
      ggplot2::theme_classic() +
      ggplot2::xlab("X coordinates (micron)") +
      ggplot2::ylab("Y coordinats (micron)")

    if(getdata==TRUE){
      outlist$objectplot <- outplot
      outlist$mesh_data <- meshdat
    }

  }
  if(getdata==TRUE){return(outlist)}
  if(getdata==FALSE){return(outplot)}
}
