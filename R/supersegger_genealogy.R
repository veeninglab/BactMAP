##SuperSegger Input for BactMAP


##16/05/2018

##uses packages R.matlab, igraph, ape, ggplot2, ggtree & treeio


#' @export
extr_SuperSeggerClist <- function(matfile, trim.orphans=TRUE, cellList=FALSE){
  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    inp <- readline("Packages 'R.matlab' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){utils::install.packages("R.matlab")}else{stop("Canceled")}
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    inp <- readline("Packages 'igraph' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){utils::install.packages("igraph")}else{stop("Canceled")}
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    inp <- readline("Packages 'ape' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){utils::install.packages("ape")}else{stop("Canceled")}
  }
  clist <- R.matlab::readMat(matfile)
  datasegger <- as.data.frame(clist$data)
  colnames(datasegger) <- unlist(clist$def)
  out <- list()
  out$cellList <- datasegger
  datasegger <- datasegger[,c("Cell ID", "Cell birth time", "Cell death time", "Cell age", "Fluor1 sum", "Fluor1 mean", "Fluor1 sum death", "Fluor1 mean death", "Mother ID", "Daughter1 ID", "Daughter2 ID")]
  colnames(datasegger) <- c("cell", "birth", "death", "edgelength", "fluorsum", "fluormean", "fluorsum_D", "fluormean_D", "parent", "child1", "child2")
  #datasegger$cell <- datasegger$cell + 1
  #datasegger$parent <- datasegger$parent + 1
  #datasegger$child1 <- datasegger$child1 + 1
  #datasegger$child2 <- datasegger$child2 + 1
  datasegger$parent[is.na(datasegger$parent)] <- 0
  rownames(datasegger) <- NULL
  if(trim.orphans==TRUE){
    out$orphans <- (datasegger[datasegger$parent==0&is.na(datasegger$child1)&is.na(datasegger$child2),])
    datasegger  <- datasegger[datasegger$parent!=0|is.na(datasegger$child1)!=TRUE,]
    print(paste("dataset trimmed. percentage of cells without offspring/parents: ", round(nrow(out$orphans)/(nrow(datasegger)+nrow(out$orphans))*100, digits=0)))
  }
  net <- igraph::graph_from_data_frame(datasegger[,c(9,1,2,3,4,5,6,7,8,10,11)])
  datasegger$root <- 0
  out$network <- net
  phylos <- getphylolist_SupSeg(datasegger)
  out$generation_lists <- phylos$generation_lists
  out$data_generation_dataframes <- phylos$generation_dataframes
  if(trim.orphans==TRUE){
    out$cellList_trimmed <- phylos$genframe
  }
  if(cellList==FALSE){
    out$cellList<-NULL
  }
  message(summary(out))
  return(out)
}


#' @importFrom utils install.packages
#' @importFrom ggtree %<+%
#' @export
plotTreeBasic <- function(phylo, extradata, yscalechange = FALSE, showClade = FALSE, layout = "rectangular", ydata, cellNumber, open.angle, linesize = 1, linecolor = "black", lines=TRUE, colors=FALSE){
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    inp <- readline("Packages 'ggtree' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if(inp=="y"|inp=="Y"){
        if (!requireNamespace("BiocManager", quietly = TRUE))utils::install.packages("BiocManager")
      BiocManager::install("ggtree", version = "3.8")
    }
  }
  if(showClade==TRUE){
    if(missing(cellNumber)){stop("cellNumber missing. Please state the ancestor cell you want to follow to show it's clade.")}
    NodeNumber <- extradata$node[extradata$cell==cellNumber]
    phylo <- ggtree::groupClade(phylo, .node=NodeNumber)
    if(lines==T&colors==F){
      gP <- ggtree::ggtree(phylo, layout=layout, open.angle=open.angle, ggplot2::aes_string(linetype='group'), size=linesize, color=linecolor)
    }
    if(lines==F&colors==T){
      gP <- ggtree::ggtree(phylo, layout=layout, open.angle=open.angle, ggplot2::aes_string(color='group'), size=linesize)
    }
    if(lines==T&colors==T){
      gP <- ggtree::ggtree(phylo, layout=layout, open.angle=open.angle, ggplot2::aes_string(linetype='group', color='group'), size=linesize)
    }
  }
  if(showClade!=TRUE){gP <- ggtree::ggtree(phylo, layout=layout, open.angle=open.angle, size=linesize, color=linecolor)}
  if(!missing(extradata)){
  gP <- gP %<+% extradata
  }
  if(yscalechange==TRUE){
    if(missing(ydata)){stop("Variable 'ydata' missing. Don't know what values to put on the y axis. Please give the column name of your data as 'ydata' in the function")}
    gP$data$y <- gP$data[,ydata]
    gP <- gP + ggplot2::theme_classic() + ggplot2::xlab("Time") +  ggplot2::ylab(ydata)
  }
  return(gP)
}

getphylolist_SupSeg <- function(CLT, prep=FALSE){
  #if(prep==TRUE){
   # CLT <- prepcellListtree(CLT)
  #}
  phylolist <- list()
  fulldatlist <- list()
  z <- 0
  for(n in unique(CLT$root[!is.na(CLT$root)])){
    onephyl <- CLT[CLT$root==n&!is.na(CLT$parent),]
    if(nrow(onephyl)>1){
      z <- z+1
      Ntip <- length(onephyl$cell[onephyl$cell%in%onephyl$parent!=T])
      if(Ntip>1){
        onephyl$edge2[onephyl$cell%in%onephyl$parent!=T] <- c(1:Ntip)
        onephyl$edge2[onephyl$cell%in%onephyl$parent==T] <- c((Ntip+2):(nrow(onephyl)+1))
        onephyl$edge1 <- lapply(onephyl$parent, function(x) onephyl$edge2[onephyl$cell==x])
        onephyl$edge1[is.na(onephyl$edge1==0)] <- Ntip+1
        onephyl$edge1 <- as.numeric(as.character(onephyl$edge1))
        onephyltree <- list()
        onephyltree$edge <- matrix(c(as.integer(onephyl$edge1), as.integer(onephyl$edge2)), nrow(onephyl), 2)
        onephyltree$tip.label <- as.character(onephyl$cell[onephyl$cell%in%onephyl$parent!=T])
        onephyltree$Nnode <- max(onephyltree$edge) - Ntip
        onephyltree$root.label <- as.character(n)
        onephyltree$node.label <- c(as.character(n), as.character(onephyl$cell[onephyl$cell%in%onephyl$parent==T]))
        onephyltree$edge.length <- onephyl$edgelength
        onephyl$node <- onephyl$edge2
        onephyl$nodelabel <- onephyl$parent
        onephyl$parent <- onephyl$edge1
        onephyl$edge1 <- NULL
        onephyl$edge2 <- NULL
        q <-data.frame(x = colnames(onephyl)=="node", y = 1:ncol(onephyl))
        nodecol <- q$y[q$x==TRUE]
        colone <- c(nodecol, q$y[q$x==FALSE])
        onephyl <- onephyl[,colone]
        fulldatlist[[z]] <- onephyl
        onephyltree$root.edge <- 0
        class(onephyltree) <- "phylo"
        onephyltree <- ape::collapse.singles(onephyltree)
        phylolist[[z]] <- onephyltree
      }
    }
  }
  fulldatframe <- do.call('rbind', fulldatlist)
  if(length(phylolist)>1){
    class(phylolist) <- "multiPhylo"
  }else{
    phylolist <- phylolist[[1]]
    fulldatlist <- fulldatlist[[1]]
  }
  return(list(generation_lists = phylolist, generation_dataframes=fulldatlist, genframe = fulldatframe))
}
