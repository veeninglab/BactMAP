\name{extr_SuperSeggerCells}
\alias{extr_SuperSeggerCells}
\alias{extr_SuperSegger}
\alias{extr_SuperSeggerClist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert SuperSegger cell outlines & genealogy to standard BactMAP format
}
\description{
Part of the *extr_*-functions, extr_SuperSeggerCells takes the mask created by SuperSegger of each cell (and if available, timepoint) and saves this as a *mesh* dataframe.

extr_SuperSeggerClist() takes the clist (SuperSegger summarized output) and forms genealogy files which can be processed to network plots and trees.
}
\usage{
extr_SuperSeggerCells(loc, frames, mag, timelapse=FALSE, startframe=0, cellList=FALSE)

extr_SuperSeggerClist(matfile, trim.orphans=TRUE, cellList=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

For extr_SuperSeggercells:

  \item{loc}{
The file path where SuperSegger's output files are located
}
  \item{frames}{
The number of xy frames (so not timepoints, just different locations!) imaged
}
  \item{mag}{
magnification conversion factor name (which is part of Pixels2um)
}
 \item{timelapse}{
 set timelapse to TRUE if you are analyzing a timelapse movie.
 }
 \item{startframe}{
 default = 0. set to 1 if the first of your xy-locations is 1.
 }
 \item{cellList}{
 Default=FALSE. When TRUE, the cellList (see below) will be part of the output of the function.
 }

For extr_SuperSeggerClist:

 \item{matfile}{
 Path to the Clist output of supersegger (save this manually from the SuperSegger Viewer after analysis).
 }
 \item{trim.orphans}{
 Default=TRUE. When TRUE, cells without offspring or parents will be removed from the dataset.
 }
 \item{cellList}{
 Default=FALSE. When TRUE, the cellList (see below) will be part of the output of the function.
 }

}
\details{
Use *addPixels2um()* to add a new conversion factor.
}
\value{

From extr_SuperSeggerCells():

 \item{cellList}{dataframe with content similar to SuperSegger's "Clist"}
 \item{mesh}{dataframe containing cell coordinates and dimensions}
 \item{pixel2um}{the magnification conversion factor used}


From extr_SuperSeggerClist():

  \item{cellList}{copy of 2-dimensional cList output of SuperSegger}
  \item{network}{igraph network data of cell genealogy }
  \item{generation_lists}{phylo object describing cell genealogy}
  \item{data_generation_dataframes}{attribute data - dataframe with information on cell fluorescence, etc}
  \item{cellList_trimmed}{cellList without the "orphan" cells. in output when orphans==TRUE}
  \item{orphans}{the cells which are removed from the cellList when orphans==TRUE}

}
\references{
Stylianidou, Stella, et al. "SuperSegger: robust image segmentation, analysis and lineage tracking of bacterial cells." Molecular microbiology 102.4 (2016): 690-700.
}
\author{
Renske van Raaphorst
}

\examples{
\dontrun{
##pick supersegger output path. in this example there are 6 x/y frames.
out_supseg <- "path_to_supersegger_data"
##get cell outlines
cell_outlines <- extr_SuperSeggerCells( paste(out_supseg, "/clist.mat", sep=""),
                                        frames=6, startframe=0)
##get clist timelapse:
cell_clist <- extr_SuperSeggerClist(out_supseg)
}
}
