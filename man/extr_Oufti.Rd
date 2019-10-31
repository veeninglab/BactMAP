\name{extr_Oufti}
\alias{extr_Oufti}
\alias{Oufti}
\alias{extr_OuftiCSV}
\alias{OuftiCSV}

\title{Function to import Oufti data.
}
\description{
Import oufti meshes, objects and spot localizations directly from .Mat, .txt or .CSV files to R. Note that oufti objects can only be imported from .csv or .txt files.
}
\usage{
extr_Oufti(matfile, mag="No_PixelCorrection", phylo=FALSE, cellList=FALSE)
}

\arguments{
  \item{matfile}{.MAT, .CSV or .TXT of oufti segmentation. Spot data, object data and fluorescence intensity will automatically be detected. Note: object data can only be extracted from .CSV or .TXT files for now!}
  \item{mag}{The pixel-micron conversion factor. See bactMAP::getPixels2um() for more information. When not given, the factor will be 1.}
  \item{phylo}{When TRUE, a list of child-parent nodes & phylo objects containing cell genealogy will be made. Takes very long, so default is set to FALSE}
  \item{cellList}{When TRUE, the cellList (see Value below) will be part of the output list.}
%%     ~~Describe \code{x} here~~
}

\value{
A list containing:
\item{$cellList}{The complete cellList as it is displayed in oufti's output.}
\item{$mesh}{A dataframe containing the following columns: \code{X, Y, cell, frame, num, max.width, length, steplength, max.length, area, Xmid, Ymid, X_rot, Y_rot, angle, Xrot_micron, Yrot_micron, max_um, maxwum}. See standard mesh output for more information.}

When spot data is included:
\item{$spotframe}{A dataframe containing the following columns: \code{l, d, x, y, position, adj_Rsquared, CI_xy, frame, cell}. See \url{www.oufti.org} for more information on the spot data output.}
\item{$spots_relative}{A dataframe containing from the mesh data: max.width, max.length, area, cell, frame. from the spot data: l, d, x, y. Added: spot, totalspot, Lmid, Dum, max_um, maxwum, pole1, pole2.}
\item{$pixel2um}{the pixel to micron conversion factor used}

When a CSV or TXT file containing object data:
\item{$objectframe}{A dataframe containing the coordinates of detected objects, cell ID, frame ID and object ID:  ob_x, ob_y, obnum, obpath, frame, cell, obID}
\item{$object_relative}{A dataframe reminescent of spots_relative containing the following variables: "frame, cell, max.length, obnum, num, obpath, obID, max.width, Dum, Lmid, ob_out_x, ob_out_y, pole1, pole2}

When "phylo" is set to TRUE:
\item{$timelapsedata}{a list of: generation_dataframes: a dataframe with cells with common ancestor, generation_lists: phylo objects which can be used for tree plotting, spot_relative_list: list of spot localization dataframes of cells with common ancestor, and meshdata: a list of cell outline data of cells with common ancestor}
}

\references{
\url{www.oufti.org}

Paintdakhi A, et al. (2016) Oufti: an integrated software package for high-accuracy, high-throughput quantitative microscopy analysis. \emph{Molecular Microbiology} 99(4):767-777.
}
\author{
Renske van Raaphorst
}


\examples{
\dontrun{
##upload matlab oufti file
matfile <- file.choose()
##extract data using extr.Oufti()
output <- extr_Oufti(matfile)
##get spot dataframe
spots <- output$spotframe
##plot spots of frame no.1
ggplot2::ggplot(spots[spots$frame==1,],
    ggplot2::aes(x=x, y=y, color=as.factor(cell))) +
  ggplot2::geom_point() +
  ggplot2::theme_minimal()
}
}

