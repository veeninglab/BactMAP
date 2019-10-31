\name{extr_Morphometrics}
\alias{extr_Morphometrics}
\alias{Morphometrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Upload cell outlines segmented using Morphometrics
}
\description{
The 'CONTOURS'-output from Morphometrics is uploaded and converted to a *mesh*-dataframe compatible with bactMAPs plotting functions. Find more examples of usage of \code{extr_Morphometrics} on the BactMAP \href{https://github.com/vrrenske/bactMAP/wiki}{wiki}.
}
\usage{
extr_Morphometrics(morphpath, mag, turncells = TRUE, cellList=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{morphpath}{
  the path to the 'CONTOURS'- output of Morphometrics (for instance "mycells_01-01-2019_CONTOURS.mat").
}
  \item{mag}{
  pixel conversion factor (see \link{pixels2um}). if added, the coordinates and cell length/width will be converted to micron
}
  \item{turncells}{
  if set to TRUE (default), X and Y coordinates will be turned so that mid-cell is at [0,0], while the length axis of the cell's bounding box is parallel to the x axis.
}
  \item{cellList}{
  when cellList==TRUE, the value \code{cellList} (see below) is part of the output.
  }
}

\value{

  \item{cellList }{dataframe containing the original data extracted from Morphometrics.}
  \item{mesh }{dataframe containing the coordinates of the outlines of the cell, cell dimensions and if *turncells==TRUE*, the coordinates of each cell turned and repositioned to [0,0]}
  \item{pixel2um}{if *mag* was defined in the arguments, the conversion factor is saved in the output list.}
%% ...
}
\references{
 Ursell T, Lee TK, Shiomi D, et al. Rapid, precise quantification of bacterial cellular dimensions across a genomic-scale knockout library. BMC Biol. 2017;15(1):17. Published 2017 Feb 21. doi:10.1186/s12915-017-0348-8
}
\author{
Renske van Raaphorst
}

\examples{

\dontrun{
#Upload Morphometrics output from bactMAP example files

Example_MM <-  file.choose()
Example_output <- extr_Morphometrics(Example_MM)$mesh

#check output:
summary(Example_output)

#plot cells
plotRaw(Example_output, frameN=1)
}
}

