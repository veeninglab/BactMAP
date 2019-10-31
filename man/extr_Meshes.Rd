\name{extr_Meshes}
\alias{extr_Meshes}

\title{
Generic function to upload cell outlines from a .csv or .txt file.
}
\description{
Output needs to have at least the following columns:

\code{X} - The x-coordinate of one cell outline point.

\code{Y} - The y-coordinate of one cell outline point.

\code{frame} - The image frame number.

\code{cell} - The cell identifyer (numeric).

Optional:

\code{max.length} - The cell length (in pixels)

\code{max.width} - The cell diameter (in pixels)

Any other parameters will be saved in the output dataframe cellList.
}
\usage{
extr_Meshes(dataloc, sep = ",", turn = TRUE, mag, cellList=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataloc}{
The path to the .csv or .txt - file containing the mesh data.
}
  \item{sep}{
The separator of the columns in the .txt or .csv file, default = ","
}
  \item{turn}{
When \code{turn==TRUE}, the cells will be turned so that the mid-point is on coordinate [0,0] and the length axis parallel to the x axis.
}
  \item{mag}{
The magnification convertion factor (see addPixels2um). When given, the output will be also given in microns.
}
  \item{cellList}{
When TRUE, cellList will be part of the output of this function (see Value below).
  }
}

\value{
\item{cellList}{the original input as dataframe}
\item{mesh}{the mesh data compatible with BactMAP}
}

\author{
Renske van Raaphorst
}

\seealso{
\code{\link{pixel2um}}, \code{\link{mesh}}
}
\examples{
\dontrun{
meshdata <- file.choose()

mesh_BM <- extr_Meshes(meshdata, sep=",", mag="100x_FRAP", turn=FALSE)
}
}
