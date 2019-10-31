\name{extr_ISBatch}
\alias{extr_ISBatch}
\alias{peakfitter}
\alias{ISBatch}

\title{
Import ISBatch/Peakfitter data into R
}
\description{
Import ISBatch/Peakfitter spot coordinates and/or trajectories from a .txt or .csv into R.
}
\usage{
extr_ISBatch(dataloc, seperator=",", cellList=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataloc}{
  Path to the results table from the "peakfitter" option of the ISBatch imageJ plugin saved as .CSV or tab-delimited .TXT.
}
  \item{seperator}{
  When .CSV is not comma-separated but by something else, indicate this here.
  }
    \item{cellList}{
  when TRUE, the value cellList is part of the returned output (see below).
  }
}

\value{
  \item{cellList}{
  copy of original data as dataframe
  }
  \item{spotframe}{
  when data is made using ISBatch's peak fitter only, this data frame contains only x, y -coordinates and frame numbers. When trajectories are recorded, it also contains the squared-displacement, trajectory ID and trajectory length.
  }

}
\references{
\url{http://singlemolecule.github.io/iSBatch/}

Caldas, V.E., Punter, C.M., Ghodke, H., Robinson, A. and van Oijen, A.M., 2015. iSBatch: a batch-processing platform for data analysis and exploration of live-cell single-molecule microscopy images and other hierarchical datasets. \emph{Molecular BioSystems}, 11(10), pp.2699-2708.
}
\author{
Renske van Raaphorst
}

\seealso{
Use \href{https://github.com/vrrenske/BactMAP/wiki/spotsInBox}{spotsInBox} to combine the spot data with segmentation data from another program.
}
\examples{
\dontrun{
##Choose ISBatch peakfitter output file
dataloc <- file.choose()
##Convert to dataframe
spotframe <- extr_ISBatch(dataloc)
}
}

