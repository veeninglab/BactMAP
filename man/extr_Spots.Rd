\name{extr_Spots}
\alias{extr_Spots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generic function for extraction of spot data
}
\description{
To upload spot localization data to R in the standard BactMAP format, save your spot data as a .txt or .csv file with (at least) the following columns:

\code{x} - the x coordinate of each spot

\code{y} - the y coordinate of each spot

\code{frame} - the number of the image frame

The function also recognises the following columns:

\code{trajectory} - ID variable for a spot trajectory

\code{trajectory_length} - the recorded length of each trajectory

\code{displacement_sq} - the squared displacement of each spot

Other columns will be saved in the output "cellList" (see below).
}

\usage{
extr_Spots(dataloc, seperator = ",", cellList=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataloc}{
  The path to the .txt or .csv file containing the spot data.
}
  \item{seperator}{
 Indicates the column seperator of the .txt or .csv. Default = ","
}
  \item{cellList}{
  When TRUE, the value 'cellList' is returned (see below).
  }
}

\value{
\item{cellList}{the original input data as dataframe}
\item{spotframe}{dataframe containing the minimal required columns (x, y & frame) and if given, the extra columns (trajectory, trajectory_length & displacement_sq)}
}

\author{
Renske van Raaphorst
}

\examples{
\dontrun{
spotpath <- file.choose()

spot_output <- extr_Spots(spotpath, seperator="/t")$spotframe
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
