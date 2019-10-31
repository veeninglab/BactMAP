\name{checkVersionCompatible}

\alias{checkVersionCompatible}

\title{check whether your data frames are compatible with the current version of BactMAP}

\description{In the new version of BactMAP, the variables Xrotum & Yrotum are renamed to Xrot_micron and Yrot_micron. While plotting functions are adapted to both old & new variable names, it is recommended to update your old dataframes so the variable names in all your datasets are the same. Run checkVersionCompatible() to check for compatibility with the new BactMAP version & change the old variable names to the new ones.}

\usage{
checkVersionCompatible(oldDataFrame, returnDataFrame=TRUE)
}

\arguments{
  \item{oldDataFrame}{Any dataframe created with BactMAP of which you wish to know whether it's compatible with the current version of bactMAP.
}
  \item{returnDataFrame}{Set this variable to FALSE if you just want to check for issues, and to TRUE if you want to get a fixed (or, if no problems are found, untouched) dataframe returned.}
}
\details{
The function standardly fixes found compatibility problems (between version 0.1.0.0 and 0.1.0.1, see \link{bactmapnews} ) and returns the fixed dataframe (or the same dataframe if no problems are found). If you just want to check for issues, set returnDataFrame to FALSE.
}
\value{
If returnDataFrame is set to \code{TRUE}, a dataframe with the correct variable names is returned.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Renske van Raaphorst
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
#using TurnedCell4 dataframe from the example dataset (\url{https://veeninglab.com/f/example_datasets.zip})


checkVersionCompatible(TurnedCell4, returnDataFrame=FALSE)
}
}
