\name{combineDataframes}
\alias{combineDataframes}

\title{
Function to combine data by condition (e.g. mutant vs wild type) and/or channel (e.g. GFP/RFP), meant for BactMAP imported data.
}
\description{
The function adds a column "condition" and/or "channel" to each dataframe in the given *listofdataframes*, after which it concatenates all dataframes. If not all columns in the dataframes match, it discards the columns which don't match. The function returns the concatenated dataframe "finalframe" and the original list of dataframes with the "condition" and/or "channel" column added to each dataframe.
}
\usage{
combineDataframes(listofdataframes, listofconditions, listofchannels, output="finalframe")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listofdataframes}{
list of dataframes of which the common columns need to be combined
}
  \item{listofconditions}{
list of conditions corresponding to each dataframe. Needs to have the same order and length as listofdataframes.
}
  \item{listofchannels}{
list of channels corresponding to each dataframe. Needs to have the same order and length as listofdataframes.
}
  \item{output}{
when output = "all", both 'finalframe' and 'originaldata' will be in the output list. when output = "finalframe" or "originaldata", only the given data will be in the output list.
  }
}

\value{
\item{finalframe}{concatenated dataframe with columns "condition" and/or "channel" added}
\item{originaldata}{the input list of dataframes with columns "condition" and/or "channel" added}
}

\author{
Renske van Raaphorst
}

\examples{
\dontrun{

##download example data from https://veeninglab.com/f/example_data.zip:

##VanFL_objDM.rda
##VanFL_objWT.rda

load("VanFL_objDM.rda")
load("VanFL_objWT.rda")

##find the data documentation on https://veeninglab.com/bactmap

combined <- combineDataframes(list(VanFL_objDM, VanFL_objWT),
                              listofconditions = list("mutant", "WT""))
}
}
