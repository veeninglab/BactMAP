\name{createPlotList}
\alias{createPlotList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to create descriptive plots of spot localization per group of cells
}
\description{
\code{createPlotList()} is a plotting function which combines detected fluorescent spot localization with mesh outline information to plot demographs, projections, histograms and a data summary in one go. These functions are all combined because they need similar data preparation. The idea is that after using this function, you have a good overview of the localization of your protein inside the cells. Sometimes this analysis can be enough to show where a protein localizes, but in other cases (when proteins are very mobile or localize heterogeneously, for instance), this may rather serve as a starting point for more detailed analysis. Therefore, the data generated to make these plots is also returned by the function.
}
\usage{
createPlotList(spotdata,
               meshdata,
               groups=4,
               colorpalette="GreenYellow",
               mag="No_PixelCorrection",
               AllPlot=F,
               Xm="X",
               Ym="Y",
               viridis=FALSE,
               showPlot=TRUE,
               getData=FALSE,
               getSummary=TRUE)
}
\arguments{
  \item{spotdata}{Data frame containing -at least- column with spot x-localizations, y-localizations and frame number. When the variable meshdata is not given, additional required columns are cell ID, relative x-localization ("l" or "Lmid") and relative y-localization ("d" or "Dum").}
  \item{meshdata}{Data frame containing -at least- the x and y coordinates of the cell meshes, the maximum cell width, maximum cell length, cell ID and frame. Optional.}
  \item{groups}{Number of groups the data will be divided in. Default = 4.}
  \item{colorpalette}{Color palette the heatmap plots will be plotted in. Default is "GreenYellow". Check which palettes are loaded by \code{showColorPaletes()} and add your own palette using \code{addPallete()}}
  \item{mag}{The name of the conversion factor for pixel to micrometer. Default is "No_PixelCorrection" which is 1. Check \code{\link{getPixels2um}} for the loaded conversionfactors and use \code{addPixels2um("ConversionName", conversionfactor)} to add a conversion factor to the list.}
  \item{AllPlot}{When \code{TRUE}, the output x/y plots will have sideplots with histograms of the spot localization. When \code{FALSE}, only the x/y plots will be plotted.}
  \item{Xm}{By default, the column names of the x/y coordinates in the MESH data frames are "X" and "Y". When this is not the case, one can tell the function by using \code{Xm} and \code{Ym}. Enter a character string with the name of the x-coordinate column for Xm, and the name of the y-coordinate column for Ym.}
  \item{Ym}{See Xm}
  \item{viridis}{when \code{viridis==TRUE}, the \code{colorpalette} needs to be one of the \code{viridis} color palettes ("A", "B", "C", "D", "E", "magma", "cividis", "viridis", "plasma" or "inferno")}
  \item{showPlot}{when \code{showPlot==TRUE}, a menu will be shown after the function is finished, which enables the user to view the plot output of the function.}
  \item{getData}{when \code{getData==TRUE}, the output list will include the input dataframes \code{spotdata} and (if given) \code{meshdata}, which include the quantile classification and the relative localization of the spots.}
  \item{getSummary}{when \code{getSummary==TRUE}, the output list will include the \code{data_summary} (see below).}
}

\value{
A list consisting of the following items:
\item{$lengthplot}{Heatmap histogram plot of the length-axis localizations of the spots ordered by the length of the cell.}
\item{$widthplot}{Heatmap histogram plot of the width-axis localizations of the spots ordered by the length of the cell.}
\item{$qplots_separate}{If the outlines of the cell are not added to the function (dataframe \code{MESH} not indicated), the output will be a list of plots of length \code{inp}, with relative x/y coordinates of the spots plotted as a density function. The data is grouped by the variable \code{bycol}, which is \code{length} by default.If the outlines of the cells are added to the function (dataframe \code{MESH} indicated), the average cell outline of each group will be overlayed on top of each plot. If \code{AllPlot==TRUE}, histograms of the x and y coordinates of the spots will be added as sideplots.}
\item{$qplots}{for easy use, qplots is a \code{grob} arrangement of the group-plots \code{qplots_separate}: using \code{plot()} you can immediately view them in the Rstudio plots window.}
\item{$plottotal}{Relative x/y coordinates of all spots with average cell outline overlayed. If \code{AllPlot==TRUE}, the x/y coordinates will be plotted as histograms on the x and y axis of the plot.}
\item{$histograms}{1D density plots of the spot localization on the length axis of the cells, grouped by \code{bycol} in the amount of groups defined by \code{groups}.}
\item{$spotdata}{data frame containing the groups, relative x/y localizations of the spots and cell length/width}
\item{$meshdata}{if \code{MESH} is included, the function returns a dataframe including the rotated meshes and groups.}
\item{$data_summary}{a character vector, or when \code{groups>1}, a list of character vectors summarizing the most important data features of each size group and the full dataset.}
}

\references{
These plotting functions were developed for:

\href{https://www.pnas.org/content/114/29/E5959}{van Raaphorst, Renske, Morten Kjos, and Jan-Willem Veening. "Chromosome segregation drives division site selection in Streptococcus pneumoniae." Proceedings of the National Academy of Sciences 114.29 (2017): E5959-E5968.}
}
\author{
Renske van Raaphorst
}
\note{
This function is meant to give an overview of the data. Check the tutorials on the \href{https://github.com/veeninglab/BactMAP/wiki/tutorials}{wiki page} for examples of modifying individual plots.
}


\examples{

\dontrun{
##use data from example dataset

##download example datasets DnaX_tracks.rda & DnaX_mesh.rda
##from https://veeninglab.com/f/example_datasets.zip

load("DnaX_tracks.rda")
load("DnaX_mesh.rda")

##create list of plots, ordered by length, grouped in 6 groups.
GFPlist <- createPlotList(DnaX_tracks, DnaX_mesh, groups=6, colorpalette="GreenYellow", AllPlot=T)

#have a look at the "plottotal" plot
plot(GFPlist$plottotal)

#save PDF of the 6 grouped plots (qplots)
ggplot2::ggsave(GFPlist$qplots, filename="GroupPlots.PDF")

#see summary:
GFPlist$data_summary
}
}
