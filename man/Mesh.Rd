\name{mesh}
\alias{mesh}
\docType{data}
\title{
Format of the Mesh dataframe
}
\description{
When extracting segmentation data, the main output is a dataframe called 'mesh'. The standard structure of this dataframe is shown below:
}

\format{
  A data frame with at least the following 7 variables:

  \describe{
    \item{\code{X}}{a numeric vector indicating the X-coordinate of one mesh outline point.}
    \item{\code{Y}}{a numeric vector indicating the Y-coordinate of one mesh outline point.}
    \item{\code{cell}}{a numeric vector giving the cell ID (unique per image frame)}
    \item{\code{frame}}{a numeric vector giving the ID of the image frame}
    \item{\code{max.length}}{a numeric vector giving the cell length (in pixels)}
    \item{\code{max.width}}{a numeric vector giving the cell diameter (in pixels)}
    \item{\code{num}}{a numeric vector indicating the order of the cell outline points.}
  }

  When the cells are turned (default) such, that the mid-point is at [0,0] and the length axis is parallel to the x axis, the following 5 variables are added:
  \describe{
    \item{\code{angle}}{the angle of the cell length axis towards the 0 x-axis}
    \item{\code{Xmid}}{the X-coordinate of the cell's midpoint}
    \item{\code{Ymid}}{the Y-coordinate of the cell's midpoint}
    \item{\code{X_rot}}{the X-coordinate of the turned cell}
    \item{\code{Y_rot}}{the Y-coordinate of the turned cell}
  }

  When a magnification conversion factor is indicated, the coordinates and lengths are calculated in micron. The following 4 variables are added:
  \describe{
    \item{\code{Xrot_micron}}{the X-coordinate of the turned cell in micron}
    \item{\code{Yrot_micron}}{the Y-coordinate of the turned cell in micron}
    \item{\code{max_um}}{the cell length in micron}
    \item{\code{maxwum}}{the cell width in micron}
    }

  In most cases, the cell area is also indicated:
  \describe{
    \item{\code{area}}{The cell area in pixel^2}
  }

  Finally, depending on the input programs, some other variables can be added to the output. See the original programs output documentation for clarification.
}


