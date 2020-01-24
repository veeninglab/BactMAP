\name{bactmapnews}
\alias{bactmapnews}

\title{What's new in BactMAP 0.1.0.2?}

\description{

1. \strong{Bug fixes with importing SuperSegger data}{problems with importing SuperSegger data through extr_SuperSeggerCells() now fixed.}

2. \strong{Variable name update:} \emph{Xrotum & Yrotum} are now named \emph{Xrot_micron & Yrot_micron}. Functions using these variables have been updated and should recognize both the old and the new variable names. In case something doesn't work with old datasets, please contact us via \email{veeninglab@@gmail.com}. We also added a new function \code{\link{checkVersionCompatible}} to check your data and solve possible compatibility problems.

3. \strong{Installation of ggtree:} Now, ggtree is automatically installed upon installation of BactMAP.

4. \strong{Leaner functions:} Import functions & plotting functions are now set to return only the most commonly used datasets, unless asked for more, to keep more working memory free. In practice, this means that: import functions only give back the dataframe \code{cellList} when asked for by setting \code{cellList=TRUE} in the function. Other functions where returning of redundant data is turned off are: \code{\link{combineDataframes}}, \code{\link{createPlotList}}, Check each individual function documentation for more information.

5. \strong{Updated tutorials & example code:} Check \href{https://github.com/veeninglab/bactmap/wiki}{the BactMAP wiki} for updated tutorials, installation instructions & example workflows.
}


