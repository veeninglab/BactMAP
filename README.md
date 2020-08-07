BacTMAP
================

## BacTMAP - Bacteria Tool for Microscopy Analysis & Plotting

<img src="man//figures//logo_bactmap_pink.png" width="200" height="200" />

This package is meant to make it easier for microbiologists to combine
and analyse segmentation & fluorescence data derived from custom
software like Oufti, Morphometrics, MicrobeJ or Supersegger. There are
standard functions for importing several popular programs, but also
options for importing CSVs or Matlab files with a specific structure. 

BactMAP uses ggplot and many ggplot-extensions to make visual summaries of the
segmentation & fluorescence data, show the cell genealogy, plot protein
trajectories in the cell and many other things. This can give you a
quick overview of your data to make a decision on further analysis and
custom visualization.


If you have any comments, requests or find bugs, please report it on the [Issues page](https://github.com/veeninglab/bactmap/issues).**
=======
## Manual & Tutorials

In the [wiki](https://github.com/vrrenske/BactMAP/wiki) you will find a
manual and some examples on how to use BactMAP as a gateway between your
segmentation & fluorescence data.

## Download and install package

    #install devtools if not done yet
    install.packages("devtools")
    
    #install bactmap from my github repository
    devtools::install_github("veeninglab/bactMAP")
    
    #or from the stable build:
    devtools::install_github("veeninglab/bactMAP")
    
    #or install it with all dependencies
    devtools::install_github("veeninglab/bactMAP", dependencies=TRUE)
    
    #load package
    library(bactMAP)

## Development version

Find the development version [here](https://github.com/vrrenske/bactMAP).

## Version 0.1.0.2

This new version contains bug fixes of `extr_MicrobeJ` and `extr_Oufti`. It also has a few new functions: `orientCells` and `onePerCell`. The first function helps you flip cells based on non-symmetrical fluorescence localization (for instance, a polar fluorescent spot) - the latter is useful if you want to make histograms of cell lengths, etc. based on `mesh` or `object` files. A short description of the functions can be found in the manual inside the package (type `?orientCells` to find the help of this function, for instance). I will add a small tutorial and documentation to the wiki shortly.

