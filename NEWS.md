## News

### Version 0.1.0.3

In this version, the following improvements are implemented:

1. extr_SuperSeggerCells() import bug in Linux and MacOS fixed (version 0.1.0.2)

2. This version contains a few bug-fixes when it comes to importing data from MicrobeJ and Oufti. Check the wiki and/or documentation for more information on the import functions! 

3. It also contains two new functions for analysis after import: orientCells and onePerCell. These functions can be useful to flip cells based on spot or object localization (orientCells) and to summarize data when you're not very familiar with R (onePerCell). Check out the help (by typing ?orientCells or ?onePerCell) to find out how they work.

### Version 0.1.0.1

In this version, the following improvements are implemented:

1. Variable name update: `Xrotum` & `Yrotum` are now named `Xrot_micron` & `Yrot_micron`. Functions using these variables have been updated and should recognize both the old and the new variable names. In case something doesn't work with old datasets, please contact us via [email](mailto::veeninglab@gmail.com). We also added a new function `checkVersionCompatible()` to check your data and solve possible compatibility problems.

2. Installation of *ggtree*: Now, ggtree is automatically installed upon installation of BactMAP.

3. Leaner functions: Import functions & plotting functions are now set to return only the most commonly used datasets, unless asked for more, to keep more working memory free. In practice, this means that: import functions only give back the dataframe cellList when asked for by setting `cellList=TRUE` in the function. Other functions where returning of redundant data is turned off are: combineDataframes, createPlotList, Check each individual function documentation for more information.

4. Updated tutorials & example code: Check the BactMAP wiki for updated tutorials, installation instructions & example workflows.

### Version 0.1.0.0

This is the first official version of BactMAP. For those who have been using the development versions for a while, there are a few changes:

1. Naming consistencies: *createPlotlist* was changed into *createPlotList* ; *plottreebasic* to *plotTreeBasic*.
2. All *extr*-functions now have a **_** (underscore) behind *extr*, instead of a **.** (dot).
3. The amount of dependencies for *BactMAP* grew. If it is not possible to immediately install *ggtree* during installation of *BactMAP*, manually install it via *BiocManager*. I'm working on a solution to make this package not mandatory.
4. Find a full manual, including tutorials, on https://github.com/vrrenske/BactMAP/wiki
5. Please report any issues on https://github.com/vrrenske/BactMAP/issues.

