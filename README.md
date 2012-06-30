# SPADE -- Spanning Tree Progression of Density Normalized Events

SPADE is a visualization and analysis tool for high-dimensional flow cytometry data. SPADE is implemented as an R package and can be installed via R's packaging facilities. Additionally a GUI is provided, as a [Cytoscape](http://www.cytoscape.org) plugin, for setting-up and interactively visualizing the results of SPADE analyses. Please see the project homepage at [cytospade.org](http://www.cytospade.org) or the github wiki for the primary documentation. This README is primarily targeted at developers working on SPADE itself.

## Prerequisites
1. Latest version of [R](http://www.r-project.org/) and the following packages: igraph0, flowCore

## Setup
The SPADE package has a C++ core that must be built before use. SPADE successfully builds on Linux, OSX and Windows (with Rtools), although Windows users will not be able to take advantage of the OpenMP parallelization used within SPADE. SPADE can be installed from the command line via

    $ R CMD INSTALL <SPADE PATH>


### Building on Windows
You will need to install the [Rtools](http://www.murdoch-sutherland.com/Rtools) that matches your R distribution. (Important note: if you're using R 2.14, you should use Rtools 2.15. Even though Rtools 2.14 is listed as compatible, it requires extra work to make g++ work. Rtools 2.15 works out of the box.) After installation make sure that your `PATH` contains the neccessary Rtools binary directories, e.g.:

1. Open *Control Panel -> System*
1. Click on *Advanced Tab* and then on *Environment Variables*
1. Highlight *PATH* and click *Edit*
1. In the character string in *Variable Value*, make sure the following appear (before other/older compilers):

    c:\Rtools\bin;c:\Rtools\perl\bin;c:\Rtools\gcc-4.6.3\bin;c:\Rtools\MinGW\bin;c:\Program Files\R\<your R version>\bin;


## Building Packages
Source packages can be built with

    $ R CMD build <SPADE PATH>

and binary packages with

    $ R CMD build --no-vignettes --binary <SPADE PATH>

## Tips and Resources
* [R manual for writing extensions](http://cran.r-project.org/doc/manuals/R-exts.html)

## Citations
SPADE was developed in the Plevritis and Nolan Labs at Stanford University, and is described in the following publications:

* Peng Qiu, Erin F. Simonds, Sean C. Bendall, Kenneth D. Gibbs, Robert V. Bruggner, Michael D. Linderman, Karen Sachs, Garry P. Nolan, Sylvia K. Plevritis, ["Phenotypically determined self-organization of flow cytometry data with spanning-tree progression analysis of density normalized events"](http://dx.doi.org/doi%3A10.1038%2Fnbt.1991). *Nature Biotechnolgy*, 2011.
* Sean C. Bendall, Erin F. Simonds, Peng Qiu, El-ad D. Amir, Peter O. Krutzik, Rachel Finck, Robert V. Bruggner, Rachel Melamed, Angelica Trejo, Olga I. Ornatsky, Robert S. Balderas, Sylvia K. Plevritis, Karen Sachs, Dana Pe’er, Scott D. Tanner, Garry P. Nolan, ["Single-Cell mass cytometry of differential immune and drug responses across a human hematopoietic continuum"](http://dx.doi.org/10.1126%2Fscience.1198704), *Science*, 332 (6030): 687-696.
* Michael D. Linderman, Erin F. Simonds, Peng Qiu, Zach Bjornson, Nikesh Kotecha, Teresa H. Meng, Sylvia Plrevritis, Garry P. Nolan, "Algorithmically recovering the hematopoietic lineage from high-dimensional cytometry data using SPADE", *Congress of the Intl. Society for Advancement of Cytometry*, 2010.
