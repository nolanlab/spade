# SPADE: Spanning Tree Progression of Density Normalized Events

SPADE is a visualization and analysis tool for high-dimensional flow cytometry data. SPADE is implemented as an R package that takes in FCS files and settings, and outputs graph files (GML), PDFs of the graphs and FCS files with the "cluster" column appended. Please see the project homepage at [cytospade.org](http://www.cytospade.org) or the [wiki](https://github.com/nolanlab/spade/wiki) for the primary documentation.

We no longer support the Cytoscape GUI. An interactive viewer is available from www.cytobank.org. The PDF files can be edited (e.g. graphs can be rearranged) with software such as Adobe Illustrator. For more complex usages, many languages are able to read GML files, including R, Mathematica and Python.

## User Setup

We no longer support distribution via CRAN or Bioconductor. The instructions below install directly from Github and will ensure you have the latest release.

1. If you do not have version 3.0 or later of [R](http://www.r-project.org/), install it ([OSX](http://cran.rstudio.com/bin/macosx/), [Windows](http://cran.rstudio.com/bin/windows/), [Linux](http://cran.rstudio.com/bin/linux/)). Unless you have a compelling reason to do otherwise, we suggest the 64-bit version of R. 

1. Install the devtools, Rclusterpp (latest) and SPADE packages:

        R> install.packages("devtools")
        R> library(devtools)
        R> devtools::install_github("nolanlab/Rclusterpp")
        R> devtools::install_github("nolanlab/spade")

(The version of Rclusterpp in CRAN (0.2.3) is out of date; the above installs the latest (i.e. >= 0.2.4).)

1. Test the installation by typing at the R prompt:

        R> library(spade)

    You should see output like the following:

        R> library('spade')
        Loading required package: igraph
        Loading required package: Rclusterpp
        Loading required package: Rcpp
        Loading required package: RcppEigen

## Usage
Check out the [example usage wiki page](https://github.com/nolanlab/spade/wiki/Example-Usage). Once you are familiar with the workflow, you can edit the file R/inst/runSPADE.R to setup your analysis, then run the file. For additional documentation about the R package, you can view the package vignette with `vignette("SPADE")` at the R prompt. Additionally all of the functions in the SPADE R package are documented; view their manual pages with `?<function>`, e.g., `?SPADE.driver`, at the R prompt.

## Developer Setup

Please refer to the [wiki pages](https://github.com/nolanlab/spade/wiki).

- - -

## Citations
SPADE was developed in the Plevritis and Nolan Labs at Stanford University, and is described in the following publications:

* Linderman MD, Bjornson Z, Simonds EF, Qiu P, Bruggner RV, Sheode K, Meng TH, Plevritis SK, Nolan GP, ["CytoSPADE: high-performance analysis and visualization of high-dimensional cytometry data"](http://www.ncbi.nlm.nih.gov/pubmed/22782546). *Bioinformatics*, 2012.
* Peng Qiu, Erin F. Simonds, Sean C. Bendall, Kenneth D. Gibbs, Robert V. Bruggner, Michael D. Linderman, Karen Sachs, Garry P. Nolan, Sylvia K. Plevritis, ["Phenotypically determined self-organization of flow cytometry data with spanning-tree progression analysis of density normalized events"](http://dx.doi.org/doi%3A10.1038%2Fnbt.1991). *Nature Biotechnolgy*, 2011.
* Sean C. Bendall, Erin F. Simonds, Peng Qiu, El-ad D. Amir, Peter O. Krutzik, Rachel Finck, Robert V. Bruggner, Rachel Melamed, Angelica Trejo, Olga I. Ornatsky, Robert S. Balderas, Sylvia K. Plevritis, Karen Sachs, Dana Pe’er, Scott D. Tanner, Garry P. Nolan, ["Single-Cell mass cytometry of differential immune and drug responses across a human hematopoietic continuum"](http://dx.doi.org/10.1126%2Fscience.1198704), *Science*, 332 (6030): 687-696.
* Michael D. Linderman, Erin F. Simonds, Peng Qiu, Zach Bjornson, Nikesh Kotecha, Teresa H. Meng, Sylvia Plrevritis, Garry P. Nolan, "Algorithmically recovering the hematopoietic lineage from high-dimensional cytometry data using SPADE", *Congress of the Intl. Society for Advancement of Cytometry*, 2010.
