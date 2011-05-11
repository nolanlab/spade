\name{SPADE.markerMedians}
\alias{SPADE.markerMedians}
\alias{SPADE.annotateMarkers}
\title{
Compute marker medians and counts for clusters
}
\description{
Compute the marker medians and observations counts for cluster annoated FCS files.
}
\usage{
SPADE.markerMedians(files, num.clusters, cols = NULL, arcsinh_cofactor = 5, cluster_cols=NULL, comp=TRUE)
SPADE.annotateMarkers(files, cols = NULL, arcsinh_cofactor = 5)
}
\arguments{
  \item{files}{
Name of input FCS file or vector of input FCS file names. FCS files must have "cluster" column.
}
  \item{num.clusters}{
Number of clusters. Note not all clusters need to be present in all files.
  }
  \item{cols}{
Usually a vector of strings specifying the columns to be used in the density
calculation, e.g., c("(Cd110)D","(Cs111)D"). Strings will be matched against
the parameter names extracted from the FCS file. The default=NULL will
use all parameters. 
}
  \item{arcsinh_cofactor}{
Cofactor used in the arcsinh transform \code{asinh(data/arcsinh_cofactor)} of data
}
  \item{cluster_cols}{
A vector of strings specifying columns that should be marked as having been used in clustering
}
  \item{comp}{
Apply compensation matrix if present in SPILL or SPILLOVER keywords
}
}
\details{
SPADE.annotateMarkers is deprecated.
}
\value{
List with:
\item{count}{Matrix of observation count for clusters}
\item{percenttotal}{Matrix of percent of total number of cells [0-100] in each cluster}
\item{medians}{Matrix of medians for specified columns}
}
\author{
Michael Linderman
}
\seealso{
\code{\link{SPADE.addClusterToFCS}}, \code{\link{SPADE.annotateGraph}}
}
