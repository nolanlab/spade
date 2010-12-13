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
SPADE.markerMedians(files, cols = NULL, arcsinh_cofactor = 5, cluster_cols=NULL)
SPADE.annotateMarkers(files, cols = NULL, arcsinh_cofactor = 5)
}
\arguments{
  \item{files}{
Name of input FCS file or vector of input FCS file names. FCS files must have "cluster" column.
}
  \item{cols}{
Usually a vector of strings specifying the columns to be used in the density calculation, e.g., c("115-CD45","110-CD3"). Strings will be matched against the parameter descriptions extracted from the FCS file. The default=NULL will use all parameters.
}
  \item{arcsinh_cofactor}{
Cofactor used in the arcsinh transform \code{asinh(data/arcsinh_cofactor)} of data
}
  \item{cluster_cols}{
A vector of strings specifying columns that should be marked as having been used in clustering
}
}
\details{
SPADE.annotateMarkers is deprecated.
}
\value{
List with:
\item{count}{Matrix of observation count for clusters}
\item{medians}{Matrix of medians for specified columns}
}
\references{
}
\author{
Michael Linderman
}
\note{
}
\seealso{
\code{\link{SPADE.addClusterToFCS}}, \code{\link{SPADE.annotateGraph}}
}
\examples{
}
