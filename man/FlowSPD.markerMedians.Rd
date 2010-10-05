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
SPADE.markerMedians(infilename, cols = NULL, arcsinh_cofactor = 5)
SPADE.annotateMarkers(infilename, cols = NULL, arcsinh_cofactor = 5)
}
\arguments{
  \item{infilename}{
Name of the input FCS file. Must have "cluster" column.
}
  \item{cols}{
Usually a vector of strings specifying the columns to be used in the density calculation, e.g., c("115-CD45","110-CD3"). Strings will be matched against the parameter descriptions extracted from the FCS file. The default=NULL will use all parameters.
}
  \item{arcsinh_cofactor}{
Cofactor used in the arcsinh transform \code{asinh(data/arcsinh_cofactor)} of data
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
