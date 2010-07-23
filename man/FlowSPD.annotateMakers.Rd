\name{FlowSPD.annotateMarkers}
\alias{FlowSPD.annotateMarkers}
\title{
Compute marker medians and counts for clusters
}
\description{
Compute the marker medians and observations counts for cluster annoated FCS files.
}
\usage{
FlowSPD.annotateMarkers(infilename, cols = NULL, arcsinh_cofactor = 5)
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
}
\value{
List with:
\item{counts}{Matrix of observation count for clusters}
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
\code{\link{FlowSPD.assignToCluster}}}
\examples{
}
