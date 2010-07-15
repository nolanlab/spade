\name{FlowSPD.addClusterToFCS}
\alias{FlowSPD.addClusterToFCS}
\title{
Annotate observations in FCS file with cluster assignment
}
\description{
Annotate observations in a FCS file with cluster assignment
}
\usage{
FlowSPD.addClusterToFCS(infilename, outfilename, clusterfilename, cols = NULL, arcsinh_cofactor = 5)
}
\arguments{
  \item{infilename}{
Name of input FCS file
}
  \item{outfilename}{
Name of output FCS file
}
  \item{clusterfilename}{
Name of file containing table of cluster centers with header line
}
  \item{cols}{
Usually a vector of strings specifying the columns to be used in the density
calculation, e.g., c("115-CD45","110-CD3"). Strings will be matched against
the parameter descriptions extracted from the FCS file. The default=NULL will
use all parameters. 
}
  \item{arcsinh_cofactor}{
Cofactor used in arcsinh transform \code{asinh(data/arcsinh_cofactor)} of data
}
}
\details{
}
\value{
The name of the written file is returned.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Michael Linderman
}
\note{
Underlying implementations have been parallelized with OpenMP. Set
\code{OMP_NUM_THREADS} in environment to control the number of threads used.
}
\seealso{
\code{\link{FlowSPD.FCSToTree}}
}
\examples{
}
