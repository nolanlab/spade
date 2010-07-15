\name{FlowSPD.FCSToTree}
\alias{FlowSPD.FCSToTree}
\title{
Cluster and build minimum spanning tree from data in FCS files 
}
\description{
Hierarchically cluster observations in a set of FCS files and build a minimum
spanning tree connecting those clusters.
}
\usage{
FlowSPD.FCSToTree(infilenames, graphfilename, clusterfilename,
		  cols = NULL, k = 200, arcsinh_cofactor = 5, 
		  desired_samples = 40000)
}
\arguments{
  \item{infilenames}{
Vector of FCS file names that should be used as input
}
  \item{graphfilename}{
Name of file to write gml graph description
}
  \item{clusterfilename}{
Name of file to write table of cluster centers
}
  \item{cols}{
Usually a vector of strings specifying the columns to be used in the density
calculation, e.g., c("115-CD45","110-CD3"). Strings will be matched against
the parameter descriptions extracted from the FCS file. The default=NULL will
use all parameters. 
}
  \item{k}{
Desired number of clusters. Algorithm might create between [k/2,3k/2]
clusters.
}
  \item{arcsinh_cofactor}{
Cofactor used in arcsinh transform \code{asinh(data/arcsinh_cofactor)} of data
}
  \item{desired_samples}{
Desired number of samples to be used in clustering. Usually leave at default.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
None.
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
Implementation can be very memory intensive.
}
\seealso{
\code{\link{FlowSPD.downsampleFCS}}
}
\examples{
}
