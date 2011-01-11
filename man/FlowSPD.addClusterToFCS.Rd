\name{SPADE.addClusterToFCS}
\alias{SPADE.addClusterToFCS}
\title{
Annotate observations in FCS file with cluster assignment
}
\description{
Annotate observations in a FCS file with cluster assignment
}
\usage{
SPADE.addClusterToFCS(infilename, outfilename, clusterfilename, cols = NULL, arcsinh_cofactor = 5, comp=TRUE)
}
\arguments{
  \item{infilename}{
Name of input FCS file
}
  \item{outfilename}{
Name of output FCS file
}
  \item{clusterfilename}{
Name of FCS file with subset of cells used in clustering
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
  \item{comp}{
Apply compensation matrix if present in SPILL or SPILLOVER keywords
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
\code{\link{SPADE.FCSToTree}}
}
\examples{
}
