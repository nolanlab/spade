\name{SPADE.FCSToTree}
\alias{SPADE.FCSToTree}
\title{
Cluster and build minimum spanning tree from data in FCS files 
}
\description{
Hierarchically cluster observations in a set of FCS files and build a minimum
spanning tree connecting those clusters.
}
\usage{
SPADE.FCSToTree(infilenames, outfilename, graphfilename, clusterfilename,
		  cols = NULL, k = 200, arcsinh_cofactor = 5, 
		  desired_samples = 50000, comp=TRUE)
}
\arguments{
  \item{infilenames}{
Vector of FCS file names that should be used as input
}
  \item{outfilename}{
Name of FCS file to write subset of cells used for clustering along with their cluster assignment
}
  \item{graphfilename}{
Name of file to write gml graph description
}
  \item{clusterfilename}{
Name of file to write table of cluster centers
}
  \item{cols}{
Usually a vector of strings specifying the columns to be used in the density
calculation, e.g., c("(Cd110)D","(Cs111)D"). Strings will be matched against
the parameter names extracted from the FCS file. The default=NULL will
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
  \item{comp}{
Apply compensation matrix if present in SPILL or SPILLOVER keywords
}

}
\value{
None.
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
\code{\link{SPADE.downsampleFCS}}
}
\examples{
	## Load two-parameters sample data included in package
	#data_file_path = paste(installed.packages()["spade","LibPath"],"spade","extdata","SimulatedRawData.fcs",sep=.Platform$file.sep)

	#output_dir <- tempdir()
	#
	## Compute and annotate FCS file with density
	#density_file_path <- paste(output_dir,.Platform$file.sep,basename(data_file_path),".density.fcs",sep="")
	#SPADE.addDensityToFCS(data_file_path, density_file_path, cols=c("marker1","marker2"))

	## Downsample FCS file based on density
	#downsample_file_path <- paste(output_dir,.Platform$file.sep,basename(data_file_path),".density.fcs",sep="")
	#SPADE.downsampleFCS(density_file_path, downsample_file_path)

	## Create tree from downsampled FCS file
	#cells_file_path <- paste(output_dir,"clusters.fcs",sep="")
	#clust_file_path <- paste(output_dir,"clusters.table",sep="")
	#graph_file_path <- paste(output_dir,"mst.gml",sep="")
	#SPADE.FCSToTree(downsample_file_path, cells_file_path, graph_file_path, clust_file_path, cols=c("marker1","marker2"))
}
