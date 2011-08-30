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
calculation, e.g., c("(Cd110)D","(Cs111)D"). Strings will be matched against
the parameter names extracted from the FCS file. The default=NULL will
use all parameters. 
}
  \item{arcsinh_cofactor}{
Cofactor used in arcsinh transform \code{asinh(data/arcsinh_cofactor)} of data
}
  \item{comp}{
Apply compensation matrix if present in SPILL or SPILLOVER keywords
}

}
\value{
The name of the written file is returned.
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
	# Load two-parameters sample data included in package
	data_file_path = paste(installed.packages()["spade","LibPath"],"spade","extdata","SimulatedRawData.fcs",sep=.Platform$file.sep)

	output_dir <- tempdir()
	
	# Compute and annotate FCS file with density
	density_file_path <- paste(output_dir,.Platform$file.sep,basename(data_file_path),".density.fcs",sep="")
	SPADE.addDensityToFCS(data_file_path, density_file_path, cols=c("marker1","marker2"))

	# Downsample FCS file based on density
	downsample_file_path <- paste(output_dir,.Platform$file.sep,basename(data_file_path),".density.fcs",sep="")
	SPADE.downsampleFCS(density_file_path, downsample_file_path)

	# Create tree from downsampled FCS file
	cells_file_path <- paste(output_dir,"clusters.fcs",sep="")
	clust_file_path <- paste(output_dir,"clusters.table",sep="")
	graph_file_path <- paste(output_dir,"mst.gml",sep="")
	SPADE.FCSToTree(downsample_file_path, cells_file_path, graph_file_path, clust_file_path, cols=c("marker1","marker2"))

	# Add cluster to FCS files (known as "upsampling")
	upsample_file_path <- paste(density_file_path,"cluster.fcs",sep=".")
	SPADE.addClusterToFCS(density_file_path, upsample_file_path, cells_file_path, cols = c("marker1","marker2"))
}
