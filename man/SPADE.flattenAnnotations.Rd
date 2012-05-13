\name{SPADE.flattenAnnotations}
\alias{SPADE.flattenAnnotations}
\title{
Flatten list of annotations to matrix
}
\description{
Helper function for flattening list of annotations 
}
\usage{
SPADE.flattenAnnotations(annotations)
}

\arguments{
  \item{annotations}{
A list of annotation matrices. All matrices must have the same number of rows.
}
}
\value{
Single matrix of annotations
}
\author{
Michael Linderman
}
\examples{
	# Not run
	## Load two-parameters sample data included in package
	#data_file_path = paste(installed.packages()["spade","LibPath"],"spade","extdata","SimulatedRawData.fcs",sep=.Platform$file.sep)

	## Run basic SPADE analyses, clustering on two parameters. 
	#output_dir <- tempdir()
	#SPADE.driver(data_file_path, out_dir=output_dir, cluster_cols=c("marker1","marker2"))

	## Compute medians, counts and other parameters from processed files
	#upsampled_file_path <- paste(output_dir,.Platform$file.sep,basename(data_file_path),".density.fcs.cluster.fcs",sep="")
	#mst_graph <- igraph0:::read.graph(paste(output_dir,"mst.gml",sep=.Platform$file.sep),format="gml")
	#anno <- SPADE.markerMedians(upsampled_file_path, igraph0:::vcount(mst_graph), cols = c("marker1","marker2"))

	## Flatten annotations so they can easily be saved in table
	#flat_anno <- SPADE.flattenAnnotations(anno)
}
