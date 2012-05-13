\name{SPADE.write.graph}
\alias{SPADE.write.graph}
\title{
Writing the graph to a file in some format
}
\description{
General function for exporting graphs to foreign file formats, however at present only the GML format is implemented.
}
\usage{
SPADE.write.graph(graph, file = "", format = c("gml"))
}
\arguments{
  \item{graph}{
The graph to export
}
  \item{file}{
A connection or a string giving the file name to write the graph to.
}
  \item{format}{
Character string giving the file format.
}
}
\details{
GML is general textual format for graphs.

The vertex and edge attributes are written to the file if they are numeric or strings. Currently only the graphics struct is supported, and only for vertices; \code{graphics.x} indicates an \code{x} attribute in the \code{graphics} struct.
}
\value{
A NULL, invisibly
}
\author{
Michael Linderman
}
\seealso{
\code{\link{write.graph}}
}
\examples{
	# Not run
	## Load two-parameters sample data included in package
	#data_file_path = paste(installed.packages()["spade","LibPath"],"spade","extdata","SimulatedRawData.fcs",sep=.Platform$file.sep)

	## Run basic SPADE analyses, clustering on two parameters. 
	#output_dir <- tempdir()
	#SPADE.driver(data_file_path, out_dir=output_dir, cluster_cols=c("marker1","marker2"))

	## Read and write minimum spanning tree graph 
	#mst_graph <- igraph0:::read.graph(paste(output_dir,"mst.gml",sep=.Platform$file.sep),format="gml")
	#SPADE.write.graph(mst_graph, file = paste(output_dir,"new_mst.gml",sep=.Platform$file.sep), format = c("gml"))
}
