\name{SPADE.layout.arch}
\alias{SPADE.layout.arch}
\title{
Generate coordinates for plotting graphs
}
\description{
Performing "arch" layouts of graph vertices
}
\usage{
SPADE.layout.arch(mst_graph)
}
\arguments{
  \item{mst_graph}{
The graph to layout. Must be acyclic and undirected.
}
}
\details{
These functions calculate the coordinates of the vertices for a graph.

\code{layout.arch} lays out the longest chain of the graph, the "backbone", on an arch, and the "side chains" as trees normal to that backbone.
}
\value{
All these functions return a numeric matrix with at least two columns, x and y positions, and the same number of lines as the number of vertices.
}
\author{
Michael Linderman
}
\seealso{
\code{\link{SPADE.annotateGraph}}
}
\examples{
	# Load two-parameters sample data included in package
	data_file_path = paste(installed.packages()["spade","LibPath"],"spade","extdata","SimulatedRawData.fcs",sep=.Platform$file.sep)

	# Run basic SPADE analyses, clustering on two parameters. 
	output_dir <- tempdir()
	SPADE.driver(data_file_path, out_dir=output_dir, cluster_cols=c("marker1","marker2"))

	# Generate PDFs of annotated graphs (into output_dir/pdf) using arch layout
	mst_graph <- igraph:::read.graph(paste(output_dir,"mst.gml",sep=.Platform$file.sep),format="gml")
	SPADE.plot.trees(mst_graph, output_dir, out_dir=paste(output_dir,"pdf",sep=.Platform$file.sep), layout=SPADE.layout.arch(mst_graph))
}
