\name{SPADE.annotateGraph}
\alias{SPADE.annotateGraph}
\title{
Add attributes to graph
}
\description{
Add specific and arbitrary attributes to a graph
}
\usage{
SPADE.annotateGraph(graph, layout = NULL, anno)
}
\arguments{
  \item{graph}{
The graph object to work on. Note that the original graph is never modified, a new graph object is returned instead; if you don't assign it to a variable your modifications will be lost!
}
  \item{layout}{
Optional numeric matrix with vertex x,y positions with the same number of rows as vertices and at least two columns, the x and y positions.
}
  \item{anno}{
List of annotations to add to the graph. Each entry in list must have a name and must be a matrix. All matrices must have the same number of rows as vertices. List entry name plus column names are used as attribute names (unless they match, then just the column name is used).
}
}
\details{
Add specific arbitrary attributes to a graph. 
}
\value{
A new graph object with the attributes added.
}
\author{
Michael Linderman
}
\seealso{
\code{\link{set.graph.attribute}}, \code{\link{set.vertex.attribute}}, \code{\link{set.edge.attribute}}
}
\examples{
	## Load two-parameters sample data included in package
	#data_file_path = paste(installed.packages()["spade","LibPath"],"spade","extdata","SimulatedRawData.fcs",sep=.Platform$file.sep)

	## Run basic SPADE analyses, clustering on two parameters. Annotated graphs will be
	## in output_dir. See SPADE.plot.trees to generate PDFs of annotated graphs.
	#output_dir <- tempdir()
	#SPADE.driver(data_file_path, out_dir=output_dir, cluster_cols=c("marker1","marker2"))

	## Add additional parameters to output graphs using SPADE.annotateGraph
	#old_graph <- igraph0:::read.graph(paste(output_dir,"SimulatedRawData.fcs.density.fcs.cluster.fcs.medians.gml",sep=.Platform$file.sep),format="gml")
	#new_graph <- SPADE.annotateGraph(old_graph, layout=igraph0:::layout.kamada.kawai(old_graph), anno=list(demo=matrix(1:igraph0:::vcount(old_graph))))
}
