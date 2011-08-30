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
