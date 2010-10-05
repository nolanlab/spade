\name{SPADE.layout.arch}
\alias{SPADE.layout.arch}
\alias{SPADE.layout.arch_layout}
\title{
Generate coordinates for plotting graphs
}
\description{
Performing "arch" layouts of graph vertices
}
\usage{
SPADE.layout.arch(mst_graph)
SPADE.layout.arch_layout(mst_graph)
}
\arguments{
  \item{mst_graph}{
The graph to layout. Must be acyclic and undirected.
}
}
\details{
These functions calculate the coordinates of the vertices for a graph.

\code{layout.arch} lays out the longest chain of the graph, the "backbone", on an arch, and the "side chains" as trees normal to that backbone.

\code{layout.arch_layout} similarly lays out the backbone on an arch, but lays out the side chains using a modified Fruchterman and Reingold force directed algorith.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
All these functions return a numeric matrix with at least two columns, x and y positions, and the same number of lines as the number of vertices.
}
\author{
Michael Linderman
}
\note{
}
\seealso{
\code{\link{SPADE.annotateGraph}}
}
\examples{
}
