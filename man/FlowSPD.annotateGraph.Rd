\name{FlowSPD.annotateGraph}
\alias{FlowSPD.annotateGraph}
\title{
Add attributes to graph
}
\description{
Add specific and arbitrary attributes to a graph
}
\usage{
FlowSPD.annotateGraph(graph, layout = NULL, ...)
}
\arguments{
  \item{graph}{
The graph object to work on. Note that the original graph is never modified, a new graph object is returned instead; if you don't assign it to a variable your modifications will be lost!
}
  \item{layout}{
Optional numeric matrix with vertex x,y positions with the same number of rows as vertices and at least two columns, the x and y positions.
}
  \item{\dots}{
Data matrices to add as attributes. All arguments must have the same number of rows as vertices. Column names are used as attribute names.
}
}
\details{
Add specific arbitrary attributes to a graph. 
}
\value{
A new graph object with the attributes added.
}
\references{
}
\author{
Michael Linderman
}
\note{
}
\seealso{
\code{\link{set.graph.attribute}}, \code{\link{set.vertex.attribute}}, \code{\link{set.edge.attribute}}
}
\examples{
}
