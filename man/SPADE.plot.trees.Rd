\name{SPADE.plot.trees}
\alias{SPADE.plot.trees}
\title{
Plot trees with annotated vertices
}
\description{
Plot trees for each vertex annotation setting vertex size and color based on the particular annotation.
}
\usage{
SPADE.plot.trees(graph, files, file_pattern = "*anno.Rsave", out_dir = ".", layout = SPADE.layout.arch, attr_pattern = "percent|medians|fold|cvs", scale = NULL, pctile_color=c(0.02,0.98), normalize="global", size_scale_factor=1, edge.color="grey", bare=FALSE, palette="bluered")
}
\arguments{
  \item{graph}{
iGraph graph object
}
  \item{files}{
Either a vector of save annotation files or a directory. If a directory, all of the files matching the \code{pattern} wildcard pattern are processesd.
}
  \item{file_pattern}{
Wildcard pattern to match files if \code{files} is a directory.
}
  \item{out_dir}{
Directory where output files are written. Will be created if it does not exist.
}
  \item{layout}{
Either a function or a numeric matrix specifying how vertices are placed on plot. If it is a matrix, the matrix must have two columns, x and y position, and as many rows as vertices. If \code{layout} is a function, it will be called with an igraph graph as the single parameter.
}
  \item{attr_pattern}{
A regular expression that matches the attributes that should be plotted for each graph.  Parameter names matching regex "median|fraction|cvs" will be plotted with a scale range set to [min, max] for that attribute, while all other parameters will be plotted on a centered scale with the range [abs(min(parameter values),max(parameter values)), abs(min(parameter values),max(parameter values))]
}
  \item{scale}{
A two element vector, e.g. c(-1,1), specifying low and upper bound for color scale. Values below and above these bounds will be forced to the 'smallest' and 'largest' color respectively. If specified, overrides 'pctile_color'.
}
  \item{pctile_color}{
A two element vector specifying lower and upper percentiles that should be used to set the color scale. Values below and above these percentiles will be forced to the 'smallest' and 'largest' color respectively. Not in effect if 'scale' is specified.
}
  \item{normalize}{
A string (either "global" or "local"), specifying color scale normalization.  Setting to "global" will set the scale range to the global min/max of all GML files in the folder, while "local" will set the scale range to the min/max of the particular GML file being plotted.
}
  \item{size_scale_factor}{
A scale factor for node size in drawing. Current function for node size: \code{percenttotal[i]/(max(percenttotal)^(1/size_scale_factor)) * 3 + 2}
}
  \item{edge.color}{
Set the edge color. See \code{igraph.plotting} for more details.
}
  \item{bare}{
Boolean specifying whether to omit titles and gradient legend.
}
  \item{palette}{
A string (either "jet" or "bluered"), specifying color palette for nodes. "bluered" tends to show up better on LCD projectors.
}
}
\author{
Michael Linderman
}
\seealso{
	\code{\link{SPADE.driver}}
}
\examples{
	# Load two-parameters sample data included in package
	data_file_path = paste(installed.packages()["spade","LibPath"],"spade","extdata","SimulatedRawData.fcs",sep=.Platform$file.sep)

	# Run basic SPADE analyses, clustering on two parameters. 
	output_dir <- tempdir()
	SPADE.driver(data_file_path, out_dir=output_dir, cluster_cols=c("marker1","marker2"))

	# Generate PDFs of annotated graphs (into output_dir/pdf)
	mst_graph <- igraph0:::read.graph(paste(output_dir,"mst.gml",sep=.Platform$file.sep),format="gml")
	SPADE.plot.trees(mst_graph, output_dir, out_dir=paste(output_dir,"pdf",sep=.Platform$file.sep), layout=igraph0:::layout.kamada.kawai(mst_graph))
}
