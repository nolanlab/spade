\name{SPADE.normalize.trees}
\alias{SPADE.normalize.trees}
\title{
Normalize tree attributes
}
\description{
Normalize SPADE attributes, e.g. medians, fold-change, within a file or across files.
}
\usage{
SPADE.normalize.trees(files, file_pattern = "*.gml", out_dir = ".", layout = SPADE.layout.arch, attr_pattern = "percent|median|fold|ratio", normalize = "global")
}
\arguments{
  \item{files}{
Either a vector of gml graph files or a directory. If a directory, all of the files matching the \code{pattern} wildcard pattern are processesd.
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
A regular expression that matches the attributes that should be plotted for each graph.
}
  \item{normalize}{
One of 'local' or 'global'. If 'local', attributes are normalized on a per-file, per-attributed basis into the range [-1,1] with zero remaining zero (e.g., [-5,10] will normalize to [-5,1]). If 'global', the same per-attribute normalization is performed, but across all input files.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Michael Linderman
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
