\name{set_scale.range}
\alias{set_scale.range}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set Scale Range of Heatmaps
}
\description{
Set and override the default scale to be used in heatmaps. Stored in the params list object, the default is c(-1,1). How well the data fits into these limits can be assessed using \link[dataVisEasy]{assessScale}
}
\usage{
set_scale.range(range)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{range}{
numeric of length 2
}
}

\details{
scale of heatmaps
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}

%\references{
%% ~put references to the literature/web site here ~
%}
\author{
~~Alison Moss~~
}

\note{
In \link[dataVisEasy]{myHeatmap} and \link[dataVisEasy]{myHeatmapByAnnotation}, breaks and colors are set using these and other parameters.
  breaks <- seq(params$scale.range[1], params$scale.range[2],length.out=params$n.colors.range)
  my_cols=colorRampPalette(params$scale.colors)(n=params$n.colors.range-1)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[dataVisEasy]{set_scale.colors} and \link[dataVisEasy]{set_n.colors.range}

See \link[dataVisEasy]{assessScale} to see how the data fits into the set scale
}

\examples{
#initiate parameters
initiate_params()

##view heatmap with default range of -1 to 1
myHeatmap(RAGP_norm)

##set scale to -2 to 2, view with new range
set_scale.range(c(-2,2))
myHeatmap(RAGP_norm)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
