\name{set_scale.colors}
\alias{set_scale.colors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set Colors to Be Used in Heatmaps Heatmaps
}
\description{
Set and override the default colors to be used in heatmaps. Stored in the params list object, the default is c("blue","black","yellow")
}
\usage{
set_scale.colors(colors)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{colors}{
vector of colors
}
}
\details{
colors used to generate heatmaps

Please see \link[dataVisEasy]{params} for more information about accessing and setting up these variables.
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
\link[dataVisEasy]{set_scale.range} and \link[dataVisEasy]{set_n.colors.range}
}
\examples{
##initiate parameters
initiate_params()

#view heatmap with default colors
myHeatmap(RAGP_norm)

##switch to green-black-red heatmap
set_scale.colors(c("green","black","red"))
myHeatmap(RAGP_norm)

##custom color range
set_scale.colors(c("purple","blue","green","yellow","orange","red"))
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
