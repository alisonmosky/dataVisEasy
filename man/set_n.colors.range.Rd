\name{set_n.colors.range}
\alias{set_n.colors.range}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set Number of Colors Used in Heatmaps
}
\description{
A necessary parameter for the generation of heatmaps with custom breaks and color scale as is done in \link[dataVisEasy]{myHeatmap} and \link[dataVisEasy]{myHeatmapByAnnotation}. See note for how this is used.
}
\usage{
set_n.colors.range(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
a single numeric
}
}
\details{
this value is used in creating breaks and color scale for heatmaps. A sequence is created along the scale range set by \link[dataVisEasy]{set_scale.range} of length n.colors.range set in this function. n.colors.range - 1 colors will then be use to generate the colors in the heatmap. (See notes for the code where this is used)

Please see \link[dataVisEasy]{params} for more information about accessing and setting up these variables.
}
\value{

}
\references{
%% ~put references to the literature/web site here ~
}
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
\link[dataVisEasy]{set_scale.range} and \link[dataVisEasy]{set_scale.colors}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n)
{
    params$n.colors.range <<- n
  }
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
