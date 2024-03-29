\name{set_annot_genes}
\alias{set_annot_genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set Gene Annotations to Appear Along Left Side of Heatmaps
}
\description{
Specify which gene annotations should appear along the left side of heatmaps
}
\usage{
set_annot_genes(annotations = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annotations}{
a character vector including the column names in the annotations.genes parameter (\link[dataVisEasy]{set_annotations.genes}) that should be shown as annotations in the heatmap. If no argument is supplied, assumes that all columns in the annotations.genes parameter should be shown in the heatmap. For showing multiple annotations, the order in which the annotations are listed is the order they will appear on the heatmap from the closest to the heatmap on the right to those farthest on the left. To remove all annotation tracks and reset, simply set to NA.
}
}
\details{
Draws from the annotations.genes parameter stored in the params list object to hold only the annotations that should be shown in heatmaps. See \link[dataVisEasy]{set_annotations.genes}, must be run prior to this function.

To specify the colors that will accompany each level of each annotation, see \link[dataVisEasy]{set_annot_cols}. If colors are not specified, default color values will be applied when necessary. (Not all annotations need to have custom colors, those that do not will be given default colors while those that are specified will be shown as such.)
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
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \link[dataVisEasy]{set_annotations.genes} and \link[dataVisEasy]{set_annot_cols}
Works the same as for sample annotations, see \link[dataVisEasy]{set_annot_samps} for examples
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
