\name{set_annot_samps}
\alias{set_annot_samps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set Sample Annotations to Appear Along Top of Heatmaps
}
\description{
Specify which sample annotations should appear along the top of heatmaps
}
\usage{
set_annot_samps(annotations = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annotations}{
a character vector including the column names in the annotations parameter (\link[dataVisEasy]{set_annotations}) that should be shown as annotations in the heatmap. If no argument is supplied, assumes that all columns in the annotations should be shown in the heatmap.
}
}
\details{
Draws from the annotations parameter stored in the params list object to hold only the annotations that should be shown in heatmaps. See \link[dataVisEasy]{set_annotations}, must be run prior to this function.

To specify the colors that will accompany each level of each annotation, see \link[dataVisEasy]{set_annot_cols}. If colors are not specified, default color values will be applied when necessary. (Not all annotations need to have custom colors, those that do not will be given default colors while those that are specified will be shown as such.)
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
~~Alison Moss~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \link[dataVisEasy]{set_annotations} and \link[dataVisEasy]{set_annot_cols}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (annotations = NULL)
{
    if (is.null(annotations) == TRUE) {
        params$annot_samps <<- params$annotations
    }
    else {
        suppressWarnings(if (is.na(annotations)) {
            params$annot_samps <<- NA
        }
        else {
            params$annot_samps <<- params$annotations[, which(colnames(params$annotations) \%in\%
                annotations), drop = FALSE]
            params$annot_samps <<- params$annot_samps[, match(annotations,
                colnames(params$annot_samps)), drop = FALSE]
        })
    }
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
