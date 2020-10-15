\name{update_annotations}
\alias{update_annotations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Update Sample Annotations
}
\description{
Update sample annotations stored in annotations parameter of the params list object. See \link[dataVisEasy]{set_annotations} for intially setting the annotations.
}
\usage{
update_annotations(annotation, values)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annotation}{
title of annotation to be updated or added. If already found in the colnames of the annotations data frame, the previous values for that column will be replaced with the supplied values here. If not found in the colnames of the annotations data frame, values will be added as a new annotation column under the name indicated.
}
  \item{values}{
values corresponding to the annotation to be updated or added.
}
}
\details{
Length of values must be equal to number of rows in original annotations dataframe. See all.samples option in \link[dataVisEasy]{extractClusters} if added annotations are identifications of cluster groups.  If new or updated values should be shown in the heatmap, must rerun \link[dataVisEasy]{set_annot_samps}. To add specific colors for the new annotations, see \link[dataVisEasy]{set_annot_cols} and \link[dataVisEasy]{update_annot_cols}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (annotation, values)
{
    if (annotation \%in\% colnames(annotations)) {
        params$annotations[, which(colnames(params$annotations) ==
            annotation)] <<- values
    }
    else {
        params$annotations$V1 <<- values
        colnames(params$annotations)[colnames(params$annotations) ==
            "V1"] <<- annotation
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
