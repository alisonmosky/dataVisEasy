\name{update_annot_cols}
\alias{update_annot_cols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Update Annotation Colors
}
\description{
Update sample annotations stored in annotations parameter of the params list object. See \link[dataVisEasy]{set_annotations} for intially setting the annotations.

Update colors corresponding to both sample and gene annotations stored in the annotations and annotations.genes parameters of the params list object that are passed to a variety of functions in the dataVisEasy package. See \link[dataVisEasy]{set_annot_cols} for intially setting the annotations.
}
\usage{
update_annot_cols(annotation, values.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annotation}{
title of annotation to be updated or added. If already found in the names of the annot_cols list, the previous values for that annotation will be replaced with the supplied values here. If not found in the names of the annot_cols list, values will be added as a new item under the name indicated.
}
  \item{values.list}{
values corresponding to the annotation to be updated or added.
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

#genemod.cols <- brewer.pal(7,"Set3"); names(genemod.cols) <- LETTERS[1:7]
#update_annot_cols("GeneMods",c(genemod.cols))


##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (annotation, values.list)
{
    if (annotation \%in\% names(params$annot_cols)) {
        params$annot_cols[which(names(params$annot_cols) == annotation)] <<- values.list
    }
    else {
        params$annot_cols$V1 <<- values.list
        names(params$annot_cols)[names(params$annot_cols) ==
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
