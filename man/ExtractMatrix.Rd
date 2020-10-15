\name{ExtractMatrix}
\alias{ExtractMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the Data Matrix from an Output Heatmap
}
\description{
Allows the user to extract the exact order of a heatmap from the output. Useful for preserving the ordering of samples or genes after clustering or often helpful for supplying supplemental information in papers.
}
\usage{
ExtractMatrix(data, heatmap, clustered.cols = TRUE, clustered.rows = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows. Should be the same data supplied to the heatmap function
}
  \item{heatmap}{
Output of a heatmap generated from \link[dataVisEasy]{myHeatmap} or \link[dataVisEasy]{myHeatmapByAnnotation}.  The heatmap generated from these functions can be saved to a variable which should serve as the input here. IMPORTANT. The heatmap supplied here MUST have show.colnames and show.rownames = TRUE.
}
  \item{clustered.cols}{
logical. should the order of the columns be extracted
}
  \item{clustered.rows}{
logical. should the order of the rows be extracted
}
}
\details{
The heatmap supplied here MUST have show.colnames and show.rownames = TRUE. If the generated heatmap was a subset of the data provided as the input (a list of genes was given) this function will subset for those genes in the output.
}
\value{
a data matrix of the same dimensions as the input heatmap
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
function (data, heatmap, clustered.cols = TRUE, clustered.rows = TRUE)
{
    attributes <- unlist(lapply(heatmap$gtable$grobs, function(x) (x$name)))
    text.attributes <- grep("GRID.text", attributes)
    ordered <- data
    if (clustered.cols == TRUE) {
        colord <- heatmap$gtable[["grobs"]][[text.attributes[2]]][["label"]]
        ordered <- ordered[, match(colord, colnames(ordered))]
    }
    if (clustered.rows == TRUE) {
        roword <- heatmap$gtable[["grobs"]][[text.attributes[3]]][["label"]]
        ordered <- ordered[match(roword, rownames(ordered)),
            ]
    }
    return(ordered)
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
