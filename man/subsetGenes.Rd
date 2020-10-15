\name{subsetGenes}
\alias{subsetGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subet Genes or Variables from a Dataset
}
\description{
A quick function to pull out genes or variables from a data matri with the option to order the columns by expression of a given gene.
}
\usage{
subsetGenes(data, list, order.by = NULL, exact = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows
}
  \item{list}{
character vector of genes/variables to be pulled out of the data matrix for viewing
}
  \item{order.by}{
optional string equal to one of the rownames of the supplied data frame. Columns of output data matrix will be ordered according to the supplied gene from lowest to highest expression value.
}
  \item{exact}{
whether or not to search for exact or inexact matches of 'list' in 'data'. If exact = T (default) heatmap will plot genes/variables that exactly match the list supplied. If set to FALSE, will search for inexact matches.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
data matrix with the same number of columns as the input data matrix but subset for the specified gene/variable list
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
function (data, list, order.by = NULL, exact = TRUE)
{
    if (exact == TRUE) {
        subset <- data[which(rownames(data) \%in\% list), ]
        if (length(subset) == 0) {
            stop("exact matches for list not found in rownames data")
        }
    }
    if (exact == FALSE) {
        subset <- data[grep(paste(list, collapse = "|"), rownames(data)),
            ]
        if (length(subset) == 0) {
            stop("inexact matches for list not found in rownames data")
        }
    }
    if (is.null(order.by) == FALSE) {
        subset <- subset[, order(data[which(rownames(data) \%in\%
            order.by), ], na.last = F)]
    }
    return(subset)
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
