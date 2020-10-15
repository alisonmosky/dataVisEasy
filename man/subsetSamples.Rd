\name{subsetSamples}
\alias{subsetSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subset Samples/observations from a data matrix
}
\description{
Quick function to subset Samples/observations from a data matrix based on an annotation or characteristic
}
\usage{
subsetSamples(data, group, take.out)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows
}
  \item{group}{
character string pointing to a column name in the annotations stored in params. Can also supply a custom vector the same length as the columns of the input matrix.
}
  \item{take.out}{
level or levels of the group vector supplied indicating the samples to be extracted
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
data matrix with the same number of rows as the input data matrix but subset for the specified samples
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
function (data, group, take.out)
{
    temp.annotations <- params$annotations
    if (group \%in\% colnames(temp.annotations)) {
        if (sum(colnames(data) \%notin\% rownames(temp.annotations)) !=
            0) {
            stop("colnames of input data do not match rownames of annotations, cannot link annotations to data")
        }
        temp.annotations <- temp.annotations[match(colnames(data),
            rownames(temp.annotations)), , drop = FALSE]
        groupings <- as.factor(temp.annotations[, group])
        if (sum(take.out \%in\% groupings) != length(take.out)) {
            stop("provided arguments to take.out not found in indicated group")
        }
        set <- data[, which(groupings \%in\% take.out)]
    }
    else {
        if (sum(take.out \%in\% group) != length(take.out)) {
            stop("provided arguments to take.out not found in indicated group")
        }
        set <- data[, which(group \%in\% take.out)]
    }
    return(set)
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
