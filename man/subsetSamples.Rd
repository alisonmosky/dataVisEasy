\name{subsetSamples}
\alias{subsetSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subset Samples/observations from a data matrix
}
\description{
Quick function to subset samples/observations from a data matrix based on an annotation or characteristic. Can also pull out samples/observations by name with the option to order the rows by expression of a given sample.
}
\usage{
subsetSamples(data, group, take.out, list, order.by = NULL, exact = F,
subset.by = "annotation")
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
  \item{list}{
character vector of samples/observations to be pulled out of the data matrix for viewing
}
  \item{order.by}{
optional string equal to one of the colnames of the supplied data frame. Rows of output data matrix will be ordered according to the supplied sample from lowest to highest expression value.
}
  \item{exact}{
whether or not to search for exact or inexact matches of 'list' in 'data'. If exact = T (default) output will contain samples/observations that exactly match the list supplied. If set to FALSE, will search for inexact matches.
}
  \item{subset.by}{
whether to subset by an annotation stored in params$annotations ("annotation", default) or by the a list of names to search for in column names ("name")
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
data matrix with the same number of rows as the input data matrix but subset for the specified samples
}
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##initiate parameters and set up annotations
initiate_params()
set_annotations(RAGP_annots)

##subset samples for SAN-Projecting cells
SAN.proj <- subsetSamples(RAGP_norm, group = "Connectivity", take.out = "SAN-Projecting")

##subset samples for cells with projection information
proj.ino <- subsetSamples(RAGP_norm, group = "Connectivity",
  take.out = c("SAN-Projecting", "Non-SAN-Projecting"))
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
