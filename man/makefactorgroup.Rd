\name{makefactorgroup}
\alias{makefactorgroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a Factor Based On Ordered Annotations
}
\description{
This function is called in \link[dataVisEasy]{myHeatmapByAnnotation} to create a factor with proper factor levels based on one or more annotations. If more than one level is provided, the function will order the provided sample annotations by the first level first, followed by subsequent levels to obtain an ordered grouping.
}
\usage{
makefactorgroup(annots, levels, specify.gaps = NULL, return.gaps = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annots}{
a data frame of annotations
}
  \item{levels}{
a character vector of length 1-3 where each element must be found in the colnames of annots. The first element will be the first annotation to be sorted. Ex. c("Age","Strain","Treatment") will sort first by Age, then by Strain, then by Treatment.
}
  \item{specify.gaps}{
numeric vector the same length as a character vector supplied to levels. The default of NULL will result in a gap in between each annotation specified. Mostly useful in scenarios where 2 or 3 levels are provided and allows flexibility in how they will be displayed. Numeric values provided also work in tandem with each other. For example a vector of c(1,1) for a supplied groupings vector of length 2 will place 1 space in between levels of the first annotation and an additional space in between levels of the second annotaiton. A vector of c(0,1) will still order the samples by the annotations as indicated but will only place gaps in between the levels of the second annotation supplied.
}
  \item{return.gaps}{
logical. Should gaps a vector of gaps be returned. Necessary when called upon in the \link[dataVisEasy]{myHeatmapByAnnotation}.
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
if return.gaps is set to FALSE (default) will return a factor the same length as the number of rows in the provided annots

if return.gaps is set to TRUE, will return a list object where "factor.group" will contain the factor indicated above and "gaps" will contain a numeric vector of where gaps should be placed in a heatmap according to the specified annotations and spacing.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
~~Alison Moss~~
}
\note{
This function is mainly used internally by \link[dataVisEasy]{myHeatmapByAnnotation} to organize heatmaps by multiple annotations.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##initiate_parameters
initiate_params()

factorgroupings <- makefactorgroup(RAGP_annots, c("State","Connectivity"))
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
