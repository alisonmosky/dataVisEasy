\name{set_expression_gradient.colors}
\alias{set_expression_gradient.colors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set Colors Used In Gene Expression Gradients
}
\description{
Set and override the default colors to be used by \link[dataVisEasy]{myColorRamp5} when generating colors based on expression levels. Default values are c("blue","lightblue","gray","indianred","firebrick")
}
\usage{
set_expression_gradient.colors(colors)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{colors}{
vector of five colors
}
}
\details{
Storing the colors in the params list object allows them to be accessed by all functions to keep consistent across data analysis types. Utilized in \link[dataVisEasy]{scatterGenes} and \link[dataVisEasy]{beeswarmGenes} in particular as well as \link[dataVisEasy]{myPCA}

Please see \link[dataVisEasy]{params} for more information about accessing and setting up these variables.
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
Example code where these colors are used:
   cols <- myColorRamp5(params$expression_gradient.colors,genedat, percent.mad = percent.mad)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[dataVisEasy]{myColorRamp5}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (colors)
{
    if (length(colors) == 5) {
        is.color = sapply(colors, function(x) {
            tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE)
        })
        if (any(is.na(names(is.color))))
            is.color[is.na(names(is.color))] = FALSE
        if (sum(is.color) == 5) {
            params$expression_gradient.colors <<- colors
        }
        else {
            stop("must contain valid colors")
        }
    }
    else {
        stop("must be a vector of colors of length 5")
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
