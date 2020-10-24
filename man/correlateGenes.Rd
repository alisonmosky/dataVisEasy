\name{correlateGenes}
\alias{correlateGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Examine Gene-Gene Correlations in Datamatrix
}
\description{
Provides a quick way to check for the range of correlations in a given data matrix. If no limits are supplied, will return a histogram of all correlation present. When limits are supplied, a data frame is returned detailing all gene pair correlations that are more negatively correlated than the lower limit and more positively correlated than the upper limit.
}
\usage{
correlateGenes(data, limits = NULL, nbreaks = 20, method = "pearson",
  NA.handling = "pairwise.complete.obs")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows
}
  \item{limits}{
vector of two numbers, giving the lower and upper bounds of the returned correlations. If NULL (the default) will return a histogram of all correlations in the data. When limits are provided, correlations below the lower limit and above the upper limit will be returned.
}
  \item{nbreaks}{
number of breaks for the histogram if limits is set to NULL
}
  \item{method}{
method to be used to calculate the correlation. Accepts values to be passed to cor() such as "pearson" (default), "spearman", and "kendall".
}
  \item{NA.handling}{
how missing values should be handled in the case of correlations, passed to the "use" argument of cor()
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
If limits are not provided, returns a plot histogram. If limits are provided, returns a data frame with three columns where the first two columns indicate the gene pair and the third column gives the correlation between them.
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
function (data, limits = NULL, nbreaks = 20, method = "pearson",
    NA.handling = "pairwise.complete.obs")
{
    main = paste("Gene Correlations")
    cor.dat <- cor(t(data), use = NA.handling, method = method)
    cor.dat[!upper.tri(cor.dat)] <- NA
    if (is.null(limits) == TRUE) {
        hist(cor.dat, breaks = nbreaks, main = main)
    }
    if (is.null(limits) == FALSE) {
        below <- limits[1]
        above <- limits[2]
        picked.cors <- which((cor.dat > above & cor.dat < 0.99) |
            cor.dat < below, arr.ind = T)
        cor.genes <- as.matrix(apply(picked.cors, 2, function(x) (colnames(t(data))[x])))
        genes.and.cors <- as.data.frame(cbind(cor.genes, Correlation = cor.dat[picked.cors]))
        colnames(genes.and.cors) <- c("Gene1", "Gene2", "Correlation")
        return(genes.and.cors)
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
