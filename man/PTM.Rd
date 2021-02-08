\name{PTM}
\alias{PTM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pavlidis Template Matching
}
\description{
a simple implementation of Pavlidis Template Matching that searches for correlations above a specified threshold within the supplied input data. Templates can be annotations, genes or samples, or a cusotm template can be provided.
}
\usage{
PTM(data, match.template, annotation.level.set.high, custom.template = NA,
  Find.Match.For = "genes",cutoff = 0.05, cut.by = "pvals", method = "pearson",
  NA.handling = "pairwise.complete.obs", return.vals = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows
}
  \item{match.template}{
template to be matched against. If Find.Match.For = "genes", can supply a gene to be matched against (rownames in input data) or a string pointing to an annotation stored in the annotations dataframe in the params list object. If Find.Match.For = "samples", can supply a sample to be matched against (colnames in input data) or a string pointing to an annotation stored in the annotations.genes dataframe in the params list object.
}
  \item{annotation.level.set.high}{
a character vector indicating which levels of an annotation factor to set high. match.template will point to an annotation stored in the params list object, levels supplied here will be set to 1 while all other levels will be set to 0, searching for matches that have high expression for the given levels and low expression for all others
}
  \item{custom.template}{
a custom template to be matched against. Must be the same length as columns of input data if Find.Match.For = "genes" or the same length as rows of input data if Find.Match.For = "samples"
}
  \item{Find.Match.For}{
whether matches should be found for genes or samples, accepts "genes" (default) or "samples"
}
  \item{cutoff}{
the cutoff for which matches are determined to be significant. should be set according to cut.by. Default of 0.05 assumes cut.by = "pvals".
}
  \item{cut.by}{
whether to use the pvalues ("pvals", default) or correlations ("rvals") as the cutoff for significant matches. If "pvals", will return matches that have a pvalue below "cutoff" and with a correlation > 0 (only positive correlations not negative correlation). If "rvals", will return matches that have correlation > cutoff.
}
  \item{method}{
method to be used to correlation. Accepts values to be passed to cor() such as "pearson" (default), "spearman", and "kendall".
}
  \item{NA.handling}{
how missing values should be handled in the case of correlations, passed to the "use" argument of cor()
}
  \item{return.vals}{
logical of whether or not to return pvalues and correlation values for all comparisons. If FALSE (default), will return a list of genes or samples that passed the threshold set by cutoff. If TRUE, will return a dataframe containing the pvalues and correlations for each item with respect to the template.
}
}
\details{
Will compute all correlations for the rows/columns of the input data compared to a provided template using the indicated correlation method. Custom templates the same length as the rows/columns of the data can be provided or matched can be found to a specific gene (see \link[dataVisEasy]{corrs2Gene}) or samples. Additionally, templates can be created based off of sample or gene annotations stored in the params list object.If pointing to an annotation stored in the params list object, levels supplied to "set.high" will be set to 1 while all other levels will be set to 0, searching for matches that have high expression for the given levels and low expression for all others.
}
\value{
if return.vals == FALSE, will return a character vector of genes or samples that passed the pvalue or correlation cutoff indicated.

if return.vals == TRUE, will return a dataframe giving the pvalue and correlation value for each variable compared to the supplied template
}
\references{
Pavlidis and Noble 2001
}
\author{
~~Alison Moss~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \link[dataVisEasy]{params},\link[dataVisEasy]{set_annotations}, \link[dataVisEasy]{set_annotations.genes}, for more information on setting up annotations.
}
\examples{
##initiate parameters and annotations
initiate_params()
set_annotations(RAGP_annots)

##using the "State" annotation, find genes were State A has high expression, return passes
A.high <- PTM(RAGP_norm, match.template = "State", annotation.level.set.high = "A")
head(A.high)
A.high.01 <- PTM(RAGP_norm, match.template = "State", annotation.level.set.high = "A",
  cutoff = 0.01)
head(A.high.01)

##return dataframe with pvalues and correlations to template
A.high.results <- PTM(RAGP_norm, match.template = "State", annotation.level.set.high = "A",
  , return.vals = TRUE)
head(A.high.results)


##find matches to a gene of interest
Th.match <- PTM(RAGP_norm, match.template = "Th", cut.by = "rvals", cutoff = 0.8)
Th.match

##the same can be applied to find matches for samples
##can also supply a custom vector to be matched against
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
