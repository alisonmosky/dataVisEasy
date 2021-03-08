\name{AOV1way}
\alias{AOV1way}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
One-Way ANOVA
}
\description{
A one-way ANOVA for the samples in the supplied data matrix based on annotations stored in the params list object.
}
\usage{
AOV1way(data.to.aov, category, pthreshold = 0.05, additional.report = "NONE")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.to.aov}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows
}
  \item{category}{
Character string pointing to a column of the annotations data frame stores in params$annotations
}
  \item{pthreshold}{
pvalue threshold for which genes will be considered significant
}
  \item{additional.report}{
Option to include full output of aov or TukeyHSD, options are "NONE", "AOV","TUKEY", or "ALL"
}
}
\details{
Runs a standard one-way ANOVA and subsequent TukeyHSD on the data provided using stored annotations as the categories for ANOVA groupings
}
\value{
~List object

\item{AOV.Results}{dataframe giving the F statistic and pvalue of each gene in the original input matrix }

\item{Sig.Genes}{character vector of the genes with pvalues below the specified pthreshold}

\item{NonSig.Genes}{character vector of the genes with pvalues that did not pass the specified pthreshold}

\item{Tukey.pVals}{dataframe containing the adjusted pvalues from TukeyHSD for each contrast for each gene in the original input matrix}

\item{Tukey.diffs}{dataframe containing the difference values from TukeyHSD for each contrast for each gene in the original input matrix}

optional items

\item{AOV.output}{list object containing the aov output for all genes in the original input matrix }

\item{Tukey.output}{list object containing the TUkeyHSD output for all genes in the original input matrix }
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
 ~~Alison Moss~~
}
\note{
The function itself determines significant and non-significant genes based on the pthreshold
provided. However, the AOV.Results provdied in the output supply the F value and p value,
additional information can also be retrieved by setting additional report = TRUE, user can
use any of this information for custom filtering of results.

For cases where there are no non-missing arguments for one or more of the groups in the supplied categories, a warning will be generated indicating which genes or observations this is true for. If there is only one group with no non-missing arguments, that gene or observation will be removed from the analysis (and will be indicated with a warning). Please pay attention to these warnings as they may affect how the results should be interpretted.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##initiate parameteres and set up annotations

initiate_params()
data(RAGP_norm)
data(RAGP_annots)
set_annotations(RAGP_annots)

##AOV for Connectivity
state.aov <- AOV1way(RAGP_norm, "State")
head(state.aov$AOV.Results)
head(state.aov$Sig.Genes)

##visualize significant genes in heatmap
set_annot_samps("State")
myHeatmapByAnnotation(RAGP_norm, state.aov$Sig.Genes, groupings = "State")


##AOV for Animal
animal.aov <- AOV1way(RAGP_norm, "Animal", pthreshold = 0.01)
summary(animal.aov)
head(animal.aov$Tukey.pVals) ##look at pvalues from tukey post hoc analysis

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
