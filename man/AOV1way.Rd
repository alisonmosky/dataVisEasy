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
function (data.to.aov, category, pthreshold = 0.05, additional.report = "NONE")
{
    if (sum(colnames(data.to.aov) \%notin\% rownames(params$annotations)) !=
        0) {
        stop("colnames of input data do not match rownames of annotations, cannot link annotations to data and assign groupings for ANOVA")
    }
    temp.annotations <- params$annotations[match(colnames(data.to.aov),
        rownames(params$annotations)), ]
    groupings <- temp.annotations[, category]
    aov.all <- apply(data.to.aov, 1, function(x) (summary(aov(x ~
        groupings))))
    aov.results <- data.frame(FVal = unlist(lapply(aov.all, function(x) ((x[[1]]$`F value`[1])))),
        pVal = unlist(lapply(aov.all, function(x) ((x[[1]]$`Pr(>F)`[1])))),
        row.names = names(aov.all))
    sig.genes <- rownames(aov.results)[which(aov.results$pVal <=
        pthreshold)]
    nonsig.genes <- rownames(aov.results)[which(aov.results$pVal >
        pthreshold)]
    sig.set <- data.to.aov[which(rownames(data.to.aov) \%in\% sig.genes),
        ]
    tukey.all <- apply(sig.set, 1, function(x) (TukeyHSD(aov(x ~
        groupings))))
    tukey.pvals <- data.frame(Reduce(rbind, lapply(tukey.all,
        function(x) (x$groupings[, 4]))), row.names = names(tukey.all))
    colnames(tukey.pvals) <- gsub("\\.", "-", colnames(tukey.pvals))
    tukey.diffs <- data.frame(Reduce(rbind, lapply(tukey.all,
        function(x) (x$groupings[, 1]))), row.names = names(tukey.all))
    colnames(tukey.diffs) <- gsub("\\.", "-", colnames(tukey.diffs))
    if (toupper(additional.report) == "ALL") {
        return(list(AOV.output = aov.all, AOV.Results = aov.results,
            Sig.Genes = sig.genes, NonSig.Genes = nonsig.genes,
            Tukey.output = tukey.all, Tukey.pVals = tukey.pvals),
            Tukey.diffs = tukey.diffs)
    }
    if (toupper(additional.report) == "TUKEY") {
        return(list(AOV.Results = aov.results, Sig.Genes = sig.genes,
            NonSig.Genes = nonsig.genes, Tukey.output = tukey.all,
            Tukey.pVals = tukey.pvals, Tukey.diffs = tukey.diffs))
    }
    if (toupper(additional.report) == "AOV") {
        return(list(AOV.output = aov.all, AOV.Results = aov.results,
            Sig.Genes = sig.genes, NonSig.Genes = nonsig.genes,
            Tukey.pVals = tukey.pvals, Tukey.diffs = tukey.diffs))
    }
    if (toupper(additional.report) == "NONE") {
        return(list(AOV.Results = aov.results, Sig.Genes = sig.genes,
            NonSig.Genes = nonsig.genes, Tukey.pVals = tukey.pvals,
            Tukey.diffs = tukey.diffs))
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
