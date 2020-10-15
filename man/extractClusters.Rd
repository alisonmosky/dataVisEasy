\name{extractClusters}
\alias{extractClusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the clusters from a Clustered Heatmap
}
\description{
Provided a clustered heatmap, and a number of groups the dendrogram should be split into, will return a vector detailing which genes/samples are in each cluster. Particularly helpful after setting row.groups or column.groups in the heatmap functions and the groupings of the genes or samples in these groups need to be identified.
}
\usage{
extractClusters(data, heatmap, extractGenes = FALSE, GeneGroup_Name = NULL, extractSamples = TRUE, SampleGroup_Name = NULL, num_Gene.groups = 5, num_Sample.groups = 5, all.genes = NULL, all.samples = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows. Should be the same data supplied to the heatmap function
}
  \item{heatmap}{
Output of a heatmap generated from \link[dataVisEasy]{myHeatmap} or \link[dataVisEasy]{myHeatmapByAnnotation}.  The heatmap generated from these functions can be saved to a variable which should serve as the input here. IMPORTANT. The heatmap supplied here MUST have show.colnames and show.rownames = TRUE.
}
  \item{extractGenes}{
logical. should clusters of genes (clustered rows) be extracted. At least one of extractGenes and extractSamples must be set to TRUE
}
  \item{GeneGroup_Name}{
optional name for the groupings to be extracted
}
  \item{extractSamples}{
logical. should clusters of samples (clustered columns) be extracted. At least one of extractGenes and extractSamples must be set to TRUE
}
  \item{SampleGroup_Name}{
optional name for the groupings to be extracted
}
  \item{num_Gene.groups}{
number of clusters along the rows that should be extracted. Default set to 5.
}
  \item{num_Sample.groups}{
number of clusters along the columns that should be extracted. Default set to 5.
}
  \item{all.genes}{
gene names, if annotations should be in placed in larger gene list than matrix in question. Default assumes that the rownames of the input data represent all genes in question (this is true even if the genes in the heatmap are a subset of the input data). Genes not represented in the heatmap or input data that do not have a cluster designation will be assigned "_".
}
  \item{all.samples}{
sample names, if annotations should be in placed in larger sample list than matrix in question. Default assumes that the colnames of the input data represent all samples in question. If these annotations should be placed within a larger list (for example if the samples in the input data are a subset of all samples). Samples not represented in the heatmap or input data that do not have a cluster designation will be assigned "_".
}
}
\details{
This function is helpful for assigned annotations to samples or genes based on how they cluster. The options of all.samples and all.genes allow these cluster assignments to be placed within a larger list of variables (useful if combining these annotations with other pre-existing annotations, where the vectors need to be of the same length to merge correctly).

If GeneGroup_Name and/or SampleGroup_Name are provided, these vectors can be passed directly to \link[dataVisEasy]{update_annotations} or \link[dataVisEasy]{update_annotations.genes}
}
\value{
a named vector or a list object of two named vectors (in a case where both genes and samples are extracted) where the assignments of each gene/sample are designated by letter A:num.groups. Importantly, these groupings are assigned as they appear in the heatmap itself where cluster assignment A is located on the top left (going across from left to right in the case of extractSamples and going down from top to bottom in the case of extractGenes)
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
See \link[dataVisEasy]{set_annotations} and \link[dataVisEasy]{update_annotations.genes} or \link[dataVisEasy]{update_annotations} and \link[dataVisEasy]{update_annotations.genes} for how the output of this function can be passed directly to these functions to set or update annotations stored in the params list object
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, heatmap, extractGenes = FALSE, GeneGroup_Name = NULL,
    extractSamples = TRUE, SampleGroup_Name = NULL, num_Gene.groups = 5,
    num_Sample.groups = 5, all.genes = NULL, all.samples = NULL)
{
    if (is.null(all.genes) == TRUE) {
        all.genes <- rownames(data)
    }
    if (is.null(all.samples) == TRUE) {
        all.samples <- colnames(data)
    }
    if (extractGenes == TRUE) {
        report.gene.groups <- rep("_", length(all.genes))
        gene.groups <- cutree(heatmap$tree_row, k = num_Gene.groups)
        alph <- LETTERS[1:num_Gene.groups]
        gene.tree.order <- heatmap$tree_row[["labels"]][heatmap$tree_row[["order"]]]
        ind = 1
        for (i in 1:length(unique(gene.groups))) {
            which.group <- gene.groups[which(names(gene.groups) ==
                gene.tree.order[ind])]
            names <- names(which(gene.groups == which.group))
            report.gene.groups[which(all.genes \%in\% names)] <- alph[i]
            ind <- ind + length(names)
        }
        report.gene.groups <- data.frame(Gene.Groups = report.gene.groups,
            stringsAsFactors = T)
        rownames(report.gene.groups) <- all.genes
        if (is.null(GeneGroup_Name) == FALSE) {
            colnames(report.gene.groups) <- GeneGroup_Name
        }
    }
    if (extractSamples == TRUE) {
        report.sample.groups <- rep("_", length(all.samples))
        sample.groups <- cutree(heatmap$tree_col, k = num_Sample.groups)
        alph <- LETTERS[1:num_Sample.groups]
        sample.tree.order <- heatmap$tree_col[["labels"]][heatmap$tree_col[["order"]]]
        ind = 1
        for (i in 1:length(unique(sample.groups))) {
            which.group <- sample.groups[which(names(sample.groups) ==
                sample.tree.order[ind])]
            names <- names(which(sample.groups == which.group))
            report.sample.groups[which(all.samples \%in\% names)] <- alph[i]
            ind <- ind + length(names)
        }
        report.sample.groups <- data.frame(Sample.Groups = report.sample.groups,
            stringsAsFactors = T)
        rownames(report.sample.groups) <- all.samples
        if (is.null(SampleGroup_Name) == FALSE) {
            colnames(report.sample.groups) <- SampleGroup_Name
        }
    }
    if (extractGenes == TRUE & extractSamples == TRUE) {
        return(c(list(Gene_Groups = report.gene.groups, Sample_Groups = report.sample.groups)))
    }
    if (extractGenes == TRUE & extractSamples == FALSE) {
        return(report.gene.groups)
    }
    if (extractGenes == FALSE & extractSamples == TRUE) {
        return(report.sample.groups)
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
