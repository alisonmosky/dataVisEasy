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
extractClusters(data, heatmap,to.extract = "genes",nclusters, GeneGroup_Name = NULL,
  SampleGroup_Name = NULL, all.genes = NULL, all.samples = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows. Should be the same data supplied to the heatmap function
}
  \item{heatmap}{
Output of a heatmap generated from \link[dataVisEasy]{myHeatmap} or \link[dataVisEasy]{myHeatmapByAnnotation}.  The heatmap generated from these functions can be saved to a variable which should serve as the input here. IMPORTANT. The heatmap supplied here MUST have show.colnames and show.rownames = TRUE.
}
  \item{to.extract}{
whether clusters should be extracted for genes (rows), samples(columns), or both. Can accept "genes", "samples", or "both"
}
  \item{nclusters}{
numeric, how many clusters should be extracted. If to.extract is set to "both", numeric vector should be supplied listing number of gene (row clusters) first
}
  \item{GeneGroup_Name}{
optional name for the groupings to be extracted
}

  \item{SampleGroup_Name}{
optional name for the groupings to be extracted
}

  \item{all.genes}{
gene names, if annotations should be in placed in larger gene list than matrix in question. Default assumes that the rownames of the input data represent all genes in question (this is true even if the genes in the heatmap are a subset of the input data). Genes not represented in the heatmap or input data that do not have a cluster designation will be assigned "No_Annot".
}
  \item{all.samples}{
sample names, if annotations should be in placed in larger sample list than matrix in question. Default assumes that the colnames of the input data represent all samples in question. If these annotations should be placed within a larger list (for example if the samples in the input data are a subset of all samples). Samples not represented in the heatmap or input data that do not have a cluster designation will be assigned "No_Annot".
}
}
\details{
This function is helpful for assigned annotations to samples or genes based on how they cluster. The options of all.samples and all.genes allow these cluster assignments to be placed within a larger list of variables (useful if combining these annotations with other pre-existing annotations, where the vectors need to be of the same length to merge correctly).

If GeneGroup_Name and/or SampleGroup_Name are provided, these vectors can be passed directly to \link[dataVisEasy]{update_annotations} or \link[dataVisEasy]{update_annotations.genes}
}
\value{
a named vector or a list object of two named vectors (in a case where both genes and samples are extracted) where the assignments of each gene/sample are designated by letter A:num.groups. Importantly, these groupings are assigned as they appear in the heatmap itself where cluster assignment A is located on the top left (going across from left to right in the case of extractSamples and going down from top to bottom in the case of extractGenes)
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

\seealso{
See \link[dataVisEasy]{set_annotations} and \link[dataVisEasy]{update_annotations.genes} or \link[dataVisEasy]{update_annotations} and \link[dataVisEasy]{update_annotations.genes} for how the output of this function can be passed directly to these functions to set or update annotations stored in the params list object
}
\examples{
##initiate parameters
initiate_params()

##create clustered heatmap, defining how many groups clusters should split into
heatmap <- myHeatmap(RAGP_norm, row.groups = 6, col.groups = 5)

##extract clusters for rows, named from top to bottom
row.clusts <- extractClusters(RAGP_norm, heatmap = heatmap, to.extract = "genes", nclusters = 6)
head(row.clusts)

##extract clusters for columns, named from left to right
col.clusts <- extractClusters(RAGP_norm, heatmap = heatmap, to.extract = "samples", nclusters = 5)
head(col.clusts)

##extract both
all.clusts <- extractClusters(RAGP_norm, heatmap = heatmap, to.extract = "both", nclusters = c(6,5))
head(all.clusts$Gene_Groups)
head(all.clusts$Sample_Groups)



###extract clusters from a heatmap with subsetted genes, catecholaminergic/cholinergic receptors
#subset data for receptors
RAGP_norm.receptors <- subsetGenes(RAGP_norm, list = c("Chrm","Chrn","Adr"),
exact = FALSE, subset.by = "name")
heatmap.receptors <- myHeatmap(RAGP_norm.receptors,
  main = "Catecholaminergic and Cholinergic Receptors", row.groups = 6)

receptor.groups <- extractClusters(RAGP_norm.receptors, heatmap = heatmap.receptors,
  to.extract = "genes", nclusters = 6, GeneGroup_Name = "Receptor_Groups")
head(receptor.groups)

##if we want the annotations to include genes not present in the heatmap, set all.genes
###genes not present in the heatmap will be assigned "No_Annot" for no annotation
receptor.groups1 <- extractClusters(RAGP_norm.receptors, heatmap = heatmap.receptors,
  to.extract = "genes", nclusters = 6, GeneGroup_Name = "Receptor_Groups",
  all.genes = rownames(RAGP_norm))
head(receptor.groups1)

##the same result would be true if a list was supplied to the heatmap to subset
heatmap.receptors2 <- myHeatmap(RAGP_norm, list = c("Chrm","Chrn","Adr"), exact = FALSE,
  main = "Catecholaminergic and Cholinergic Receptors", row.groups = 6)
receptor.groups2 <- extractClusters(RAGP_norm, heatmap = heatmap.receptors2,
  to.extract = "genes", nclusters = 6, GeneGroup_Name = "Receptor_Groups")
head(receptor.groups2)
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
