\name{myHeatmap}
\alias{myHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a customizable Heatmap of Expression Data
}
\description{
Creates a heatmap of supplied expression data with options to subset for specific genes/variables and with a variety of customizable options for aesthetics. Draws from params list object for parameters such as the range of the plotted data, the colors used, and annotations to be displayed if desired.
}
\usage{
myHeatmap(data, list = NULL, exact = TRUE, method = "pearson", linkage = "complete",
  NA.handling = "pairwise.complete.obs", clust.rows = TRUE, clust.cols = TRUE,
  scale.rows = FALSE, row.groups = NA, col.groups = NA, gaps.row = NULL, gaps.col = NULL,
  gap.width = 1, main = NULL, order.by.gene = NULL, order.by.sample = NULL,
  cell.width = NA, cell.height = NA, fontsize.row = 10, fontsize.col = 10,
  show.rownames = TRUE, show.colnames = FALSE, treeheight.row = 20, treeheight.col = 20,
  hide.plot = FALSE, na.fix = FALSE, na.offset = 2, show.legend = TRUE,
  show.annotations = TRUE, is.raw.Ct = FALSE, drop.annot.levels = TRUE,
  annotation.names.row = TRUE, annotation.names.col = TRUE, border.color = NA,
  na.color = "grey90")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows
}
  \item{list}{
character vector of genes/variables to be pulled out of the data matrix for viewing. If left as NULL, all the rownames of the supplied data matrix will be plotted
}
  \item{exact}{
whether or not to search for exact or inexact matches of 'list' in 'data'. If exact = T (default) heatmap will plot genes/variables that exactly match the list supplied. If set to FALSE, will search for inexact matches.
}
  \item{method}{
method to be used to calculate distance matrix for clustering. Accepts values to be passed to cor() such as "pearson" (default), "spearman", and "kendall" which well then be coerced to a distance matrix or any options accepted by dist()
}
  \item{linkage}{
linkage clustering method used for clustering and to be passed to hclust().  Accepts all methods accepted by hclust()
}
  \item{NA.handling}{
how missing values should be handled in the case of correlations, passed to the "use" argument of cor()
}
  \item{clust.rows}{
should rows be clustered, default = TRUE
}
  \item{clust.cols}{
should columns be clustered, default = TRUE
}
  \item{scale.rows}{
whether or not the rows should be scaled. Default is not to scale rows with options to median center, "median", or zscore, "zscore". Regardless of method used, the range of the heatmap will still be linked to the range indicated in the params list object.
}
  \item{row.groups}{
numeric to be passed to cutree(). Will split the dendrogram into the number of groups indicated
}
  \item{col.groups}{
same as row.groups but for columns
}
  \item{gaps.row}{
numeric vector specifying gaps to be inserted into rows of the data, only works if clust.rows = FALSE
}
  \item{gaps.col}{
numeric vector specifying gaps to be inserted into columns of the data, only works if clust.cols = FALSE
}
  \item{gap.width}{
numeric indicating how wide the spaces indicated in gaps.row and gaps.col should be
}
  \item{main}{
title of heatmap. Default will display "Genes of Interest:" followed by the genes supplied in the list argument with the clustering method and linkage displayed underneath.
}
  \item{order.by.gene}{
optional character equal to one of the rownames of the data to order the columns of the data by increasing levels of indicated row
}
  \item{order.by.sample}{
optional character equal to one of the colnames of the data to order the rows of the data by increasing levels of indicated column
}
  \item{cell.width}{
individual cell width in points. If left as NA, then the values depend on the size of plotting window.
}
  \item{cell.height}{
individual cell height in points. If left as NA, then the values depend on the size of plotting window.
}
  \item{fontsize.row}{
size of font for row names
}
  \item{fontsize.col}{
size of font for column names
}
  \item{show.rownames}{
logical value determining if rownames should be displayed, default = TRUE
}
  \item{show.colnames}{
logical value determining if colnames should be displayed, default = FALSE
}
  \item{treeheight.row}{
the height of a dendrogram tree for rows, if these are clustered. Default value 20 points.
}
  \item{treeheight.col}{
the height of a dendrogram tree for columns if these are clustered. Default value 20 points.
}
  \item{hide.plot}{
should the plot be displayed
}
  \item{na.fix}{
logical: should missing values be treated as NA or be set to a low value (see na.offset). Values will still be colored gray in heatmap but may aid in clustering when many missing values are present.
}
  \item{na.offset}{
option to treat missing/NA values as an offset from the minimum value. Ex a value of 2 will set missing values to min(data) - 2. Values will still be colored gray in heatmap but may aid in clustering when many missing values are present.
}
  \item{show.legend}{
logical, should legend be shown
}
  \item{show.annotations}{
logical, should annotation legend be shown
}
  \item{is.raw.Ct}{
logical. If set to TRUE, will reverse the scale of the data to indicate low values as high expression as in the case of raw Ct values from qPCR, in this case, missing values will also be set to a high value to reflect low expression level
}
  \item{drop.annot.levels}{
logial, should annotations not included in the output heatmap be shown in the annotation legend.
}
  \item{annotation.names.row}{
logial value showing if the names for row annotation tracks should be drawn
}
  \item{annotation.names.col}{
logial value showing if the names for column annotation tracks should be drawn
}
  \item{border.color}{
color for borders, default set to NA
}
  \item{na.color}{
color for cells with NA values
}
}
\details{
This function utilizes pheatmap to display a heatmap of supplied data. User can specify the genes/variables to be displayed in the heatmap which will be subset within the function itself. All samples/observations supplied will be plotted. myHeatmap() draws many of its parameters from the params object. params$scale.range indiates the range of data to be displayed in the heatmap. Values outside this range will be squished to fit this range after clustering where values above will be indicated with the highest expression levels and values below will be indicated with the lowest expression level. params$scale.colors holds the colors used in the heatmap and params$n.colors.range is a numeric value indicating how many diffrent colors should be allowed. If column and row annotations are supplied in params$annot_samps and params$annot_genes, respectivey, the annotations will be indicated along the top and left side of the heatmap respectively. If the colors corresponding to the levels of these annotations are specified in params$annot_cols, the inicated colors will be used, if not supplied, defualt colors will be assigned to show the annotation levels. For more information on the setting of parameters found in the params list object, see params and the associated functions to set and update these values with set_'parameter' where 'parameter' is equal to any of the above parameters listed above.A subtitle for the heatmap indicates the method for clustering/creation of a distance matrix and the linkage method used for clustering.
}
\value{
a pheatmap object
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
See \link[dataVisEasy]{params}, \link[dataVisEasy]{set_scale.range}, \link[dataVisEasy]{set_scale.colors}, \link[dataVisEasy]{set_n.colors.range},\link[dataVisEasy]{set_annotations}, \link[dataVisEasy]{set_annotations.genes}, \link[dataVisEasy]{set_annot_samps}, \link[dataVisEasy]{set_annot_genes}, and \link[dataVisEasy]{set_annot_cols} for more information on setting up annotations and associated colors.
}
\examples{
##initiate parameters
initiate_params()

##heatmap of all data
myHeatmap(RAGP_norm)

##heatmap of selected genes (inexact matches)
myHeatmap(RAGP_norm, c("Kcn","Scn","Hcn","Cacn"), exact = FALSE, main = "Ion Channels")

##split clusters into groups
myHeatmap(RAGP_norm, row.groups = 6, col.groups = 5)

##order by gene expression
myHeatmap(RAGP_norm, order.by.gene = "NeuN")

##set width, height, other aesthetics
myHeatmap(RAGP_norm, c("Sst","Npy","Gal","Tac","Bdk"), exact = FALSE,
  treeheight.row = 25, treeheight.col = 5, cell.width = .5, cell.height = 20, fontsize.row = 20)


##add annotation tracks
set_annotations(RAGP_annots) ##sets annotations
set_annot_samps(c("Animal","Connectivity")) ##choose which annotations show up on heatmap
myHeatmap(RAGP_norm)



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
