\name{myPCA}
\alias{myPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Principle Component Analysis
}
\description{
A simple function that performs a PCA of the supplied data for either the genes/rows or samples/columns with options to color according to annotations.
}
\usage{
myPCA(data, to.pca = "samples", nPcs = 3, color.by = "blue", custom.color.vec = FALSE,
  shape.by = NA, PCs.to.plot = c("PC1", "PC2"), legend.position = "right", main = NULL,
  point.size = 5,transparency = 1, percent.mad = 0.5, return.ggplot.input = FALSE,
  return.loadings = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows
}
  \item{to.pca}{
character indicating what the PCA should be performed on. Accepts either "samples" (the default) or "genes"
}
  \item{nPcs}{
number of PCs to calculate, default of 3
}
  \item{color.by}{
how the points in the PCA should be colored. Can set to a color value which will color all points the indicated color. If to.pca="samples", can set to a gene/rowname where it will be colored according to expression level. Also accepts a string pointing to a column name of the annotations (for samples) or annotations.genes (for genes) stored in the params list object. If annot_cols are specified in the params list object, points will be colored accordingly, otherwise default colors will be used.
}
  \item{custom.color.vec}{
custom color vector the same length as the number of columns (for to.pca = "samples") or the number of rows (for to.pca = "genes"). Order of custom color vector should be in the same order corresponding to the data.
}
 \item{shape.by}{
Option to change the shape of the points in the PCA. Accepts a string pointing to a column name of the annotations (for samples) or annotations.genes (for genes) stored in the params list object. Allows for up to five different shapes to used. User may override default shapes, see \link[ggplot2]{scale_shape_manual}. Keep in mind that the points are set to 'fill' and may cause problems depending on custom shapes desired. Please set return.ggplot.input to TRUE to allow for further customizations.
}
  \item{PCs.to.plot}{
character vector of length 2 indicating which PCs should be displayed in the plot. Default is c("PC1","PC2"). Can specify any 2 PCs from 1:nPcs.
}
  \item{legend.position}{
position of the legend
}
  \item{main}{
plot title if desired
}
  \item{point.size}{
size of points to be plotted
}
  \item{transparency}{
transparency or alpha value of the points
}
  \item{percent.mad}{
if coloring points by expression level. Passed to \link[dataVisEasy]{myColorRamp5} to determine how the data is binned
}
  \item{return.ggplot.input}{
logical. If true, will return the input dataframe to the ggplot object as well as the code to reproduce. Useful if more customization is required.
}
  \item{return.loadings}{
logical. If true, will return the loadings from the pca.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A ggplot object. The input data frame to the ggplot object contains the scores, the sample/gene names, and annotations/gene annotations if an annotation is specified. This input object can be accessed by setting return.ggplot.input to TRUE. The inclusion of the gene and sample names allows the user to add the text of gene or sample names to the plot (see examples, these are specifically designed to be added outside the function itself for proper positioning). Additional layers can be added to the returned ggplot object to further customize theme and aesthetics.

A ggplot object. Additional layers can be added to the returned ggplot object to further customize theme and aesthetics.

If return.ggplot.input is set to TRUE, will return a list with the dataframe, coloring and call to ggplot for plotting.
\item{'input_data'}{the dataframe used for plotting which will contain the scores for each PC as well as the annotations if provided in addition to the Sample/Gene names (which can be added to the plot with geom_label or geom_plot  with aes(label =Genes/Samples) if desired )}
\item{'coloring'}{the coloring parameters used in the plot the plot call will refer to items in this object and may need to be saved separately}
\item{'PCA_variability'}{details the percent variability in each PC, plot_call will refer to this for axis titles}
\item{plot_call}{the call to ggplot that generated the plot. Note that simply accessing it by $plot_call will include escape characters. The full call can be accessed by cat(plot$plot_call). Please note that many of the parameters (those in lowercase) in the call are input parameters to the original function and must be input to properly recreate the plot.}


If return.loadings is set to TRUE, will return a matrix of the loadings from the PCA.
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
##initiate parameters
initiate_params()

##view PCA, specify PCs (default is to PCA samples, can also do genes)
myPCA(RAGP_norm)
myPCA(RAGP_norm, nPcs = 5, PCs.to.plot = c("PC1","PC3"))

##color by gene expression
myPCA(RAGP_norm, color.by = "Th")

##color by annotation
set_annotations(RAGP_annots)
myPCA(RAGP_norm, nPcs = 5,  color.by = "Animal")
myPCA(RAGP_norm, PCs.to.plot = c("PC1","PC3"), color.by = "Connectivity")

##specify colors
state.cols <- RColorBrewer::brewer.pal(6,"Set1"); names(state.cols) <- LETTERS[1:6]
annot_cols <- list('Connectivity'=c("SAN-Projecting"="blue","Non-SAN-Projecting"="violet",
                  "No Info Available"="grey"),
                   'Animal'=c("PR1534"="#0571b0","PR1643"="#ca0020","PR1705"="#92c5de",
                   "PR1729"="#f4a582"),
                   "State"=c(state.cols))
set_annot_cols(annot_cols)

myPCA(RAGP_norm, nPcs = 5,  color.by = "Animal")


##Add layers onto ggplot object
###the function returns a ggplot object, therefore aesthetics can be added with additional layers
myPCA(RAGP_norm, PCs.to.plot = c("PC1","PC3"), color.by = "Connectivity") +
  ggtitle("PCA colored by Connectivity") +
  theme(legend.position = "bottom", legend.direction = "vertical",
  axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
  plot.title = element_text(size = 25))


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
