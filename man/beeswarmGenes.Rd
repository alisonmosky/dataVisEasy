\name{beeswarmGenes}
\alias{beeswarmGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a Beeswarm Plot for a Set of Genes
}
\description{
Create highly customizable beeswarm plots for a given set of genes. Includes options to color for and group by various annotations as well as to facet the genes.
}
\usage{
beeswarmGenes(data, list, exact = TRUE, is.raw.Ct = FALSE, na.fix = 2,
  squishy = FALSE, color.by = "blue", custom.color.vec = FALSE, groupby.x = NULL,
  custom.group.vec = FALSE, facet.wrap = FALSE, ncols = 2, scales = "free_y",
  legend.position = "default", axis.text.x.size = 25, point.size = 3, transparency = 1,
  percent.mad = 0.5, dodge.width =0.8, return.ggplot.input = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows
}
  \item{list}{
character vector of genes/variables to be pulled out of the data matrix for viewing
}
  \item{exact}{
whether or not to search for exact or inexact matches of 'list' in 'data'. If exact = T (default) heatmap will plot genes/variables that exactly match the list supplied. If set to FALSE, will search for inexact matches.
}
  \item{is.raw.Ct}{
logical. If set to TRUE, will reverse the scale of the data to indicate low values as high expression as in the case of raw Ct values from qPCR, in this case, missing values will also be set to a high value to reflect low expression level
}
  \item{na.fix}{
option to treat missing/NA values as an offset from the minimum value. Ex a value of 2 will set missing values to min(data) - 2. If coloring by a specific gene, will still be colored black. When na.fix is set to FALSE missing values will be removed
}
  \item{squishy}{
similar to setting limits for the plot but all values outside the range will be set to the minimum or maximum value depending on whether or not it falls below or above the range respectively.
}
  \item{color.by}{
How the points are colored. There are several different options that this argument can take. If set to a single color (the default, blue), all points will be colored by that color. A gene name (must be present in the rownames of the input data) can be supplied where the points will be colored according to the expression level of the indicated gene, see \link[dataVisEasy]{myColorRamp5}. The name of an annotation can be provided that must match the colnames of the annotations dataframe stored in the params list object. If the colors of this annotation are also specified in annot_cols, also stored in the params list object, those colors will be used for the indicated levels of the annotation. If the colors are not provided, default colors will be used.
}
  \item{custom.color.vec}{
option to provide a custom color vector not linked to annotations or gene expression level. In this case, the order of the colors should correspond to the order of the samples/columns in the input data.
}
  \item{groupby.x}{
If and how the points should be grouped. Teh default is set to NULL. When color.by is provided and linked to the annotations dataframe stored in the params list object, groupby.x will be grouped by the same annotations as color.by.  If separating the points by a group is not desired, setting this parameter to FALSE will override this default. Likewise, if groupby.x = NULL (the default) and color.by is not found in the colnames of the annotations, groupby.x will be set to FALSE. Additionally, the name of an annotation can be provided that must match the colnames of the annotations dataframe stored in the params list object. The points will then be grouped by this annotation even if they are colored by a separate annotation.
}
  \item{custom.group.vec}{
option to provide a custom group vector not linked to annotations. In this case, the order of the indicated groupings should correspond to the order of the samples/columns in the input data.
}
  \item{facet.wrap}{
logical, should genes be faceted. When set to FALSE (the default), the x axis of the plot will be the genes in question. If set to TRUE, each gene will be plotted separately in its own facet and the x axis will be whatever color.by or groupby.x is set to if provided.
}
  \item{ncols}{
if facet.wrap = TRUE, how many columns in the facet. See \link[ggplot2]{facet_wrap}
}
  \item{scales}{
freedom of scales if facet.wrap = TRUE. See \link[ggplot2]{facet_wrap}
}

  \item{legend.position}{
should the legend be shown and if so where should it be placed. If left as default, legend will be drawn at the right if colored by an annotation, otherwise will not be drawn, can be overwritten by setting the legend position to one of "top", "right","left","bottom","none".
}
  \item{axis.text.x.size}{
size of axis text on the x axis
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

  \item{dodge.width}{
  the dodge.width of the grouped points, see \link[ggbeeswarm]{geom_beeswarm}
}
  \item{return.ggplot.input}{
logical. If true, will return the input dataframe to the ggplot object as well as the code to reproduce. Useful if more customization is required.
}
}
\details{
A beeswarm plot will be generated from the input data, subset for the list of genes provided. Options to color and group the plotted points are most easily acheived through use of the annotations dataframe stored in the params list object (although this is not necessary and custom color and group vectors can be provided as well as coloring by a single color or by expression level). See \link[dataVisEasy]{params}, \link[dataVisEasy]{set_annotations}, and \link[dataVisEasy]{set_annot_cols} for more information on setting up annotations.
}
\value{
A ggplot object. Additional layers can be added to the returned ggplot object to further customize theme and aesthetics.

If return.ggplot.input is set to TRUE, will return a list with the dataframe, coloring and call to ggplot for plotting.
\item{'input_data'}{the dataframe used for plotting which will contain the expression levels of the chosen genes (if values are squished to fit the plot the values will be similarly squished) as well as the annotations if available.}
\item{'coloring'}{the coloring parameters used in the plot the plot call will refer to items in this object and may need to be saved separately}
\item{plot_call}{the call to ggplot that generated the plot. Note that simply accessing it by $plot_call will include escape characters. The full call can be accessed by cat(plot$plot_call). Please note that many of the parameters (those in lowercase) in the call are input parameters to the original function and must be input to properly recreate the plot.}


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
For more information on customizing the returned ggplot object, please see ggplot2 helpfiles, specifically those related to setting the theme.

See \link[dataVisEasy]{params}, \link[dataVisEasy]{set_annotations}, and \link[dataVisEasy]{set_annot_cols} for more information on setting up annotations and associated colors.

}
\examples{

initiate_params()
data(RAGP_norm)
data(RAGP_annots)
set_annotations(RAGP_annots)

gene.list <- c("Npy","Sst","Gal","Th","NeuN","")
beeswarmGenes(RAGP_norm, gene.list)
beeswarmGenes(RAGP_norm, gene.list, color.by = "purple")
beeswarmGenes(RAGP_rawCt, gene.list, is.raw.Ct = TRUE)

##Color by an annotation
beeswarmGenes(RAGP_norm, gene.list, color.by = "State")
beeswarmGenes(RAGP_norm, gene.list, color.by = "State", facet.wrap = TRUE) ##Facet, separate plots

##Color by one annotation, group by another annotation or no groupings
beeswarmGenes(RAGP_norm, gene.list, color.by = "State", groupby.x = "Connectivity",
  facet.wrap = TRUE)
beeswarmGenes(RAGP_norm, gene.list, color.by = "State", groupby.x = FALSE)

##limit scale
###add xlimits onto ggplot object, points outside limits will be removed
beeswarmGenes(RAGP_norm, gene.list, color.by = "State", groupby.x = FALSE) + ylim(c(-5,5))
###point outside limits will be set to the max/min
beeswarmGenes(RAGP_norm, gene.list, color.by = "State", groupby.x = FALSE, squishy = c(-5,5))


##change other asthetics within the function
beeswarmGenes(RAGP_norm, gene.list, color.by = "Connectivity", groupby.x = FALSE,
point.size = 5, transparency = 0.5)


##Add layers onto ggplot object
###the function returns a ggplot object, therefore aesthetics can be added with additional layers
beeswarmGenes(RAGP_norm, gene.list, color.by = "Connectivity", groupby.x = FALSE) +
    theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    legend.position = "bottom") + ggtitle("Beeswarm of Select Genes")
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
