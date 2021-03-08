\name{scatterGenes}
\alias{scatterGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a Scatter Plot Between Two Genes
}
\description{
Creates a scatter plot between two genes within a supplied data matrix. Many customizable options to color according to expression levels or annotations.
}
\usage{
scatterGenes(data, gene1, gene2, custom.x = FALSE, custom.y = FALSE,
  is.raw.Ct = FALSE, na.fix = 2, color.by = "blue",
  custom.color.vec = FALSE, xlimits = FALSE, ylimits = FALSE, squish1 = FALSE,
  squish2 = FALSE, point.size = 5, transparency = 1,legend.position = "default",
  percent.mad = 0.5, return.ggplot.input = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows
}
  \item{gene1}{
first gene to be plotted, must be present in rownames of input data
}
  \item{gene2}{
second gene to be plotted, must be present in rownames of input data
}
  \item{custom.x}{
custom values to be plotted on the x axis, MUST be same length and order as input data. If both custom.x and custom.y are supplied and the user wishes to color by annotation (example if tSNE coordiates are being supplied), must provide accompanying data matrix so that annotations can be properly subsetted. If wishing to use as simple scatterplot without coloring linked to annotations or expression levels, set data = NULL.
}
  \item{custom.y}{
custom values to be plotted on the y axis, MUST be same length and order as input data. If both custom.x and custom.y are supplied and the user wishes to color by annotation (example if tSNE coordiates are being supplied), must provide accompanying data matrix so that annotations can be properly subsetted. If wishing to use as simple scatterplot without coloring linked to annotations or expression levels, set data = NULL.
}
  \item{is.raw.Ct}{
logical: If set to TRUE, will reverse the scale of the data to indicate low values as high expression as in the case of raw Ct values from qPCR, in this case, missing values will also be set to a high value to reflect low expression level.
}
  \item{na.fix}{
option to treat missing/NA values as an offset from the minimum value. Ex a value of 2 will set missing values to min(data) - 2. If coloring by a specific gene, will still be colored black. In na.fix=F missing values will be removed
}
  \item{color.by}{
How the points are colored. There are several different options that this argument can take. If set to a single color (the default, blue), all points will be colored by that color. A gene name (must be present in the rownames of the input data but need not be the genes being plotted) can be supplied where the points will be colored according to the expression level of the indicated gene, see \link[dataVisEasy]{myColorRamp5}. The name of an annotation can be provided that must match the colnames of the annotations dataframe stored in the params list object. If the colors of this annotation are also specified in annot_cols, also stored in the params list object, those colors will be used for the indicated levels of the annotation. If the colors are not provided, default colors will be used.
}
  \item{custom.color.vec}{
option to provide a custom color vector not linked to annotations or gene expression level. In this case, the order of the colors should correspond to the order of the samples/columns in the input data.
}
  \item{xlimits}{
FALSE or numerical vector of length 2. Default of FALSE will allow limits to be set automatically based on the data. Supply desired limits on the x axis to override.
}
  \item{ylimits}{
FALSE or numerical vector of length 2. Default of FALSE will allow limits to be set automatically based on the data. Supply desired limits on the y axis to override.
}
  \item{squish1}{
FALSE or numerical vector of length 2. Should data corresponding to gene1 be limited to a specific range. In this case, values above and below the specified range will be set to the maximum and minimum respectively. Distinct from setting x or y limits as setting the limits will remove all points outside the specfied range. The squish option will restrict the range of the data to the specified range and will set the limits accordingly.
}
  \item{squish2}{
FALSE or numerical vector of length 2. Should data corresponding to gene2 be limited to a specific range. In this case, values above and below the specified range will be set to the maximum and minimum respectively. Distinct from setting x or y limits as setting the limits will remove all points outside the specfied range. The squish option will restrict the range of the data to the specified range and will set the limits accordingly.
}
  \item{point.size}{
size of points to be plotted
}
  \item{transparency}{
transparency or alpha value of the points
}
  \item{legend.position}{
should the legend be shown and if so where should it be placed. If left as default, legend will be drawn at the right if colored by an annotation, otherwise will not be drawn, can be overwritten by setting the legend position to one of "top", "right","left","bottom","none".
}
  \item{percent.mad}{
if coloring points by expression level. Passed to \link[dataVisEasy]{myColorRamp5} to determine how the data is binned
}
  \item{return.ggplot.input}{
logical. If true, will return the input dataframe to the ggplot object. Useful if more customization is required.
}
}
\details{
A scatter plot will be generated from the input data for the two genes provided. Options to color the plotted points are most easily acheived through use of the annotations dataframe stored in the params list object (although this is not necessary and custom color can be provided as well as coloring by a single color or by expression level). See \link[dataVisEasy]{params}, \link[dataVisEasy]{set_annotations}, and \link[dataVisEasy]{set_annot_cols} for more information on setting up annotations.
}
\value{
A ggplot object. Additional layers can be added to the returned ggplot object to further customize theme and aesthetics.

If return.ggplot.input is set to TRUE, will return a list with the dataframe, coloring and call to ggplot for plotting.
\item{'input_data'}{the dataframe used for plotting which will contain the expression levels of the chosen genes (if values are squished to fit the plot the values will be similarly squished) as well as the annotations if available.}
\item{'coloring'}{the coloring parameters used in the plot}
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
##initiate parameters
initiate_params()

scatterGenes(RAGP_norm, "Th","Chat")
scatterGenes(RAGP_norm, "Th","Chat", color.by = "hotpink")


##color by gene expression
scatterGenes(RAGP_norm, "Th","Chat", color.by = "Th")
scatterGenes(RAGP_norm, "Th","Chat", color.by = "Npy")

##color by an annotation
set_annotations(RAGP_annots) #see set_annot_cols to specify colors
scatterGenes(RAGP_norm, "Th","Chat", color.by = "State")

##Add layers onto ggplot object
###the function returns a ggplot object, therefore aesthetics can be added with additional layers
scatterGenes(RAGP_norm, "Th","Chat", color.by = "Connectivity") +
    theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    legend.position = "bottom")

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
