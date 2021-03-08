\name{DensityGenes}
\alias{DensityGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Density Plots for a Set of Genes
}
\description{
Create highly customizable density plots for a given set of genes. Includes options to color for and group by various annotations.
}
\usage{
DensityGenes(data, list, color.by = "blue", exact = TRUE, facet.annotation = FALSE,
  is.raw.Ct = FALSE, na.fix = 2, transparency = 0.5, ncols = 2, scales = "free",
  legend.position = "default", return.ggplot.input = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows
}
  \item{list}{
character vector of genes/variables to be pulled out of the data matrix for viewing
}
  \item{color.by}{
How the plots are colored and grouped. There are several different options that this argument can take. If set to a single color (the default, blue), one density plot will be shown in the indicated color. The name of an annotation can be provided that must match the colnames of the annotations dataframe stored in the params list object. If the colors of this annotation are also specified in annot_cols, also stored in the params list object, those colors will be used for the indicated levels of the annotation. If the colors are not provided, default colors will be used. In this case, multiple overlapping density plots will be shown separating the annotations. Useful for seeing if the range of data differs between different groups of samples.
}
  \item{exact}{
whether or not to search for exact or inexact matches of 'list' in 'data'. If exact = T (default) heatmap will plot genes/variables that exactly match the list supplied. If set to FALSE, will search for inexact matches.
}
  \item{facet.annotation}{
By default, when multiple genes are plotted the function will facet to separate the genes. If facet.annotation is set to "wrap" or "grid", will also separate by the annotation supplied for coloring. See \link[ggplot2]{facet_wrap} and \link[ggplot2]{facet_grid} for more information on faceting.
}
  \item{is.raw.Ct}{
logical. If set to TRUE, will reverse the scale of the data to indicate low values as high expression as in the case of raw Ct values from qPCR, in this case, missing values will also be set to a high value to reflect low expression level
}
  \item{na.fix}{
option to treat missing/NA values as an offset from the minimum value. Ex a value of 2 will set missing values to min(data) - 2. If coloring by a specific gene, will still be colored black. In na.fix=F missing values will be removed
}
  \item{transparency}{
transparency or alpha value of the fill for the density plots
}
  \item{ncols}{
how many columns in the facet. See \link[ggplot2]{facet_wrap}
}
  \item{scales}{
freedom of scales in facets. See \link[ggplot2]{facet_wrap}
}
  \item{legend.position}{
should the legend be shown and if so where should it be placed. If left as default, legend will be drawn at the right if colored by an annotation, otherwise will not be drawn, can be overwritten by setting the legend position to one of "top", "right","left","bottom","none".
}
  \item{return.ggplot.input}{
logical. If true, will return the input dataframe to the ggplot object. Useful if more customization is required.
}
}
\details{
A density plot will be generated from the input data for each gene provided. Options to group and color the plots according to sample groups or annotations are acheived through use of the annotations dataframe stored in the params list object. See \link[dataVisEasy]{params}, \link[dataVisEasy]{set_annotations}, and \link[dataVisEasy]{set_annot_cols} for more information on setting up annotations.
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
##initiate_parameters
initiate_params()

##View density plots for selection of genes
DensityGenes(RAGP_norm, "Sst", exact = FALSE)  ##Sst and its receptors
DensityGenes(RAGP_rawCt, "Sst", exact = FALSE, is.raw.Ct = TRUE) #in raw Ct

##set up annotations and view distribution by annotation
set_annotations(RAGP_annots)
DensityGenes(RAGP_norm, "Sst", exact = FALSE, color.by = "State", ncols = 1)
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
