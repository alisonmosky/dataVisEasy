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
DensityGenes(data, list, color.by = "blue", exact = TRUE, is.raw.Ct = FALSE,
  na.fix = 2, transparency = 0.5, ncols = 2, scales = "free", legend.position = "right",
  return.ggplot.input = FALSE)
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
should the legend be shown and if so where should it be placed. No legend will be shown if colored by gene or a custom color vector is supplied.
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

If return.ggplot.input is set to TRUE, will return the dataframe used for plotting which will contain the expression levels of the chosen genes (if values are squished to fit the plot the values will be similarly squished) as well as the annotations if available.
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
For more information on customizing the returned ggplot object, please see ggplot2 helpfiles, specifically those related to setting the theme.

See \link[dataVisEasy]{params}, \link[dataVisEasy]{set_annotations}, and \link[dataVisEasy]{set_annot_cols} for more information on setting up annotations and associated colors.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, list, exact = TRUE, color.by = "blue", is.raw.Ct = FALSE,
    na.fix = 2, transparency = 0.5, ncols = 2, scales = "free",
    legend.position = "right")
{
    if (exact == TRUE) {
        dat <- data[which(rownames(data) \%in\% list), ]
        if (length(dat) == 0) {
            stop("exact matches for list not found in rownames data")
        }
        if (is.raw.Ct == F & na.fix != F) {
            dat[which(is.na(dat))] <- (min(dat, na.rm = T) -
                na.fix)
        }
        if (is.raw.Ct == T & na.fix != F) {
            dat[which(is.na(dat))] <- (max(dat, na.rm = T) +
                na.fix)
        }
    }
    if (exact == FALSE) {
        dat <- data[grep(paste(list, collapse = "|"), rownames(data)),
            ]
        if (length(dat) == 0) {
            stop("inexact matches for list not found in rownames data")
        }
        if (is.raw.Ct == F & na.fix != F) {
            dat[which(is.na(dat))] <- (min(dat, na.rm = T) -
                na.fix)
        }
        if (is.raw.Ct == T & na.fix != F) {
            dat[which(is.na(dat))] <- (max(dat, na.rm = T) +
                na.fix)
        }
    }
    temp.annotations <- params$annotations
    if (color.by \%in\% colnames(temp.annotations)) {
        if (sum(colnames(dat) \%notin\% rownames(temp.annotations)) !=
            0) {
            stop("colnames of input data do not match rownames of annotations, cannot link annotations to data")
        }
        temp.annotations <- temp.annotations[match(colnames(dat),
            rownames(temp.annotations)), , drop = FALSE]
        dat.to.plot <- data.frame(t(dat))
        dat.to.plot <- cbind(dat.to.plot, temp.annotations)
        dat.to.plot <- melt(dat.to.plot, id.vars = colnames(temp.annotations))
        if (color.by \%in\% names(params$annot_cols)) {
            cols <- as.factor(dat.to.plot[, which(colnames(dat.to.plot) ==
                color.by)])
            colors <- params$annot_cols[[which(names(params$annot_cols) ==
                color.by)]]
        }
        else {
            cols <- as.factor(dat.to.plot[, which(colnames(dat.to.plot) ==
                color.by)])
            colors <- hue_pal()(length(levels(cols)))
        }
        if (is.raw.Ct == T) {
            p <- ggplot(dat.to.plot, aes(x = value, fill = cols,
                group = eval(parse(text = color.by)))) + geom_density(alpha = transparency) +
                facet_wrap(~variable, ncol = ncols, scales = scales) +
                scale_fill_manual(values = colors) + theme_bw() +
                theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5,
                  size = 40), strip.text = element_text(size = 25),
                  strip.background.x = element_blank(), legend.position = legend.position,
                  axis.title = element_text(size = 20)) + xlab("Raw Ct Value") +
                ylab("Denstiy") + scale_y_reverse()
        }
        else {
            p <- ggplot(dat.to.plot, aes(x = value, fill = cols,
                group = eval(parse(text = color.by)))) + geom_density(alpha = transparency) +
                facet_wrap(~variable, ncol = ncols, scales = scales) +
                scale_fill_manual(values = colors) + theme_bw() +
                theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5,
                  size = 40), strip.text = element_text(size = 25),
                  strip.background.x = element_blank(), legend.position = legend.position,
                  axis.title = element_text(size = 20), axis.text.x = element_text(size = 15)) +
                xlab("Normalized Expression Level") + ylab("Density")
        }
    }
    else {
        cols <- color.by
        colors <- color.by
        dat.to.plot <- data.frame(t(dat))
        dat.to.plot <- cbind(dat.to.plot, temp.annotations)
        dat.to.plot <- melt(dat.to.plot, id.vars = colnames(temp.annotations))
        suppressWarnings(if (is.na(temp.annotations) == TRUE) {
            dat.to.plot <- dat.to.plot[-which(dat.to.plot$variable ==
                "temp.annotations"), ]
        })
        if (is.raw.Ct == T) {
            p <- ggplot(dat.to.plot, aes(x = value, fill = cols)) +
                geom_density(alpha = transparency) + facet_wrap(~variable,
                ncol = ncols, scales = scales) + scale_fill_manual(values = colors) +
                theme_bw() + theme(panel.grid = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 40),
                strip.text = element_text(size = 25), strip.background.x = element_blank(),
                legend.position = legend.position, axis.title = element_text(size = 20)) +
                xlab("Raw Ct Value") + ylab("Denstiy") + scale_y_reverse()
        }
        else {
            p <- ggplot(dat.to.plot, aes(x = value, fill = cols)) +
                geom_density(alpha = transparency) + facet_wrap(~variable,
                ncol = ncols, scales = scales) + scale_fill_manual(values = colors) +
                theme_bw() + theme(panel.grid = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 40),
                strip.text = element_text(size = 25), strip.background.x = element_blank(),
                legend.position = legend.position, axis.title = element_text(size = 20),
                axis.text.x = element_text(size = 15)) + xlab("Normalized Expression Level") +
                ylab("Density")
        }
    }
    return(p)
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
