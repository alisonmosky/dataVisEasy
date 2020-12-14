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
  PCs.to.plot = c("PC1", "PC2"), legend.position = "right", main = NULL,point.size = 5,
  transparency = 1, percent.mad = 0.5)
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A ggplot object. Within the ggplot object are additional options that allow the user to add the text of gene or sample names to the plot (see examples, these are specifically designed to be added outside the function itself for proper positioning). Additional layers can be added to the returned ggplot object to further customize theme and aesthetics.
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
function (data, to.pca = "samples", nPcs = 3, color.by = "blue",
    custom.color.vec = FALSE, PCs.to.plot = c("PC1", "PC2"),
    legend.position = "right", main = NULL, percent.mad = 0.5)
{
    if (to.pca == "samples") {
        temp.annotations <- params$annotations
        pca <- pcaMethods::pca(t(data), nPcs = nPcs)
        pca.scrs <- pcaMethods::scores(pca)
        pca.ldgs <- pcaMethods::loadings(pca)
        pca.data <- data.frame(pca.scrs, temp.annotations)
        if (color.by \%in\% rownames(data) | sum(custom.color.vec !=
            FALSE) > 0) {
            if (color.by \%in\% rownames(data)) {
                genedat <- data[which(rownames(data) == color.by),
                  ]
                cols <- myColorRamp5(params$expression_gradient.colors,
                  genedat, percent.mad = percent.mad)
            }
            else {
                cols <- custom.color.vec
            }
            p <- ggplot(pca.data, aes(x = eval(parse(text = PCs.to.plot[1])),
                y = eval(parse(text = PCs.to.plot[2])), fill = cols)) +
                geom_point(pch = 21, color = "black", size = 5) +
                scale_fill_identity() + labs(x = paste(PCs.to.plot[1]),
                y = paste(PCs.to.plot[2])) + ggtitle(main) +
                theme_bw() + theme(panel.grid = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 40),
                axis.text = element_text(size = 25), axis.title = element_text(size = 30),
                legend.position = legend.position)
        }
        else {
            if (color.by \%in\% colnames(temp.annotations)) {
                suppressWarnings(if (is.na(temp.annotations) ==
                  FALSE) {
                  if (sum(colnames(data) \%notin\% rownames(temp.annotations)) !=
                    0) {
                    stop("colnames of input data do not match rownames of annotations, cannot link annotations to data")
                  }
                  temp.annotations <- temp.annotations[match(colnames(data),
                    rownames(temp.annotations)), ]
                })
                if (color.by \%in\% names(params$annot_cols)) {
                  cols <- as.factor(pca.data[, which(colnames(pca.data) ==
                    color.by)])
                  colors <- params$annot_cols[[which(names(params$annot_cols) ==
                    color.by)]]
                }
                else {
                  cols <- as.factor(pca.data[, which(colnames(pca.data) ==
                    color.by)])
                  colors <- hue_pal()(length(levels(cols)))
                }
            }
            else {
                cols <- color.by
                colors <- color.by
            }
            p <- ggplot(pca.data, aes(x = eval(parse(text = PCs.to.plot[1])),
                y = eval(parse(text = PCs.to.plot[2])), fill = cols)) +
                geom_point(pch = 21, color = "black", size = 5) +
                scale_fill_manual(values = colors) + labs(x = paste(PCs.to.plot[1]),
                y = paste(PCs.to.plot[2])) + ggtitle(main) +
                theme_bw() + theme(panel.grid = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 40),
                axis.text = element_text(size = 25), axis.title = element_text(size = 30),
                legend.position = legend.position)
        }
    }
    if (to.pca == "genes") {
        temp.annotations.genes <- params$annotations.genes
        pca <- pcaMethods::pca((data), nPcs = nPcs)
        pca.scrs <- pcaMethods::scores(pca)
        pca.ldgs <- pcaMethods::loadings(pca)
        pca.data <- data.frame(pca.scrs, temp.annotations.genes)
        if (sum(custom.color.vec != FALSE) > 0) {
            cols <- custom.color.vec
            p <- ggplot(pca.data, aes(x = eval(parse(text = PCs.to.plot[1])),
                y = eval(parse(text = PCs.to.plot[2])), fill = cols)) +
                geom_point(pch = 21, color = "black", size = 5) +
                scale_fill_identity() + labs(x = paste(PCs.to.plot[1]),
                y = paste(PCs.to.plot[2])) + ggtitle(main) +
                theme_bw() + theme(panel.grid = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 40),
                axis.text = element_text(size = 25), axis.title = element_text(size = 30),
                legend.position = legend.position)
        }
        else {
            if (color.by \%in\% colnames(temp.annotations.genes)) {
                suppressWarnings(if (is.na(temp.annotations.genes) ==
                  FALSE) {
                  if (sum(rownames(data) \%notin\% rownames(temp.annotations.genes)) !=
                    0) {
                    stop("rownames of input data do not match rownames of annotations, cannot link annotations to data")
                  }
                  temp.annotations.genes <- temp.annotations.genes[match(rownames(data),
                    rownames(temp.annotations.genes)), ]
                })
                if (color.by \%in\% names(params$annot_cols)) {
                  cols <- as.factor(pca.data[, which(colnames(pca.data) ==
                    color.by)])
                  colors <- params$annot_cols[[which(names(params$annot_cols) ==
                    color.by)]]
                }
                else {
                  cols <- as.factor(pca.data[, which(colnames(pca.data) ==
                    color.by)])
                  colors <- hue_pal()(length(levels(cols)))
                }
            }
            else {
                cols <- color.by
                colors <- color.by
            }
            p <- ggplot(pca.data, aes(x = eval(parse(text = PCs.to.plot[1])),
                y = eval(parse(text = PCs.to.plot[2])), fill = cols)) +
                geom_point(pch = 21, color = "black", size = 5) +
                scale_fill_manual(values = colors) + labs(x = paste(PCs.to.plot[1]),
                y = paste(PCs.to.plot[2])) + ggtitle(main) +
                theme_bw() + theme(panel.grid = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 40),
                axis.text = element_text(size = 25), axis.title = element_text(size = 30),
                legend.position = legend.position)
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
