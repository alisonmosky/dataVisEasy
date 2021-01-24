\name{find.silhouette}
\alias{find.silhouette}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find and Compare the Silhouette Coefficient for Clustering Assignments
}
\description{
Compare the silhouette coefficient for a given number of clusters against other cluster numbers to determine the optimal number of clustes and against randomized trials to determine how the silhouette coefficients of the identified clusters compare to those of randomized data.
}
\usage{
find.silhouette(data, ngroups, maxgroups = 12, max.iter = 10, method = "pearson",
  NA.handling = "pairwise.complete.obs", linkage = "complete", to.sil = "samples",
  to.view = "rand.to.clust", main = "Average Silhouette Width",
  axis.label = "Silhouette Width", main.label.size = 30, axis.label.size = 20,
  legend.position = "bottom")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows
}
  \item{ngroups}{
number of clusters to compare (when to.view is set to "rand.to.clust")
}
  \item{maxgroups}{
maximum number of clusters to compare (when to.view is set to "all.clust" or "rand.all.clusts"). Will test all numbers of clusters from 2:maxgroups
}
  \item{max.iter}{
number of randomized trials to compare against (when to.view is set to "rand.to.clust" or "rand.all.clusts")
}
  \item{method}{
method to be used to calculate distance matrix for clustering. Accepts values to be passed to cor() such as "pearson" (default), "spearman", and "kendall" which well then be coerced to a distance matrix or any options accepted by dist()
}
  \item{NA.handling}{
how missing values should be handled in the case of correlations, passed to the "use" argument of cor()
}
  \item{linkage}{
linkage clustering method used for clustering and to be passed to hclust().  Accepts all methods accepted by hclust()
}
  \item{to.sil}{
whether the clustering and comparison will be performed on the columns ("samples", default) or the rows ("genes")
}
  \item{to.view}{
there are three options for which comparisons to show. "rand.to.clust" (the default) will show a density plot comparing the clustering of the provided data for ngroups clusters to max.iter randomized trials. "all.clusts" will compare the silhouette coefficients for 2:maxgroups clusters for the provided data. "rand.all.clusts" combines the above two options, showing 2:maxgroups clusters for the provided data in addition to max.iter randomized trials for each number of clusters
}
  \item{main}{
title for generated plot, default is "Average Silhouette Width"
}
  \item{axis.label}{
axis label for the silhouette widt, default is "Silhouette Width"
}
  \item{main.label.size}{
size of plot title
}
  \item{axis.label.size}{
size of axis title
}
  \item{legend.position}{
position of legend when to.view is set to "rand.all.clusts"
}
}
\details{
Assesses the sihouette coefficient of clustered data using correlation or distance metrics to be passed to hclust which is then split based on cutree, corresponds to same specifications of a heatmap with the same inputs using \link[dataVisEasy]{myHeatmap} or \link[dataVisEasy]{myHeatmapByAnnotation}.
}
\value{
A ggplot object. Additional layers can be added to the returned ggplot object to further customize theme and aesthetics.
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
function (data, ngroups, maxgroups = 12, max.iter = 10, method = "pearson",
    NA.handling = "pairwise.complete.obs", linkage = "complete",
    to.sil = "samples", to.view = "rand.to.clust", main = "Average Silhouette Width",
    axis.label = "Silhouette Width", main.label.size = 30, axis.label.size = 20,
    legend.position = "bottom")
{
    if (method \%in\% c("spearman", "pearson", "kendall")) {
        clust.genes <- (as.dist(1 - cor(t(data), method = method,
            use = NA.handling)))
        clust.samps <- (as.dist(1 - cor(data, method = method,
            use = NA.handling)))
    }
    if (method \%in\% c("euclidean", "maximum", "manhattan", "canberra",
        "binary", "minkowski")) {
        clust.genes <- dist(data, method = method)
        clust.samps <- dist(t(data), method = method)
    }
    if (to.sil == "samples") {
        clust.mat <- clust.samps
    }
    if (to.sil == "genes") {
        clust.mat <- clust.genes
    }
    clusts <- hclust(clust.mat, method = linkage)
    if (to.view == "rand.to.clust") {
        groupings <- cutree(clusts, k = ngroups)
        sil.coef <- summary(silhouette(groupings, clust.mat))$avg.width
        rands <- replicate(max.iter, summary(silhouette(permute(groupings),
            clust.mat))$avg.width, simplify = "vector")
        data.to.plot <- data.frame(vals = c(rands, sil.coef))
        p <- ggplot(data.to.plot, aes(x = vals)) + geom_density(size = 1) +
            geom_point(x = sil.coef, y =
            density(data.to.plot$vals)$y[which.min(abs(density(data.to.plot$vals)$x -
                sil.coef))], size = 3) + geom_text(x = sil.coef,
            y = (density(data.to.plot$vals)$y[which.min(abs(density(data.to.plot$vals)$x -
                sil.coef))]) + 5, label = round(sil.coef, 3)) +
            geom_segment(x = sil.coef, xend = sil.coef, y = 0,
                yend = density(data.to.plot$vals)$y[which.min(abs(density(data.to.plot$vals)$x -
                  sil.coef))], size = 1) + theme_bw() + theme(panel.grid = element_blank(),
            axis.title = element_text(size = axis.label.size),
            plot.title = element_text(size = main.label.size,
                hjust = 0.5)) + xlab(paste(axis.label)) + ylab("Density") +
            ggtitle(paste(main))
    }
    if (to.view != "rand.to.clust") {
        sil.coefs <- c()
        for (nclust in 2:maxgroups) {
            groupings <- cutree(clusts, k = nclust)
            sil.coefs <- c(sil.coefs, summary(silhouette(groupings,
                clust.mat))$avg.width)
        }
        if (to.view == "all.clusts") {
            to.plot <- data.frame(clusts = 2:maxgroups, coefs = sil.coefs,
                type = "Original Data")
            p <- ggplot(to.plot, aes(x = clusts, y = sil.coefs)) +
                geom_line(size = 1) + geom_point(size = 3) +
                theme_bw() + theme(panel.grid = element_blank()) +
                xlab("Clusters") + ylab(paste(axis.label)) +
                ggtitle(paste(main)) + theme(axis.title = element_text(size = axis.label.size),
                plot.title = element_text(size = main.label.size,
                  hjust = 0.5)) + scale_x_continuous(breaks = scales::pretty_breaks())
        }
        if (to.view == "rand.all.clusts") {
            rands <- NULL
            for (nclust in 2:maxgroups) {
                groupings <- cutree(clusts, k = nclust)
                rands <- cbind(rands, replicate(max.iter, summary(silhouette(permute(groupings),
                  clust.mat))$avg.width, simplify = "vector"))
            }
            to.plot <- data.frame(clusts = 2:maxgroups, coefs = sil.coefs,
                type = "Original Data")
            colnames(rands) <- 2:maxgroups
            rands.melt <- melt(rands)
            rands.melt$type <- "Randomized Data"
            p <- ggplot() + geom_line(data = to.plot, aes(x = clusts,
                y = sil.coefs), size = 1) + geom_point(data = to.plot,
                aes(x = clusts, y = sil.coefs, shape = type),
                size = 3) + geom_boxplot(data = rands.melt, aes(x = Var2,
                y = value, group = Var2)) + geom_point(data = rands.melt,
                aes(x = Var2, y = value, group = Var2, shape = type),
                size = 1) + stat_summary(data = rands.melt, aes(x = Var2,
                y = value), fun = median, geom = "line") + stat_summary(data = rands.melt,
                aes(x = Var2, y = value), fun = median, geom = "point",
                size = 4, shape = 17) + theme_bw() + theme(panel.grid = element_blank()) +
                xlab("Clusters") + scale_shape_manual(values = c(19,
                17), labels = c("Original Data", "Randomized Data")) +
                ylab(paste(axis.label)) + ggtitle(paste(main)) +
                theme(axis.title = element_text(size = axis.label.size),
                  plot.title = element_text(size = main.label.size,
                    hjust = 0.5), legend.position = legend.position,
                  legend.direction = "horizontal", legend.title = element_blank(),
                  legend.background = element_rect(color = "black")) +
                scale_x_continuous(breaks = scales::pretty_breaks()) +
                ylim(c(min(rands), max(sil.coefs)))
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
