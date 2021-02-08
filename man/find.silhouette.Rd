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
find.silhouette(data,to.sil = "samples", to.view = "rand.to.clust", ngroups,
  maxgroups = 12, max.iter = 10, method = "pearson",
  NA.handling = "pairwise.complete.obs", linkage = "complete",
  main = "Average Silhouette Width",
  axis.label = "Silhouette Width", main.label.size = 30, axis.label.size = 20,
  legend.position = "bottom")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows
}
  \item{to.sil}{
whether the clustering and comparison will be performed on the columns ("samples", default) or the rows ("genes")
}
  \item{to.view}{
there are three options for which comparisons to show. "rand.to.clust" (the default) will show a density plot comparing the clustering of the provided data for ngroups clusters to max.iter randomized trials. "all.clusts" will compare the silhouette coefficients for 2:maxgroups clusters for the provided data. "rand.all.clusts" combines the above two options, showing 2:maxgroups clusters for the provided data in addition to max.iter randomized trials for each number of clusters
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
##initiate_parameters
initiate_params()

##view heatmap of data
myHeatmap(RAGP_norm)

###use find.silhouette to help determine optimal number of clusters


##compare a given cluster number to the same cluster number for randomized data
find.silhouette(RAGP_norm, ngroups = 6, to.sil = "genes", to.view = "rand.to.clust")
find.silhouette(RAGP_norm, ngroups = 6, to.sil = "genes", to.view = "rand.to.clust", max.iter = 20)
find.silhouette(RAGP_norm, ngroups = 6, to.sil = "genes", to.view = "rand.to.clust", max.iter = 100)

##for the data provided, determine optimal number of clusters
find.silhouette(RAGP_norm, to.sil = "genes", to.view = "all.clusts", maxgroups = 12)
find.silhouette(RAGP_norm, to.sil = "genes", to.view = "all.clusts", maxgroups = 15)

##for data provided, compare multiple cluster options to randomized data
find.silhouette(RAGP_norm, to.sil = "genes", to.view = "rand.all.clusts", maxgroups = 12)
find.silhouette(RAGP_norm, to.sil = "genes", to.view = "rand.all.clusts",
  maxgroups = 9, max.iter = 100)
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
