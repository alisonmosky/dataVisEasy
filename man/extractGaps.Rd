\name{extractGaps}
\alias{extractGaps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the Positions of Spacing from Clustered Groups
}
\description{
Very similar to \link[dataVisEasy]{extractClusters} but for the gaps that appear between the clusters when row.groups or col.groups are set in the heatmap functions. Useful in conjunction with \link[dataVisEasy]{ExtractMatrix} if spaces between the groups want to be maintained or widened without actually clustering the data. Output can be supplied to gaps.row and gaps.col in \link[dataVisEasy]{myHeatmap} or gaps.row.spec and gaps.col.spec in \link[dataVisEasy]{myHeatmapByAnnotation}
}
\usage{
extractGaps(data, heatmap, extractRows = FALSE, extractCols = TRUE,
  num_Rows = 5, num_Cols = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows. Should be the same data supplied to the heatmap function
}
  \item{heatmap}{
Output of a heatmap generated from \link[dataVisEasy]{myHeatmap} or \link[dataVisEasy]{myHeatmapByAnnotation}.  The heatmap generated from these functions can be saved to a variable which should serve as the input here. IMPORTANT. The heatmap supplied here MUST have show.colnames and show.rownames = TRUE.
}
  \item{extractRows}{
logical. should gaps along clustered rows be extracted. At least one of extractRows and extractCols must be set to TRUE
}
  \item{extractCols}{
logical. should gaps along clustered columns be extracted. At least one of extractRows and extractCols must be set to TRUE
}
  \item{num_Rows}{
number of clusters along the rows for which the gaps should be extract. Default set to 5.
}
  \item{num_Cols}{
number of clusters along the columns for which the gaps should be extract. Default set to 5.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a numericvector or a list object of two numeric vectors (in a case where both rows and columns are extracted) where the positions of each gap along the rows/columns is indicated.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
~~Alison Moss~~
}
\note{
a similar output can be obtained from using \link[dataVisEasy]{extractClusters} where gaps would equal the cumulative sum of the frequency of each cluster from the output (cumsum(table(output)))
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \link[dataVisEasy]{extractClusters}, \link[dataVisEasy]{ExtractMatrix}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, heatmap, extractRows = FALSE, extractCols = TRUE,
    num_Rows = 5, num_Cols = 5)
{
    if (extractRows == TRUE) {
        gene.groups <- cutree(heatmap$tree_row, k = num_Rows)
        gene.tree.order <- heatmap$tree_row[["labels"]][heatmap$tree_row[["order"]]]
        rowvec <- c(0)
        ind = 1
        for (i in 1:length(unique(gene.groups))) {
            which.group <- gene.groups[which(names(gene.groups) ==
                gene.tree.order[ind])]
            names <- names(which(gene.groups == which.group))
            rowvec <- c(rowvec, rowvec[i] + length(names))
            ind <- ind + length(names)
        }
        rowvec <- rowvec[-length(rowvec)]
        rowvec <- rowvec[-1]
    }
    if (extractCols == TRUE) {
        sample.groups <- cutree(heatmap$tree_col, k = num_Cols)
        sample.tree.order <- heatmap$tree_col[["labels"]][heatmap$tree_col[["order"]]]
        colvec <- c(0)
        ind = 1
        for (i in 1:length(unique(sample.groups))) {
            which.group <- sample.groups[which(names(sample.groups) ==
                sample.tree.order[ind])]
            names <- names(which(sample.groups == which.group))
            colvec <- c(colvec, colvec[i] + length(names))
            ind <- ind + length(names)
        }
        colvec <- colvec[-length(colvec)]
        colvec <- colvec[-1]
    }
    if (extractRows == TRUE & extractCols == TRUE) {
        return(c(list(Row_Groups = rowvec, Col_Groups = colvec)))
    }
    if (extractRows == TRUE & extractCols == FALSE) {
        return(rowvec)
    }
    if (extractRows == FALSE & extractCols == TRUE) {
        return(colvec)
    }
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
