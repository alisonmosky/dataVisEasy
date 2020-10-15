\name{reportGenes}
\alias{reportGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Report Expression Information for a given set of genes
}
\description{
Based on supplied expression data set, will return a report detailing expression information for the list of genes supplied to the function. Will report the percent of samples with detectable expression (non-missing data), Percent Samples NA (missing data), and percent of samples high, mid-range, and low expression with respect to the median of each gene to detail the range of expression of the genes in question.
}
\usage{
reportGenes(data, list, ranges = "fixed", fixed.range = 2, weight = 1.25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
numeric data matrix with samples/observations in the columns and genes/variables in the rows
}
  \item{list}{
character vector of genes/variables to be pulled out of the data matrix for report of expression ranges
}
  \item{ranges}{
option dictating how the range the ranges of expression are determined. The default, "fixed" sets the threshold for high and low expression a fixed value away from the median. Will also accept "mad" which sets these thresholds as a function of the median average deviation from the median.
}
  \item{fixed.range}{
if ranges = "fixed". Indicates the threshold away from the median where expression is considered high or low. For the default case where fixed.range = 2, expression > median() + 2 will be considered high expression while expression < median() - 2 will be considered low expression with all values in between considered mid range.
}
  \item{weight}{
if ranges = "mad". Indicates the weight applied to the mad that sets the threshold away from the median where expression is considered high or low. For the default case where weight = 1.25, expression > median() + mad*1.25 will be considered high expression while expression < median() - mad*1.25 will be considered low expression with all values in between considered mid range.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a data frame with six columns: Gene, Percent Samples Detected, Percent Samples NA, Percent High Expressing, Percent Mid-Range Expression, and Percent Low Expressing
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
function (data, list, ranges = "fixed", fixed.range = 2, weight = 1.25)
{
    report <- matrix(nrow = length(list), ncol = 6)
    if (sum(list \%in\% rownames(data)) == 0) {
        stop("matches for list not found in rownames data")
    }
    for (i in 1:length(list)) {
        gene <- data[which(rownames(data) \%in\% list[i]), ]
        percentNA <- sum(is.na(gene))/length(gene)
        percentdetected <- sum(!is.na(gene))/length(gene)
        gene <- gene[which(!is.na(gene))]
        median <- median(gene, na.rm = T)
        if (ranges == "fixed") {
            upper <- median - fixed.range
            lower <- median + fixed.range
        }
        if (ranges == "mad") {
            upper <- median - (mad(gene, na.rm = T) * weight)
            lower <- median + (mad(gene, na.rm = T) * weight)
        }
        high <- round(sum(gene < upper)/length(gene), 2)
        middle <- round(sum(gene > upper & gene < lower)/length(gene),
            2)
        low <- round(sum(gene > lower)/length(gene), 2)
        report[i, 1] <- list[i]
        report[i, 2] <- percentdetected
        report[i, 3] <- percentNA
        report[i, 4] <- high
        report[i, 5] <- middle
        report[i, 6] <- low
    }
    colnames(report) <- c("Gene", "Percent Samples Detected",
        "Percent Samples NA", "Percent High Expressing", "Percent Mid-Range Expression",
        "Percent Low Expressing")
    return(as.data.frame(report))
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
