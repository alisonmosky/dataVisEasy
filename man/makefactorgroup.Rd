\name{makefactorgroup}
\alias{makefactorgroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a Factor Based On Ordered Annotations
}
\description{
This function is called in \link[dataVisEasy]{myHeatmapByAnnotation} to create a factor with proper factor levels based on one or more annotations. If more than one level is provided, the function will order the provided sample annotations by the first level first, followed by subsequent levels to obtain an ordered grouping.
}
\usage{
makefactorgroup(annots, levels, specify.gaps = NULL, return.gaps = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annots}{
a data frame of annotations
}
  \item{levels}{
a character vector of length 1-3 where each element must be found in the colnames of annots. The first element will be the first annotation to be sorted. Ex. c("Age","Strain","Treatment") will sort first by Age, then by Strain, then by Treatment.
}
  \item{specify.gaps}{
numeric vector the same length as a character vector supplied to levels. The default of NULL will result in a gap in between each annotation specified. Mostly useful in scenarios where 2 or 3 levels are provided and allows flexibility in how they will be displayed. Numeric values provided also work in tandem with each other. For example a vector of c(1,1) for a supplied groupings vector of length 2 will place 1 space in between levels of the first annotation and an additional space in between levels of the second annotaiton. A vector of c(0,1) will still order the samples by the annotations as indicated but will only place gaps in between the levels of the second annotation supplied.
}
  \item{return.gaps}{
logical. Should gaps a vector of gaps be returned. Necessary when called upon in the \link[dataVisEasy]{myHeatmapByAnnotation}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
if return.gaps is set to FALSE (default) will return a factor the same length as the number of rows in the provided annots

if return.gaps is set to TRUE, will return a list object where "factor.group" will contain the factor indicated above and "gaps" will contain a numeric vector of where gaps should be placed in a heatmap according to the specified annotations and spacing.
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
function (annots, levels, specify.gaps = NULL, return.gaps = FALSE)
{
    if (length(levels) == 1) {
        Level.1 <- annots[, which(colnames(annots) == levels[1])]
        factor.group <- factor(Level.1)
        names(factor.group) <- rownames(annots)
        if (is.null(specify.gaps) == FALSE) {
            if (length(specify.gaps) != length(levels)) {
                stop("length of gap specifications not equal to number of levels provided")
            }
            gaps <- rep(cumsum(rev(rev(rle(as.vector(sort(factor.group)))$lengths)[-1])),
                specify.gaps[1])
        }
        else {
            gaps <- cumsum(rev(rev(rle(as.vector(sort(factor.group)))$lengths)[-1]))
        }
    }
    if (length(levels) == 2) {
        Level.1 <- annots[, which(colnames(annots) == levels[1])]
        Level.2 <- annots[, which(colnames(annots) == levels[2])]
        combo <- data.frame((Level.1), (Level.2), paste(Level.1,
            Level.2))
        colnames(combo) <- c("Level1", "Level2", "Combo")
        rownames(combo) <- rownames(annots)
        bylevel1 <- combo[order(combo[, 1]), ]
        bylevel2 <- bylevel1[order(bylevel1[, 2]), ]
        factor.group <- factor(combo[, which(colnames(combo) ==
            "Combo")], levels = c(unique(as.character(bylevel2[,
            which(colnames(bylevel2) == "Combo")]))))
        names(factor.group) <- rownames(annots)
        if (is.null(specify.gaps) == FALSE) {
            if (length(specify.gaps) != length(levels)) {
                stop("length of gap specifications not equal to number of levels provided")
            }
            gaps <- sort(c(rep(cumsum(rev(rev(rle(as.vector(bylevel2$Level1))$lengths)[-1])),
                specify.gaps[1]), rep(cumsum(rev(rev(rle(as.vector(bylevel2$Level2))$lengths)[-1])),
                specify.gaps[2])))
        }
        else {
            gaps <- cumsum(rev(rev(rle(as.vector(bylevel2$Combo))$lengths)[-1]))
        }
    }
    if (length(levels) == 3) {
        Level.1 <- annots[, which(colnames(annots) == levels[1])]
        Level.2 <- annots[, which(colnames(annots) == levels[2])]
        Level.3 <- annots[, which(colnames(annots) == levels[3])]
        combo <- data.frame((Level.1), (Level.2), (Level.3),
            paste(Level.1, Level.2, Level.3))
        colnames(combo) <- c("Level1", "Level2", "Level3", "Combo")
        rownames(combo) <- rownames(annots)
        bylevel1 <- combo[order(combo[, 1]), ]
        bylevel2 <- bylevel1[order(bylevel1[, 2]), ]
        bylevel3 <- bylevel2[order(bylevel2[, 3]), ]
        factor.group <- factor(combo[, which(colnames(combo) ==
            "Combo")], levels = c(unique(as.character(bylevel3[,
            which(colnames(bylevel3) == "Combo")]))))
        names(factor.group) <- rownames(annots)
        if (is.null(specify.gaps) == FALSE) {
            if (length(specify.gaps) != length(levels)) {
                stop("length of gap specifications not equal to number of levels provided")
            }
            gaps <- sort(c(rep(cumsum(rev(rev(rle(as.vector(bylevel3$Level1))$lengths)[-1])),
                specify.gaps[1]), rep(cumsum(rev(rev(rle(as.vector(bylevel3$Level2))$lengths)[-1])),
                specify.gaps[2]), rep(cumsum(rev(rev(rle(as.vector(bylevel3$Level3))$lengths)[-1])),
                specify.gaps[3])))
        }
        else {
            gaps <- cumsum(rev(rev(rle(as.vector(bylevel3$Combo))$lengths)[-1]))
        }
    }
    if (return.gaps == FALSE) {
        return(factor.group)
    }
    if (return.gaps == TRUE) {
        return(list(factor.group = factor.group, gaps = gaps))
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
