\name{myColorRamp5}
\alias{myColorRamp5}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a Color Gradient for a Numeric Vector
}
\description{
Based on a provided numeric vector, returns a color gradient based on quartiles indicating level of expression.
}
\usage{
myColorRamp5(colors, values, percent.mad = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{colors}{
a character vector of length 5 providing the colors for which the gradient will be based on
}
  \item{values}{
numeric vector for which the gradient will be based on
}
  \item{percent.mad}{
percent of the median absolute deviation away from the median value that will determine the binning of the data
}
}
\details{
The function differs from the standard colorRamp in that the input values are binned before being assigned colors to avoid the heavy influence of extreme values on the color gradient. This function separates the values into quartiles based on the median value of the data and the indicated percent.mad which determines how the upper and lower bins are assigned.

The lowest quartile set from colors[1]-colors[2] ranges from
min(values) < x <= median(values) - percent.mad*mad(values)

The lower middle quartile set from colors[2]-colors[3] ranges from
median(values) - percent.mad*mad(values) < x <= median(values)

The higher middle quartile set from colors[3]-colors[4] ranges from
median(values) < x <= median(values) + percent.mad*mad(values)

The highest quartile set from colors[4]-colors[5] ranges from
median(values) + percent.mad*mad(values) < x <= max(values)

(NAs are removed for all calculations)
}
\value{
a vector of color assignments for the values provided
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
function (colors, values, percent.mad = 0.5)
{
    out <- rep(rgb(0, 0, 0), length(values))
    for (i in 1:length(values)) {
        if (is.na(values[i])) {
        }
        else {
            if (values[i] <= (median(values, na.rm = T) - percent.mad *
                mad(values, na.rm = T))) {
                v <- (values[i] - min(values, na.rm = T))/((median(values,
                  na.rm = T) - percent.mad * mad(values, na.rm = T)) -
                  min(values, na.rm = T))
                x <- colorRamp(colors[1:2])(v)
                out[i] <- rgb(x[, 1], x[, 2], x[, 3], maxColorValue = 255)
            }
            if (values[i] > median(values, na.rm = T) - percent.mad *
                mad(values, na.rm = T) & values[i] <= median(values,
                na.rm = T)) {
                v <- (values[i] - (median(values, na.rm = T) -
                  percent.mad * mad(values, na.rm = T)))/(median(values,
                  na.rm = T) - (median(values, na.rm = T) - percent.mad *
                  mad(values, na.rm = T)))
                x <- colorRamp(colors[2:3])(v)
                out[i] <- rgb(x[, 1], x[, 2], x[, 3], maxColorValue = 255)
            }
            if (values[i] <= median(values, na.rm = T) + percent.mad *
                mad(values, na.rm = T) & values[i] > median(values,
                na.rm = T)) {
                v <- (values[i] - median(values, na.rm = T))/((median(values,
                  na.rm = T) + percent.mad * mad(values, na.rm = T)) -
                  median(values, na.rm = T))
                x <- colorRamp(colors[3:4])(v)
                out[i] <- rgb(x[, 1], x[, 2], x[, 3], maxColorValue = 255)
            }
            if (values[i] > median(values, na.rm = T) + percent.mad *
                mad(values, na.rm = T)) {
                v <- (values[i] - (median(values, na.rm = T) +
                  percent.mad * mad(values, na.rm = T)))/(max(values,
                  na.rm = T) - (median(values, na.rm = T) + percent.mad *
                  mad(values, na.rm = T)))
                x <- colorRamp(colors[4:5])(v)
                out[i] <- rgb(x[, 1], x[, 2], x[, 3], maxColorValue = 255)
            }
        }
    }
    return(out)
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
