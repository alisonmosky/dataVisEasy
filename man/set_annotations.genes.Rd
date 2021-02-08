\name{set_annotations.genes}
\alias{set_annotations.genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set Annotations for Genes
}
\description{
Set gene annotations that are passed to a variety of functions in the dataVisEasy package.
}
\usage{
set_annotations.genes(annotations.genes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annotations.genes}{
a data frame where the genes/rownames of the data matrix are the same as the rownames of the annotations matrix. Colnames of provided annotations data frame will serve as the name of the annotations when specifying a "group" or "groupings.genes" argument for a variety of functions. Helpful if the columns of the data frame are factors but not necessary.
}
}
\details{
Setting up annotations through use of this function allows for a wide range of versatility for all other functions in the package. The rownames of the annotations MUST match with the rownames of input data (an error will be returned if that is not the case) and this allows the user to supply subsets of their data without having to worry about subsetting annotations or data being in the same order as the annotations, all of which is done within the functions themselves. By simply supplying the appropriate arguments with one of the colnames of this annotaitons data frame, data can be grouped and/or colored by these annotations in addition to serving as levels for template matching. Importantly, for functionality purposes NAs are not permitted within the dataframe and if NAs are present they will be coerced to "No_Annot" and will function as their own category.

To specify the colors that will accompany each level of each annotation, see \link[dataVisEasy]{set_annot_cols}. If colors are not specified, default color values will be applied when necessary. (Not all annotations need to have custom colors, those that do not will be given default colors while those that are specified will be shown as such.)

Please see \link[dataVisEasy]{params} for more information about accessing and setting up these variables.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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
See \link[dataVisEasy]{set_annot_genes} and \link[dataVisEasy]{update_annotations.genes} for more usage and update these annotations and see \link[dataVisEasy]{set_annotations} for the corresponding annotations for the samples/columns of the data.

Works the same as for sample annotations, see \link[dataVisEasy]{set_annotations} for examples
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
