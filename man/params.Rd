\name{parameters}
\alias{params}
\docType{data}
\title{
List Object Holding Parameters Called by Other Functions
}
\description{
The functions in the dataVisEasy package often call upon variables stored within this list object for annotations, coloring, and to keep consistent throughout the use of the package without having to unnecessarily spec
}
\usage{data("parameters")}
\format{
  The format is:
List of 9
 $ scale.range               : num [1:2] -1 1
 $ scale.colors              : chr [1:3] "blue" "black" "yellow"
 $ n.colors.range            : num 300
 $ annotations               : logi NA
 $ annot_samps               : logi NA
 $ annotations.genes         : logi NA
 $ annot_genes               : logi NA
 $ annot_cols                : logi NA
 $ expression_gradient.colors: chr [1:5] "blue" "lightblue" "gray" "indianred" ...
}
\details{
All of the parameters stored in params can be set and updated using the functions listed below. Likewise, all of these objects can be accessed from the terminal using params$'parameter' where 'parameter' represents any of the list object the user may wish to check on or change.

The functions in the dataVisEasy package often call upon variables stored within this list object for annotations, coloring, and to keep consistent throughout the use of the package without having to unnecessarily specify various parameters in the functions themselves. For example, in the case of heatmaps, instead of having to specify the colors or the range of the heatmap within the function every time it is used, the functions themselves simply refer to variables stored within this list object. Similarly, once annotations for the data in question are set and stored here, all functions in the package will refer directly to that stored object, allowing for ease of use and consistency across all included functions. While on the front end this may take more time to set up the annotations and variables properly, once they are set the user need not worry about keeping track of our specifying these parameters later on. Likewise, these variables can be changed or updated at any time. For example if the user wishes to change the scale range of the heatmap, they need to set it using \link[dataVisEasy]{set_scale.range}. The new values they have set will then be used until it is changed to something else.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\seealso{
See \link[dataVisEasy]{set_scale.range}, \link[dataVisEasy]{set_scale.colors}, \link[dataVisEasy]{set_n.colors.range}, \link[dataVisEasy]{set_annotations}, \link[dataVisEasy]{set_annot_samps}, \link[dataVisEasy]{set_annotations.genes}, \link[dataVisEasy]{set_annot_genes}, \link[dataVisEasy]{set_annot_cols}, and \link[dataVisEasy]{set_expression_gradient.colors} for information on changing and setting these parameters
}

\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{

}
\keyword{datasets}
