\name{parameters}
\alias{params}
\alias{parameters}
\title{
List Object Holding Parameters Called by Other Functions
}
\description{
The functions in the dataVisEasy package often call upon variables stored within this list object for annotations, coloring, and to keep consistent throughout the use of the package without having to unnecessarily spec
}
\format{
  The format is a list of 9 items:
 \describe{
  \item{scale.range}{the scale of heatmaps, default set to c(-1,1), see \link[dataVisEasy]{set_scale.range}}
  \item{scale.colors}{default colors for heatmaps, default set to c("blue","black","yellow"), see \link[dataVisEasy]{set_scale.colors}}
  \item{n.colors.range}{number of distinct colors in a heatmap, default set to 100, see \link[dataVisEasy]{set_n.colors.range}}
  \item{annotations}{Will hold a dataframe of annotations for samples (data columns), see \link[dataVisEasy]{set_annotations}. Initialized as NA.}
  \item{annot_samps}{Will hold a dataframe of annotations to be displayed as tracks along heatmaps, see \link[dataVisEasy]{set_annot_samps}. Initialized as NA.}
  \item{annotations.genes}{Will hold a dataframe of annotations for genes (data rows), see \link[dataVisEasy]{set_annotations.genes}. Initialized as NA.}
  \item{annot_genes}{Will hold a dataframe of gene annotations to be displayed as tracks along heatmaps, see \link[dataVisEasy]{set_annot_genes}. Initialized as NA.}
  \item{annot_cols}{Specific colors to accompany annotations, will hold a list of color assignments, see \link[dataVisEasy]{set_annot_cols}. Initialized as NA. }
  \item{expression_gradiet.colors}{color range for plots colored by expression level, default set to c("blue","lightblue","gray","indianred","firebrick"), see \link[dataVisEasy]{set_expression_gradient.colors}}
  }
}
\details{
All of the parameters stored in params can be set and updated using the functions listed below. Likewise, all of these objects can be accessed from the terminal using params$'parameter' where 'parameter' represents any of the list object the user may wish to check on or change.

The functions in the dataVisEasy package often call upon variables stored within this list object for annotations, coloring, and to keep consistent throughout the use of the package without having to unnecessarily specify various parameters in the functions themselves. For example, in the case of heatmaps, instead of having to specify the colors or the range of the heatmap within the function every time it is used, the functions themselves simply refer to variables stored within this list object. Similarly, once annotations for the data in question are set and stored here, all functions in the package will refer directly to that stored object, allowing for ease of use and consistency across all included functions. While on the front end this may take more time to set up the annotations and variables properly, once they are set the user need not worry about keeping track of our specifying these parameters later on. Likewise, these variables can be changed or updated at any time. For example if the user wishes to change the scale range of the heatmap, they need to set it using \link[dataVisEasy]{set_scale.range}. The new values they have set will then be used until it is changed to something else.
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}

\seealso{
See \link[dataVisEasy]{set_scale.range}, \link[dataVisEasy]{set_scale.colors}, \link[dataVisEasy]{set_n.colors.range}, \link[dataVisEasy]{set_annotations}, \link[dataVisEasy]{set_annot_samps}, \link[dataVisEasy]{set_annotations.genes}, \link[dataVisEasy]{set_annot_genes}, \link[dataVisEasy]{set_annot_cols}, and \link[dataVisEasy]{set_expression_gradient.colors} for information on changing and setting these parameters
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
%\examples{

%}
\keyword{datasets}
