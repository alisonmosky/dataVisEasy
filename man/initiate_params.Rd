\name{initiate_params}
\alias{initiate_params}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Initiate Parameters
}
\description{
A function that is run upon loading of the package. Brings the params object inthe the global environment so that it can be accessed and updated. Can also be used to reset the params to the default. See \link[dataVisEasy]{params} for the format of these parameters, their defaults, and how to set them.
}
\usage{
initiate_params(parameters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{parameters}{the parameters}
}
\details{
should be run exactly as stated in the Usage section: "initiate_params(parameters)" where parameters is an internal data structure within the dataVisEasy package
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
See \link[dataVisEasy]{params} for the format of these parameters and their defaults.
}
\examples{

state.cols <- brewer.pal(6,"Set1"); names(state.cols) <- LETTERS[1:6]
annot_cols <- list('Cell_Type'=c("FB"="blue","CV"="violet","NONE"="grey"),
                   'Animal'=c("PR1534"="#0571b0","PR1643"="#ca0020","PR1705"="#92c5de","PR1729"="#f4a582"),
                   "State"=c(state.cols))

###set color annotations
set_annot_cols(annot_cols)

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
initiate_params <- function(){
  assign(params,params,envir =  .GlobalEnv)
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
