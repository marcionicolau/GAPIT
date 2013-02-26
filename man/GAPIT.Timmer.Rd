\name{GAPIT.Timmer}
\alias{GAPIT.Timmer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Timmer(Timmer = NULL, Infor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Timmer}{
%%     ~~Describe \code{Timmer} here~~
}
  \item{Infor}{
%%     ~~Describe \code{Infor} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
%%  ~~who you are~~
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
function (Timmer = NULL, Infor) 
{
    Time <- Sys.time()
    if (is.null(Timmer)) {
        Elapsed = 0
        Timmer = cbind(Infor, Time, Elapsed)
    }
    else {
        Elapsed = 0
        Timmer.current = cbind(Infor, Time, Elapsed)
        Timmer = rbind(Timmer, Timmer.current)
        Timmer[nrow(Timmer), 3] = as.numeric(as.matrix(Timmer[nrow(Timmer), 
            2])) - as.numeric(as.matrix(Timmer[nrow(Timmer) - 
            1, 2]))
    }
    return(Timmer)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
