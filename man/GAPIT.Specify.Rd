\name{GAPIT.Specify}
\alias{GAPIT.Specify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Specify(GI = NULL, GP = NULL, bin.size = 1e+07, inclosure.size = NULL, MaxBP = 1e+10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GI}{
%%     ~~Describe \code{GI} here~~
}
  \item{GP}{
%%     ~~Describe \code{GP} here~~
}
  \item{bin.size}{
%%     ~~Describe \code{bin.size} here~~
}
  \item{inclosure.size}{
%%     ~~Describe \code{inclosure.size} here~~
}
  \item{MaxBP}{
%%     ~~Describe \code{MaxBP} here~~
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
function (GI = NULL, GP = NULL, bin.size = 1e+07, inclosure.size = NULL, 
    MaxBP = 1e+10) 
{
    if (is.null(GP)) 
        return(list(index = NULL, BP = NULL))
    ID.GP = as.numeric(as.vector(GP[, 3])) + as.numeric(as.vector(GP[, 
        2])) * MaxBP
    bin.GP = floor(ID.GP/bin.size)
    binP = as.matrix(cbind(bin.GP, NA, NA, ID.GP, as.numeric(as.vector(GP[, 
        4]))))
    n = nrow(binP)
    binP = binP[order(as.numeric(as.vector(binP[, 5]))), ]
    binP = binP[order(as.numeric(as.vector(binP[, 1]))), ]
    binP[2:n, 2] = binP[1:(n - 1), 1]
    binP[1, 2] = 0
    binP[, 3] = binP[, 1] - binP[, 2]
    ID.GP = binP[binP[, 3] > 0, ]
    ID.GP = ID.GP[order(as.numeric(as.vector(ID.GP[, 5]))), ]
    ID.GP = ID.GP[!is.na(ID.GP[, 4]), 4]
    if (!is.null(inclosure.size)) 
        ID.GP = ID.GP[1:inclosure.size]
    theIndex = NULL
    if (!is.null(GI)) {
        ID.GI = as.numeric(as.vector(GI[, 3])) + as.numeric(as.vector(GI[, 
            2])) * MaxBP
        theIndex = ID.GI \%in\% ID.GP
    }
    return(list(index = theIndex, CB = ID.GP))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
