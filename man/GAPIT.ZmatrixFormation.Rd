\name{GAPIT.ZmatrixFormation}
\alias{GAPIT.ZmatrixFormation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.ZmatrixFormation(Z, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
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
function (Z, Y) 
{
    taxa.Z = as.matrix(Z[-1, 1])
    taxa.Y = as.matrix(Y[, 1])
    taxa.diff = setdiff(taxa.Y, taxa.Z)
    taxa.I = as.matrix(taxa.Y[match(taxa.diff, taxa.Y, nomatch = 0)])
    taxa.Z.col = as.matrix(Z[1, -1])
    Z0 = matrix(data = 0, nrow = nrow(taxa.Z), ncol = nrow(taxa.I))
    Z1 = diag(1, nrow(taxa.I))
    ZC = as.matrix(rbind(Z0, Z1))
    label.row = rbind(as.matrix(Z[, 1]), taxa.I)
    label.col = t(taxa.I)
    position = t(as.matrix(match(taxa.Z.col, taxa.I, nomatch = 0)))
    ZC[1:nrow(taxa.Z), position] = as.matrix(Z[-1, -1])
    dataPart = rbind(label.col, ZC)
    Z = data.frame(cbind(label.row, dataPart))
    return(Z)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
