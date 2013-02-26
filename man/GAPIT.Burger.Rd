\name{GAPIT.Burger}
\alias{GAPIT.Burger}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Burger(Y = NULL, CV = NULL, GK = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{CV}{
%%     ~~Describe \code{CV} here~~
}
  \item{GK}{
%%     ~~Describe \code{GK} here~~
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
function (Y = NULL, CV = NULL, GK = NULL) 
{
    myFaSTREML = GAPIT.get.LL(pheno = matrix(Y[, -1], nrow(Y), 
        1), geno = NULL, snp.pool = as.matrix(GK[, -1]), X0 = as.matrix(cbind(matrix(1, 
        nrow(CV), 1), CV[, -1])))
    REMLs = -2 * myFaSTREML$LL
    delta = myFaSTREML$delta
    vg = myFaSTREML$vg
    ve = myFaSTREML$ve
    return(list(REMLs = REMLs, vg = vg, ve = ve, delta = delta))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
