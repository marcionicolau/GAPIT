\name{GAPIT.kinship.VanRaden}
\alias{GAPIT.kinship.VanRaden}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.kinship.VanRaden(snps, hasInbred = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snps}{
%%     ~~Describe \code{snps} here~~
}
  \item{hasInbred}{
%%     ~~Describe \code{hasInbred} here~~
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
function (snps, hasInbred = TRUE) 
{
    print("Calculating kinship with VanRaden method...")
    nSNP = ncol(snps)
    nInd = nrow(snps)
    n = nInd
    snpMean = apply(snps, 2, mean)
    snps = snps - snpMean
    print("Getting X'X...")
    K = tcrossprod((snps), (snps))
    print("Adjusting...")
    i = 1:n
    j = (i - 1) * n
    index = i + j
    d = K[index]
    DL = min(d)
    DU = max(d)
    floor = min(K)
    K = (K - floor)/(DL - floor)
    MD = (DU - floor)/(DL - floor)
    if (is.na(K[1, 1])) 
        stop("GAPIT says: Missing data is not allowed for numerical genotype data")
    if (MD > 2) 
        K[index] = K[index]/(MD - 1) + 1
    if (MD < 2 & hasInbred) 
        K = 2 * K/((DU - floor)/(DL - floor))
    print("Calculating kinship with VanRaden method: done")
    return(K)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
