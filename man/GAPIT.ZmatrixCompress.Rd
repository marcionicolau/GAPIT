\name{GAPIT.ZmatrixCompress}
\alias{GAPIT.ZmatrixCompress}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.ZmatrixCompress(Z, GAU)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{GAU}{
%%     ~~Describe \code{GAU} here~~
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
function (Z, GAU) 
{
    effect.Z = as.matrix(Z[1, -1])
    effect.GAU = as.matrix(GAU[, 1])
    taxa = as.data.frame(Z[-1, 1])
    GAU0 = GAU[effect.GAU \%in\% effect.Z, ]
    order.GAU = order(GAU0[, 1])
    GAU1 <- GAU0[order.GAU, ]
    id.1 = GAU1[which(GAU1[, 3] < 2), 4]
    n = max(as.numeric(as.vector(id.1)))
    x = as.numeric(as.matrix(GAU1[, 4]))
    DS = diag(n)[x, ]
    order.Z = order(effect.Z)
    Z = Z[-1, -1]
    Z <- Z[, order.Z]
    Z <- matrix(as.numeric(as.matrix(Z)), nrow = nrow(Z), ncol = ncol(Z))
    Z = Z \%*\% DS
    Z = data.frame(cbind(taxa, Z))
    Z = Z[order(as.matrix(taxa)), ]
    return(list(Z = Z))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
