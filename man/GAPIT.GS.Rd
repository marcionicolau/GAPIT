\name{GAPIT.GS}
\alias{GAPIT.GS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.GS(KW, KO, KWO, GAU, UW)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KW}{
%%     ~~Describe \code{KW} here~~
}
  \item{KO}{
%%     ~~Describe \code{KO} here~~
}
  \item{KWO}{
%%     ~~Describe \code{KWO} here~~
}
  \item{GAU}{
%%     ~~Describe \code{GAU} here~~
}
  \item{UW}{
%%     ~~Describe \code{UW} here~~
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
function (KW, KO, KWO, GAU, UW) 
{
    UO = try(t(KWO) \%*\% solve(KW) \%*\% UW)
    if (inherits(UO, "try-error")) 
        UO = t(KWO) \%*\% ginv(KW) \%*\% UW
    n = ncol(UW)
    BLUP = data.frame(as.matrix(GAU[, 1:4]))
    BLUP.W = BLUP[which(GAU[, 3] < 2), ]
    order.W = order(as.numeric(as.matrix(BLUP.W[, 4])))
    ID.W = as.numeric(as.matrix(BLUP.W[order.W, 4]))
    n.W = max(ID.W)
    DS.W = diag(n.W)[ID.W, ]
    ind.W = DS.W \%*\% UW
    all.W = cbind(BLUP.W[order.W, ], ind.W)
    all = all.W
    BLUP.O = BLUP[which(GAU[, 3] == 2), ]
    if (nrow(BLUP.O) > 0) {
        order.O = order(as.numeric(as.matrix(BLUP.O[, 4])))
        ID.O = as.numeric(as.matrix(BLUP.O[order.O, 4]))
        n.O = max(ID.O)
        DS.O = diag(n.O)[ID.O, ]
        ind.O = DS.O \%*\% UO
        all.O = cbind(BLUP.O[order.O, ], ind.O)
        all = rbind(all.W, all.O)
    }
    colnames(all) = c("Taxa", "Group", "RefInf", "ID", "BLUP", 
        "PEV")
    return(list(BLUP = all))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
