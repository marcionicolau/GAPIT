\name{GAPIT.Block}
\alias{GAPIT.Block}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Block(Z, GA, KG)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{GA}{
%%     ~~Describe \code{GA} here~~
}
  \item{KG}{
%%     ~~Describe \code{KG} here~~
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
function (Z, GA, KG) 
{
    taxa.Z = as.matrix(Z[1, -1])
    taxa.GA = as.matrix(GA[, 1])
    position = taxa.GA \%in\% taxa.Z
    GAU = cbind(GA, 2)
    GAU[position, 3] = 1
    grp.12 = as.matrix(as.vector(unique(GAU[, 2])))
    grp.1 = as.matrix(as.vector(unique(GAU[which(GAU[, 3] == 
        1), 2])))
    grp.2 = as.matrix(as.vector(setdiff(grp.12, grp.1)))
    numWithout = length(grp.2)
    order.1 = 1:length(grp.1)
    order.2 = 1:length(grp.2)
    if (numWithout > 0) 
        grpblock = as.matrix(rbind(cbind(grp.1, 1, order.1), 
            cbind(grp.2, 2, order.2)))
    if (numWithout == 0) 
        grpblock = as.matrix(cbind(grp.1, 1, order.1), )
    order.block = order(as.matrix(GAU[, 3]))
    colnames(grpblock) = c("grp", "block", "ID")
    GAU1 <- merge(GAU[order.block, ], grpblock, by.x = "X2", 
        by.y = "grp")
    GAU1[, 4] = (as.numeric(GAU1[, 3]) + as.numeric(GAU1[, 4]))/2
    GAU = GAU1[, c(2, 1, 4, 5)]
    KW = KG[grp.1, grp.1]
    KO = KG[grp.2, grp.2]
    KWO = KG[grp.1, grp.2]
    return(list(GAU = GAU, KW = KW, KO = KO, KWO = KWO))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
