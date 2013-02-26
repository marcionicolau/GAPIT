\name{GAPIT.HapMap}
\alias{GAPIT.HapMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.HapMap(G, SNP.effect = "Add", SNP.impute = "Middle", heading = TRUE, Create.indicator = FALSE, Major.allele.zero = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
%%     ~~Describe \code{G} here~~
}
  \item{SNP.effect}{
%%     ~~Describe \code{SNP.effect} here~~
}
  \item{SNP.impute}{
%%     ~~Describe \code{SNP.impute} here~~
}
  \item{heading}{
%%     ~~Describe \code{heading} here~~
}
  \item{Create.indicator}{
%%     ~~Describe \code{Create.indicator} here~~
}
  \item{Major.allele.zero}{
%%     ~~Describe \code{Major.allele.zero} here~~
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
function (G, SNP.effect = "Add", SNP.impute = "Middle", heading = TRUE, 
    Create.indicator = FALSE, Major.allele.zero = FALSE) 
{
    print(paste("Converting HapMap format to numerical under model of ", 
        SNP.impute, sep = ""))
    if (heading) {
        GT = t(G[1, -(1:11)])
        GI = G[-1, c(1, 3, 4)]
    }
    else {
        GT = NULL
        GI = G[, c(1, 3, 4)]
    }
    if (heading) 
        colnames(GT) = "taxa"
    colnames(GI) = c("SNP", "Chromosome", "Position")
    GD = NULL
    bit = nchar(as.character(G[2, 12]))
    print("Perform numericalization")
    if (heading) {
        if (!Create.indicator) 
            GD = apply(G[-1, -(1:11)], 1, function(one) GAPIT.Numericalization(one, 
                bit = bit, effect = SNP.effect, impute = SNP.impute, 
                Major.allele.zero = Major.allele.zero))
        if (Create.indicator) 
            GD = t(G[-1, -(1:11)])
    }
    else {
        if (!Create.indicator) 
            GD = apply(G[, -(1:11)], 1, function(one) GAPIT.Numericalization(one, 
                bit = bit, effect = SNP.effect, impute = SNP.impute, 
                Major.allele.zero = Major.allele.zero))
        if (Create.indicator) 
            GD = t(G[, -(1:11)])
    }
    if (is.null(GD)) {
        GT = NULL
        GI = NULL
    }
    if (!Create.indicator) 
        print(paste("Succesfuly finished converting HapMap which has bits of ", 
            bit, sep = ""))
    return(list(GT = GT, GD = GD, GI = GI))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
