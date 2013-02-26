\name{GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure}
\alias{GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI = PWI, FDR.Rate = 0.05, FDR.Procedure = "BH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PWI}{
%%     ~~Describe \code{PWI} here~~
}
  \item{FDR.Rate}{
%%     ~~Describe \code{FDR.Rate} here~~
}
  \item{FDR.Procedure}{
%%     ~~Describe \code{FDR.Procedure} here~~
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
function (PWI = PWI, FDR.Rate = 0.05, FDR.Procedure = "BH") 
{
    if (is.null(PWI)) {
        PWIP = NULL
        number.of.significant.SNPs = 0
    }
    if (!is.null(PWI)) {
        if (dim(PWI)[1] == 1) {
            PWIP <- cbind(PWI, PWI[4])
            colnames(PWIP)[9] <- "FDR_Adjusted_P-values"
        }
        if (dim(PWI)[1] > 1) {
            res <- mt.rawp2adjp(PWI[, 4], FDR.Procedure)
            adjp <- res$adjp[order(res$index), ]
            PWIP <- cbind(PWI, adjp[, 2])
            PWIP <- PWIP[order(PWIP[, 4]), ]
            colnames(PWIP)[9] <- "FDR_Adjusted_P-values"
        }
    }
    return(list(PWIP = PWIP))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
