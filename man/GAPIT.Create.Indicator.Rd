\name{GAPIT.Create.Indicator}
\alias{GAPIT.Create.Indicator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Create.Indicator(xs, SNP.impute = "Major")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xs}{
%%     ~~Describe \code{xs} here~~
}
  \item{SNP.impute}{
%%     ~~Describe \code{SNP.impute} here~~
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
function (xs, SNP.impute = "Major") 
{
    bit = nchar(as.character(xs[1]))
    if (bit == 1) {
        xss[xss == "xs"] = "N"
        xs[xs == "-"] = "N"
        xs[xs == "+"] = "N"
        xs[xs == "/"] = "N"
        xs[xs == "K"] = "Z"
    }
    if (bit == 2) {
        xs[xs == "xsxs"] = "N"
        xs[xs == "--"] = "N"
        xs[xs == "++"] = "N"
        xs[xs == "//"] = "N"
        xs[xs == "NN"] = "N"
    }
    xs.temp <- xs[-which(xs == "N")]
    frequ <- NULL
    for (i in 1:length(unique(xs.temp))) frequ <- c(frequ, length(which(xs == 
        unique(xs)[i])))
    unique.sorted <- cbind(unique(xs.temp), frequ)
    print("unique.sorted is")
    print(unique.sorted)
    unique.sorted <- unique.sorted[order(unique.sorted[, 2]), 
        ]
    unique.sorted <- unique.sorted[, -2]
    if (SNP.impute == "Major") 
        xs[which(is.na(xs))] = unique.sorted[1]
    if (SNP.impute == "Minor") 
        xs[which(is.na(xs))] = unique.sorted[length(unique.sorted)]
    if (SNP.impute == "Middle") 
        xs[which(is.na(xs))] = unique.sorted[2]
    x.ind <- NULL
    for (i in unique.sorted) {
        x.col <- rep(NA, length(xs))
        x.col[which(xs == i)] <- 1
        x.col[which(xs != i)] <- 0
        x.ind <- cbind(x.ind, x.col)
    }
    return(x.ind)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
