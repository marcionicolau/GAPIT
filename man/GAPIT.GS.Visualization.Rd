\name{GAPIT.GS.Visualization}
\alias{GAPIT.GS.Visualization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.GS.Visualization(gsBLUP = gsBLUP, BINS = BINS, name.of.trait = name.of.trait)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gsBLUP}{
%%     ~~Describe \code{gsBLUP} here~~
}
  \item{BINS}{
%%     ~~Describe \code{BINS} here~~
}
  \item{name.of.trait}{
%%     ~~Describe \code{name.of.trait} here~~
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
function (gsBLUP = gsBLUP, BINS = BINS, name.of.trait = name.of.trait) 
{
    nBin = BINS
    BLUP = gsBLUP[, 5]
    PEV = gsBLUP[, 6]
    if (BLUP[1] == "NaN") {
        warning("It was not converged. BLUP was not created!")
    }
    if (BLUP[1] != "NaN") {
        BLUP.max = try(max(BLUP))
        BLUP.min = try(min(BLUP))
        if (inherits(BLUP.max, "try-error")) 
            return()
        range.BLUP = BLUP.max - BLUP.min
        range.PEV = max(PEV) - min(PEV)
        interval.BLUP = range.BLUP/nBin
        interval.PEV = range.PEV/nBin
        bin.BLUP = floor(BLUP/max(BLUP) * nBin) * max(BLUP)/nBin
        bin.PEV = floor(PEV/max(PEV) * nBin) * max(PEV)/nBin
        distinct.BLUP = unique(bin.BLUP)
        distinct.PEV = unique(bin.PEV)
        if ((length(distinct.BLUP) < 2) | (length(distinct.PEV) < 
            2)) 
            return()
        Position.BLUP = match(bin.BLUP, distinct.BLUP, nomatch = 0)
        Position.PEV = match(bin.PEV, distinct.PEV, nomatch = 0)
        value = matrix(1, length(Position.BLUP))
        KG <- (tapply(as.numeric(value), list(Position.BLUP, 
            Position.PEV), sum))
        rownames(KG) = round(distinct.BLUP, digits = 4)
        colnames(KG) = round(distinct.PEV, digits = 4)
        rownames(KG) <- rownames(KG)[order(as.numeric(rownames(KG)))]
        colnames(KG) <- colnames(KG)[order(as.numeric(colnames(KG)))]
        pdf(paste("GAPIT.", name.of.trait, ".GPS.BLUPvsPEV", 
            ".pdf", sep = ""), width = 9)
        par(mar = c(5, 5, 6, 5))
        nba_heatmap <- heatmap(KG, Rowv = NA, Colv = NA, col = rev(heat.colors(256)), 
            scale = "column", xlab = "PEV", ylab = "BLUP", main = " ")
        dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
