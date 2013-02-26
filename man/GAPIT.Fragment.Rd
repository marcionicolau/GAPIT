\name{GAPIT.Fragment}
\alias{GAPIT.Fragment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Fragment(file.path = NULL, file.from = NULL, file.to = NULL, file.total = NULL, file.G = NULL, file.Ext.G = NULL, seed = 123, SNP.fraction = 1, SNP.effect = "Add", SNP.impute = "Middle", genoFormat = NULL, file.GD = NULL, file.Ext.GD = NULL, file.GM = NULL, file.Ext.GM = NULL, file.fragment = NULL, file = 1, frag = 1, LD.chromosome = NULL, LD.location = NULL, LD.range = NULL, Create.indicator = FALSE, Major.allele.zero = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.path}{
%%     ~~Describe \code{file.path} here~~
}
  \item{file.from}{
%%     ~~Describe \code{file.from} here~~
}
  \item{file.to}{
%%     ~~Describe \code{file.to} here~~
}
  \item{file.total}{
%%     ~~Describe \code{file.total} here~~
}
  \item{file.G}{
%%     ~~Describe \code{file.G} here~~
}
  \item{file.Ext.G}{
%%     ~~Describe \code{file.Ext.G} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{SNP.fraction}{
%%     ~~Describe \code{SNP.fraction} here~~
}
  \item{SNP.effect}{
%%     ~~Describe \code{SNP.effect} here~~
}
  \item{SNP.impute}{
%%     ~~Describe \code{SNP.impute} here~~
}
  \item{genoFormat}{
%%     ~~Describe \code{genoFormat} here~~
}
  \item{file.GD}{
%%     ~~Describe \code{file.GD} here~~
}
  \item{file.Ext.GD}{
%%     ~~Describe \code{file.Ext.GD} here~~
}
  \item{file.GM}{
%%     ~~Describe \code{file.GM} here~~
}
  \item{file.Ext.GM}{
%%     ~~Describe \code{file.Ext.GM} here~~
}
  \item{file.fragment}{
%%     ~~Describe \code{file.fragment} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{frag}{
%%     ~~Describe \code{frag} here~~
}
  \item{LD.chromosome}{
%%     ~~Describe \code{LD.chromosome} here~~
}
  \item{LD.location}{
%%     ~~Describe \code{LD.location} here~~
}
  \item{LD.range}{
%%     ~~Describe \code{LD.range} here~~
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
function (file.path = NULL, file.from = NULL, file.to = NULL, 
    file.total = NULL, file.G = NULL, file.Ext.G = NULL, seed = 123, 
    SNP.fraction = 1, SNP.effect = "Add", SNP.impute = "Middle", 
    genoFormat = NULL, file.GD = NULL, file.Ext.GD = NULL, file.GM = NULL, 
    file.Ext.GM = NULL, file.fragment = NULL, file = 1, frag = 1, 
    LD.chromosome = NULL, LD.location = NULL, LD.range = NULL, 
    Create.indicator = FALSE, Major.allele.zero = FALSE) 
{
    genoFormat = "hapmap"
    if (!is.null(file.GD) & is.null(file.G)) 
        genoFormat = "EMMA"
    if (genoFormat == "hapmap") {
        G = NULL
        if (frag == 1) {
            skip.1 = 0
            G <- try(read.delim(paste(file.path, file.G, file, 
                ".", file.Ext.G, sep = ""), head = FALSE, skip = skip.1, 
                nrows = file.fragment + 1), silent = TRUE)
        }
        else {
            skip.1 <- (frag - 1) * file.fragment + 1
            G <- try(read.delim(paste(file.path, file.G, file, 
                ".", file.Ext.G, sep = ""), head = FALSE, skip = skip.1, 
                nrows = file.fragment), silent = TRUE)
        }
        if (inherits(G, "try-error")) {
            G = NULL
        }
        if (is.null(G)) {
            return(list(GD = NULL, GI = NULL, GT = NULL, linesRead = NULL, 
                GLD = NULL, heading = NULL))
        }
        heading = (frag == 1)
        if (heading) {
            n = nrow(G) - 1
        }
        else {
            n = nrow(G)
        }
        linesRead = n
        if (SNP.fraction < 1) {
            set.seed(seed + (file * 1000) + frag)
            mySample = sample(1:n, max(2, floor(n * SNP.fraction)))
            if (heading) {
                G = G[c(1, (1 + mySample)), ]
            }
            else {
                G = G[mySample, ]
            }
        }
        print("Call hapmap from fragment")
        hm = GAPIT.HapMap(G, SNP.effect = SNP.effect, SNP.impute = SNP.impute, 
            heading = heading, Create.indicator = Create.indicator, 
            Major.allele.zero = Major.allele.zero)
        if (!is.null(LD.chromosome) & !is.null(hm$GD)) {
            index = (G[, 3] == LD.chromosome[1]) & abs((as.numeric(G[, 
                4]) - as.numeric(LD.location[1])) < (as.numeric(LD.range[1])/2))
            GLD = G[index, ]
        }
        else {
            GLD = NULL
        }
        print("hapmap called successfuly from fragment")
        return(list(GD = hm$GD, GI = hm$GI, GT = hm$GT, linesRead = linesRead, 
            GLD = GLD, heading = heading, G = G))
        print("ERROR: It should not get here!!!")
    }
    if (genoFormat == "EMMA") {
        GD = NULL
        skip.1 <- (frag - 1) * file.fragment
        GD.temp <- try(read.table(paste(file.path, file.GD, file, 
            ".", file.Ext.GD, sep = ""), head = TRUE, nrows = 1), 
            silent = TRUE)
        num.SNP <- ncol(GD.temp) - 1
        rm(GD.temp)
        read.in <- min(file.fragment, (num.SNP - skip.1))
        skip.2 <- max((num.SNP - (skip.1 + read.in)), 0)
        GD <- try(read.table(paste(file.path, file.GD, file, 
            ".", file.Ext.GD, sep = ""), head = TRUE, colClasses = c("factor", 
            rep("NULL", skip.1), rep("numeric", read.in), rep("NULL", 
                skip.2))), silent = TRUE)
        GI <- try(read.table(paste(file.path, file.GM, file, 
            ".", file.Ext.GM, sep = ""), head = TRUE, skip = skip.1, 
            nrows = file.fragment), silent = TRUE)
        if (inherits(GD, "try-error")) {
            GD = NULL
            print("File end reached for GD!!!")
        }
        if (inherits(GI, "try-error")) {
            GI = NULL
            print("File end reached for GI!!!")
        }
        if (is.null(GD)) 
            return(list(GD = NULL, GI = NULL, GT = NULL, linesRead = NULL, 
                GLD = NULL))
        GT = GD[, 1]
        GD = GD[, -1]
        linesRead = ncol(GD)
        if (SNP.fraction == 1) 
            return(list(GD = GD, GI = GI, GT = GT, linesRead = linesRead, 
                GLD = NULL))
        if (SNP.fraction < 1) {
            n = ncol(GD)
            set.seed(seed + file)
            sample = sample(1:n, floor(n * SNP.fraction))
            return(list(GD = GD[, sample], GI = GI[sample, ], 
                GT = GT, linesRead = linesRead, GLD = NULL))
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
