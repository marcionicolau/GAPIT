\name{GAPIT.Compress}
\alias{GAPIT.Compress}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Compress(KI, kinship.cluster = "average", kinship.group = "Mean", GN = nrow(KI), Timmer, Memory)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KI}{
%%     ~~Describe \code{KI} here~~
}
  \item{kinship.cluster}{
%%     ~~Describe \code{kinship.cluster} here~~
}
  \item{kinship.group}{
%%     ~~Describe \code{kinship.group} here~~
}
  \item{GN}{
%%     ~~Describe \code{GN} here~~
}
  \item{Timmer}{
%%     ~~Describe \code{Timmer} here~~
}
  \item{Memory}{
%%     ~~Describe \code{Memory} here~~
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
function (KI, kinship.cluster = "average", kinship.group = "Mean", 
    GN = nrow(KI), Timmer, Memory) 
{
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "CP start")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "cp start")
    line.names <- KI[, 1]
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Does this change memory0")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "Does this change memory0")
    KI <- KI[, -1]
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "CP distance")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "cp distance")
    cluster.distance.matrix <- hclust(as.dist(2 - KI), method = kinship.cluster)
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "CP cluster")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "cp cluster")
    group.membership <- cutree(cluster.distance.matrix, k = GN)
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "CP cutree")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "cp cutree")
    if (kinship.group == "Mean") {
        x = as.factor(group.membership)
        n = max(as.numeric(as.vector(x)))
        b = diag(n)[x, ]
        KG = t(b) \%*\% as.matrix(KI) \%*\% b
        CT = t(b) \%*\% (0 * as.matrix(KI) + 1) \%*\% b
        KG = as.matrix(KG/CT)
        rownames(KG) = c(1:nrow(KG))
        colnames(KG) = c(1:ncol(KG))
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "CP calculation original")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "cp calculation original")
    }
    else {
        gm = as.factor(group.membership)
        kv = as.numeric(as.matrix(KI))
        kvr = rep(gm, ncol(KI))
        kvc = as.numeric(t(matrix(kvr, nrow(KI), ncol(KI))))
        kInCol = t(rbind(kv, kvr, kvc))
        rm(gm)
        rm(kv)
        rm(kvr)
        rm(kvc)
        rm(KI)
        gc()
        if (kinship.group == "Max") 
            KG <- tapply(kInCol[, 1], list(kInCol[, 2], kInCol[, 
                3]), max)
        if (kinship.group == "Min") 
            KG <- tapply(kInCol[, 1], list(kInCol[, 2], kInCol[, 
                3]), min)
        if (kinship.group == "Median") 
            KG <- tapply(kInCol[, 1], list(kInCol[, 2], kInCol[, 
                3]), median)
    }
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "CP calculation")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "cp calculation")
    GA <- data.frame(cbind(as.character(line.names), as.numeric(group.membership)))
    return(list(GA = GA, KG = KG, Timmer = Timmer, Memory = Memory))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
