\name{GAPIT.kinship.loiselle}
\alias{GAPIT.kinship.loiselle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.kinship.loiselle(snps, method = "additive", use = "all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snps}{
%%     ~~Describe \code{snps} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{use}{
%%     ~~Describe \code{use} here~~
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
function (snps, method = "additive", use = "all") 
{
    n0 <- sum(snps == 0, na.rm = TRUE)
    nh <- sum(snps == 0.5, na.rm = TRUE)
    n1 <- sum(snps == 1, na.rm = TRUE)
    nNA <- sum(is.na(snps))
    prod(dim(snps))
    if (method == "dominant") {
        flags <- matrix(as.double(rowMeans(snps, na.rm = TRUE) > 
            0.5), nrow(snps), ncol(snps))
        snps[!is.na(snps) && (snps == 0.5)] <- flags[!is.na(snps) && 
            (snps == 0.5)]
    }
    else if (method == "recessive") {
        flags <- matrix(as.double(rowMeans(snps, na.rm = TRUE) < 
            0.5), nrow(snps), ncol(snps))
        snps[!is.na(snps) && (snps == 0.5)] <- flags[!is.na(snps) && 
            (snps == 0.5)]
    }
    else if ((method == "additive") && (nh > 0)) {
        dsnps <- snps
        rsnps <- snps
        flags <- matrix(as.double(rowMeans(snps, na.rm = TRUE) > 
            0.5), nrow(snps), ncol(snps))
        dsnps[!is.na(snps) && (snps == 0.5)] <- flags[is.na(snps) && 
            (snps == 0.5)]
        flags <- matrix(as.double(rowMeans(snps, na.rm = TRUE) < 
            0.5), nrow(snps), ncol(snps))
        rsnps[!is.na(snps) && (snps == 0.5)] <- flags[is.na(snps) && 
            (snps == 0.5)]
        snps <- rbind(dsnps, rsnps)
    }
    if (use == "all") {
        mafs <- matrix(rowMeans(snps, na.rm = TRUE), nrow(snps), 
            ncol(snps))
        snps[is.na(snps)] <- mafs[is.na(snps)]
    }
    else if (use == "complete.obs") {
        mafs <- matrix(rowMeans(snps, na.rm = TRUE), nrow(snps), 
            ncol(snps))
        snps <- snps[rowSums(is.na(snps)) == 0, ]
    }
    mafs_comp <- 1 - mafs
    snps_comp <- 1 - snps
    n <- ncol(snps)
    K <- matrix(nrow = n, ncol = n)
    diag(K) <- 1
    missing <- rep(NA, dim(snps)[1])
    for (i in 1:dim(snps)[1]) {
        missing[i] <- sum(is.na(snps[i, ]))
    }
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            Num_First_Term_1 <- (snps[, i] - mafs[, i]) * (snps[, 
                j] - mafs[, j])
            Num_First_Term_2 <- (snps_comp[, i] - mafs_comp[, 
                i]) * (snps_comp[, j] - mafs_comp[, j])
            First_Term <- sum(Num_First_Term_1) + sum(Num_First_Term_2)
            Num_Second_Term_1 <- mafs[, i] * (1 - mafs[, i])
            Num_Second_Term_2 <- mafs_comp[, i] * (1 - mafs_comp[, 
                i])
            Num_Second_Term_Bias_Correction <- 1/((2 * n) - missing - 
                1)
            Num_Second_Term <- Num_Second_Term_1 + Num_Second_Term_2
            Second_Term <- sum(Num_Second_Term * Num_Second_Term_Bias_Correction)
            Third_Term <- sum(Num_Second_Term)
            f <- (First_Term + Second_Term)/Third_Term
            K[i, j] <- f
            if (K[i, j] < 0) 
                K[i, j] = 0
            K[j, i] <- K[i, j]
        }
    }
    return(K)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
