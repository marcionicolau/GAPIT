\name{GAPIT.get.LL}
\alias{GAPIT.get.LL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.get.LL(pheno, geno = NULL, snp.pool, X0 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pheno}{
%%     ~~Describe \code{pheno} here~~
}
  \item{geno}{
%%     ~~Describe \code{geno} here~~
}
  \item{snp.pool}{
%%     ~~Describe \code{snp.pool} here~~
}
  \item{X0}{
%%     ~~Describe \code{X0} here~~
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
function (pheno, geno = NULL, snp.pool, X0 = NULL) 
{
    y = pheno
    p = 0
    if (is.null(X0)) {
        X0 = matrix(1, nrow(snp.pool), 1)
    }
    X = X0
    K.X.svd = svd(snp.pool)
    d = K.X.svd$d
    d = d[d > 1e-08]
    d = d^2
    U1 = K.X.svd$u
    U1 = U1[, 1:length(d)]
    n = nrow(U1)
    I = diag(1, nrow(U1))
    U1TX = crossprod(U1, X)
    U1TY = crossprod(U1, y)
    yU1TY <- y - U1 \%*\% U1TY
    XU1TX <- X - U1 \%*\% U1TX
    IUU = (I - tcrossprod(U1, U1))
    IUX = crossprod(IUU, X)
    IUY = crossprod(IUU, y)
    for (m in seq(-5, 5, by = 0.1)) {
        p = p + 1
        delta <- exp(m)
        beta1 = 0
        for (i in 1:length(d)) {
            one = matrix(U1TX[i, ], nrow = 1)
            beta = crossprod(one, (one/(d[i] + delta)))
            beta1 = beta1 + beta
        }
        beta2 = 0
        for (i in 1:nrow(U1)) {
            one = matrix(IUX[i, ], nrow = 1)
            dim(one)
            beta = crossprod(one, one)
            beta2 = beta2 + beta
        }
        beta2 <- beta2/delta
        beta3 = 0
        for (i in 1:length(d)) {
            one1 = matrix(U1TX[i, ], nrow = 1)
            one2 = matrix(U1TY[i, ], nrow = 1)
            beta = crossprod(one1, (one2/(d[i] + delta)))
            beta3 = beta3 + beta
        }
        beta4 = 0
        for (i in 1:nrow(U1)) {
            one1 = matrix(IUX[i, ], nrow = 1)
            one2 = matrix(IUY[i, ], nrow = 1)
            beta = crossprod(one1, one2)
            beta4 = beta4 + beta
        }
        beta4 <- beta4/delta
        zw1 <- try(solve(beta1 + beta2))
        if (inherits(zw1, "try-error")) {
            zw1 <- ginv(beta1 + beta2)
        }
        zw2 = (beta3 + beta4)
        beta = crossprod(zw1, zw2)
        part11 <- n * log(2 * 3.14)
        part12 <- 0
        for (i in 1:length(d)) {
            part12_pre = log(d[i] + delta)
            part12 = part12 + part12_pre
        }
        part13 <- (nrow(U1) - length(d)) * log(delta)
        part1 <- -1/2 * (part11 + part12 + part13)
        part21 <- nrow(U1)
        part221 = 0
        for (i in 1:length(d)) {
            one1 = matrix(U1TX[i, ], nrow = 1)
            one2 = matrix(U1TY[i, ], nrow = 1)
            part221_pre = (one2 - one1 \%*\% beta)^2/(d[i] + delta)
            part221 = part221 + part221_pre
        }
        part222 = 0
        for (i in 1:n) {
            one1 = matrix(XU1TX[i, ], nrow = 1)
            one2 = matrix(yU1TY[i, ], nrow = 1)
            part222_pre = ((one2 - one1 \%*\% beta)^2)/delta
            part222 = part222 + part222_pre
        }
        part22 <- n * log((1/n) * (part221 + part222))
        part2 <- -1/2 * (part21 + part22)
        LL <- part1 + part2
        part1 <- 0
        part2 <- 0
        if (p == 1) {
            beta.save = beta
            delta.save = delta
            LL.save = LL
        }
        else {
            if (LL > LL.save) {
                beta.save = beta
                delta.save = delta
                LL.save = LL
            }
        }
    }
    beta = beta.save
    delta = delta.save
    LL = LL.save
    names(delta) = NULL
    names(LL) = NULL
    U1TX = crossprod(U1, X)
    U1TY = crossprod(U1, y)
    sigma_a1 = 0
    for (i in 1:length(d)) {
        one1 = matrix(U1TX[i, ], nrow = 1)
        one2 = matrix(U1TY[i, ], nrow = 1)
        sigma_a1_pre = (one2 - one1 \%*\% beta)^2/(d[i] + delta)
        sigma_a1 = sigma_a1 + sigma_a1_pre
    }
    IU = I - tcrossprod(U1, U1)
    IUX = crossprod(IU, X)
    IUY = crossprod(IU, y)
    sigma_a2 = 0
    for (i in 1:nrow(U1)) {
        one1 = matrix(IUX[i, ], nrow = 1)
        one2 = matrix(IUY[i, ], nrow = 1)
        sigma_a2_pre <- (one2 - one1 \%*\% beta)^2
        sigma_a2 = sigma_a2 + sigma_a2_pre
    }
    sigma_a2 <- sigma_a2/delta
    sigma_a <- 1/n * (sigma_a1 + sigma_a2)
    sigma_e <- delta * sigma_a
    return(list(beta = beta, delta = delta, LL = LL, vg = sigma_a, 
        ve = sigma_e))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
