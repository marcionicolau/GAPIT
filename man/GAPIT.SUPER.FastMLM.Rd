\name{GAPIT.SUPER.FastMLM}
\alias{GAPIT.SUPER.FastMLM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.SUPER.FastMLM(ys, xs, vg, delta, Z = NULL, X0 = NULL, snp.pool = NULL, LD = 0.01, method = "FaST")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ys}{
%%     ~~Describe \code{ys} here~~
}
  \item{xs}{
%%     ~~Describe \code{xs} here~~
}
  \item{vg}{
%%     ~~Describe \code{vg} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{X0}{
%%     ~~Describe \code{X0} here~~
}
  \item{snp.pool}{
%%     ~~Describe \code{snp.pool} here~~
}
  \item{LD}{
%%     ~~Describe \code{LD} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (ys, xs, vg, delta, Z = NULL, X0 = NULL, snp.pool = NULL, 
    LD = 0.01, method = "FaST") 
{
    ys = unlist(ys)
    if (is.null(dim(ys)) || ncol(ys) == 1) 
        ys <- matrix(ys, 1, length(ys))
    if (is.null(dim(xs)) || ncol(xs) == 1) 
        xs <- matrix(xs, 1, length(xs))
    if (is.null(X0)) 
        X0 <- matrix(1, nrow(snp.pool), 1)
    g <- nrow(ys)
    n <- nrow(xs)
    m <- ncol(xs)
    t <- nrow(xs)
    q0 <- ncol(X0)
    q1 <- q0 + 1
    dfs <- matrix(nrow = m, ncol = g)
    stats <- matrix(nrow = m, ncol = g)
    ps <- matrix(nrow = m, ncol = g)
    betavalue <- matrix(nrow = m, ncol = g)
    if (method == "SUPER") {
        LDsqr = sqrt(LD)
        for (j in 1:g) {
            for (i in 1:m) {
                if ((i > 0) & (floor(i/500) == i/500)) 
                  print(paste("SNP: ", i, " ", sep = ""))
                if (min(xs[, i]) == max(xs[, i])) {
                  dfs[i, j] <- n - q1
                  betavalue[i, j] = 0
                  stats[i, j] <- 0
                }
                if (min(xs[, i]) != max(xs[, i])) {
                  snp.corr = cor(xs[, i], snp.pool)
                  index.k = snp.corr <= LDsqr
                  K.X = snp.pool[, index.k]
                  K.X.svd = svd(K.X)
                  d = K.X.svd$d
                  d = d[d > 1e-08]
                  d = d^2
                  U1 = K.X.svd$u
                  U1 = U1[, 1:length(d)]
                  n <- nrow(U1)
                  I = diag(1, nrow(U1))
                  X <- cbind(X0, xs[, i])
                  U <- U1 * matrix(sqrt(1/(d + delta)), nrow(U1), 
                    length(d), byrow = TRUE)
                  Xt <- crossprod(U, X)
                  XX1 <- crossprod(Xt, Xt)
                  XX2 <- crossprod((I - tcrossprod(U1, U1)) \%*\% 
                    X, (I - tcrossprod(U1, U1)) \%*\% X)/delta
                  iXX <- try(solve(XX1 + XX2))
                  if (inherits(iXX, "try-error")) {
                    iXX <- ginv(XX1 + XX2)
                  }
                  U1TX = crossprod(U1, X)
                  beta1 = 0
                  for (ii in 1:length(d)) {
                    one = matrix(U1TX[ii, ], nrow = 1)
                    dim(one)
                    beta = crossprod(one, one)/(d[ii] + delta)
                    beta1 = beta1 + beta
                  }
                  IUX = (I - tcrossprod(U1, U1)) \%*\% X
                  beta2 = 0
                  for (ii in 1:nrow(U1)) {
                    one = matrix(IUX[ii, ], nrow = 1)
                    dim(one)
                    beta = t(one) \%*\% one
                    beta2 = beta2 + beta
                  }
                  beta2 <- beta2/delta
                  U1TY = crossprod(U1, ys[j, ])
                  beta3 = 0
                  for (ii in 1:length(d)) {
                    one1 = matrix(U1TX[ii, ], nrow = 1)
                    one2 = matrix(U1TY[ii, ], nrow = 1)
                    beta = crossprod(one1, one2)/(d[ii] + delta)
                    beta3 = beta3 + beta
                  }
                  IUY = (I - tcrossprod(U1, U1)) \%*\% ys[j, ]
                  beta4 = 0
                  for (ii in 1:nrow(U1)) {
                    one1 = matrix(IUX[ii, ], nrow = 1)
                    one2 = matrix(IUY[ii, ], nrow = 1)
                    beta = crossprod(one1, one2)
                    beta4 = beta4 + beta
                  }
                  beta4 <- beta4/delta
                  beta = ginv(beta1 + beta2) \%*\% (beta3 + beta4)
                  betavalue[i, j] = beta[q1, 1]
                  stats[i, j] <- beta[q1, 1]/sqrt(iXX[q1, q1] * 
                    vg)
                  dfs[i, j] <- n - q1
                }
            }
            ps[, j] <- 2 * pt(abs(stats[, j]), dfs[, j], lower.tail = FALSE)
        }
        return(list(beta = betavalue, ps = ps, stats = stats, 
            dfs = dfs, effect = betavalue))
    }
    if (method == "FaST") {
        K.X.svd = svd(snp.pool)
        d = K.X.svd$d
        d = d[d > 1e-08]
        d = d^2
        U1 = K.X.svd$u
        U1 = U1[, 1:length(d)]
        n <- nrow(U1)
        I = diag(1, nrow(U1))
        U <- U1 * matrix(sqrt(1/(d + delta)), nrow(U1), length(d), 
            byrow = TRUE)
        for (j in 1:g) {
            for (i in 1:m) {
                if ((i > 0) & (floor(i/500) == i/500)) 
                  print(paste("SNP: ", i, " ", sep = ""))
                if (min(xs[, i]) == max(xs[, i])) {
                  dfs[i, j] <- n - q1
                  betavalue[i, j] = 0
                  stats[i, j] <- 0
                }
                if (min(xs[, i]) != max(xs[, i])) {
                  K.X.svd = svd(snp.pool)
                  d = K.X.svd$d
                  d = d[d > 1e-08]
                  d = d^2
                  U1 = K.X.svd$u
                  U1 = U1[, 1:length(d)]
                  n <- nrow(U1)
                  I = diag(1, nrow(U1))
                  X <- cbind(X0, xs[, i])
                  U <- U1 * matrix(sqrt(1/(d + delta)), nrow(U1), 
                    length(d), byrow = TRUE)
                  Xt <- crossprod(U, X)
                  XX1 <- crossprod(Xt, Xt)
                  XX2 <- crossprod((I - tcrossprod(U1, U1)) \%*\% 
                    X, (I - tcrossprod(U1, U1)) \%*\% X)/delta
                  iXX <- try(solve(XX1 + XX2))
                  if (inherits(iXX, "try-error")) {
                    iXX <- ginv(XX1 + XX2)
                  }
                  U1TX = crossprod(U1, X)
                  beta1 = 0
                  for (ii in 1:length(d)) {
                    one = matrix(U1TX[ii, ], nrow = 1)
                    dim(one)
                    beta = crossprod(one, one)/(d[ii] + delta)
                    beta1 = beta1 + beta
                  }
                  IUX = (I - tcrossprod(U1, U1)) \%*\% X
                  beta2 = 0
                  for (ii in 1:nrow(U1)) {
                    one = matrix(IUX[ii, ], nrow = 1)
                    dim(one)
                    beta = crossprod(one, one)
                    beta2 = beta2 + beta
                  }
                  beta2 <- beta2/delta
                  U1TY = crossprod(U1, ys[j, ])
                  beta3 = 0
                  for (ii in 1:length(d)) {
                    one1 = matrix(U1TX[ii, ], nrow = 1)
                    one2 = matrix(U1TY[ii, ], nrow = 1)
                    beta = crossprod(one1, one2)/(d[ii] + delta)
                    beta3 = beta3 + beta
                  }
                  IUY = (I - tcrossprod(U1, U1)) \%*\% ys[j, ]
                  beta4 = 0
                  for (ii in 1:nrow(U1)) {
                    one1 = matrix(IUX[ii, ], nrow = 1)
                    one2 = matrix(IUY[ii, ], nrow = 1)
                    beta = crossprod(one1, one2)
                    beta4 = beta4 + beta
                  }
                  beta4 <- beta4/delta
                  beta = ginv(beta1 + beta2) \%*\% (beta3 + beta4)
                  betavalue[i, j] = beta[q1, 1]
                  stats[i, j] <- beta[q1, 1]/sqrt(iXX[q1, q1] * 
                    vg)
                  dfs[i, j] <- n - q1
                }
            }
            ps[, j] <- 2 * pt(abs(stats[, j]), dfs[, j], lower.tail = FALSE)
        }
        return(list(beta = betavalue, ps = ps, stats = stats, 
            dfs = dfs, effect = betavalue))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
