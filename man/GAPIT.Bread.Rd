\name{GAPIT.Bread}
\alias{GAPIT.Bread}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Bread(Y = NULL, CV = NULL, Z = NULL, KI = NULL, GK = NULL, GD = NULL, GM = NULL, method = NULL, delta = NULL, vg = NULL, ve = NULL, LD = 0.01, GTindex = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{CV}{
%%     ~~Describe \code{CV} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{KI}{
%%     ~~Describe \code{KI} here~~
}
  \item{GK}{
%%     ~~Describe \code{GK} here~~
}
  \item{GD}{
%%     ~~Describe \code{GD} here~~
}
  \item{GM}{
%%     ~~Describe \code{GM} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{vg}{
%%     ~~Describe \code{vg} here~~
}
  \item{ve}{
%%     ~~Describe \code{ve} here~~
}
  \item{LD}{
%%     ~~Describe \code{LD} here~~
}
  \item{GTindex}{
%%     ~~Describe \code{GTindex} here~~
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
function (Y = NULL, CV = NULL, Z = NULL, KI = NULL, GK = NULL, 
    GD = NULL, GM = NULL, method = NULL, delta = NULL, vg = NULL, 
    ve = NULL, LD = 0.01, GTindex = NULL) 
{
    if (method == "GLM") {
        myGAPIT <- GAPIT(Y = Y, CV = CV, Z = Z, KI = KI, GD = GD, 
            GM = GM, group.from = 0, group.to = 0, QC = FALSE, 
            GTindex = GTindex)
        GWAS = myGAPIT$GWAS
        GPS = myGAPIT$GPS
        REMLs = myGAPIT$REMLs
        delta = myGAPIT$ve/myGAPIT$va
        vg = myGAPIT$vg
        ve = myGAPIT$ve
    }
    if (method == "MLM") {
        myGAPIT <- GAPIT(Y = Y, CV = CV, Z = Z, KI = KI, GD = GD, 
            GM = GM, group.from = nrow(Y), group.to = nrow(Y), 
            QC = FALSE, GTindex = GTindex)
        GWAS = myGAPIT$GWAS
        GPS = myGAPIT$GPS
        REMLs = myGAPIT$REMLs
        delta = myGAPIT$ve/myGAPIT$va
        vg = myGAPIT$vg
        ve = myGAPIT$ve
    }
    if (method == "CMLM") {
        myGAPIT <- GAPIT(Y = Y, CV = CV, Z = Z, KI = KI, GD = GD, 
            GM = GM, group.from = 1, group.to = nrow(Y), QC = FALSE, 
            GTindex = GTindex)
        GWAS = myGAPIT$GWAS
        GPS = myGAPIT$GPS
        REMLs = myGAPIT$REMLs
        delta = myGAPIT$ve/myGAPIT$va
        vg = myGAPIT$vg
        ve = myGAPIT$ve
    }
    if (method == "FaST" | method == "SUPER" | method == "DC") {
        GWAS = NULL
        GPS = NULL
        if (!is.null(vg) & !is.null(vg) & is.null(delta)) 
            delta = ve/vg
        if (is.null(vg) & is.null(ve)) {
            myFaSTREML = GAPIT.get.LL(pheno = matrix(Y[, -1], 
                nrow(Y), 1), geno = NULL, snp.pool = as.matrix(GK[, 
                -1]), X0 = as.matrix(cbind(matrix(1, nrow(CV), 
                1), CV[, -1])))
            REMLs = -2 * myFaSTREML$LL
            delta = myFaSTREML$delta
            vg = myFaSTREML$vg
            ve = myFaSTREML$ve
        }
        mySUPERFaST = GAPIT.SUPER.FastMLM(ys = matrix(Y[, -1], 
            nrow(Y), 1), X0 = as.matrix(cbind(matrix(1, nrow(CV), 
            1), CV[, -1])), snp.pool = as.matrix(GK[-1]), xs = as.matrix(GD[GTindex, 
            -1]), vg = vg, delta = delta, LD = LD, method = method)
        GWAS = cbind(GM, mySUPERFaST$ps, mySUPERFaST$stats, mySUPERFaST$dfs, 
            mySUPERFaST$effect)
    }
    return(list(GWAS = GWAS, GPS = GPS, REMLs = REMLs, vg = vg, 
        ve = ve, delta = delta))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
