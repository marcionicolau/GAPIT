\name{GAPIT.QC}
\alias{GAPIT.QC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.QC(Y = NULL, KI = NULL, GT = NULL, CV = NULL, Z = NULL, GK = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{KI}{
%%     ~~Describe \code{KI} here~~
}
  \item{GT}{
%%     ~~Describe \code{GT} here~~
}
  \item{CV}{
%%     ~~Describe \code{CV} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{GK}{
%%     ~~Describe \code{GK} here~~
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
function (Y = NULL, KI = NULL, GT = NULL, CV = NULL, Z = NULL, 
    GK = NULL) 
{
    print("Removing duplicates...")
    Y <- GAPIT.RemoveDuplicate(Y)
    CV <- GAPIT.RemoveDuplicate(CV)
    GK <- GAPIT.RemoveDuplicate(GK)
    if (!is.null(Z)) 
        Z <- GAPIT.RemoveDuplicate(Z)
    print("Removing NaN...")
    Y <- Y[which(Y[, 2] != "NaN"), ]
    print("Remove duplicates for GT...")
    if (!is.null(GT)) 
        taxa.kept <- unique(GT[, 1])
    print("Remove duplicates for KI...")
    if (!is.null(KI)) {
        taxa.all <- KI[, 1]
        taxa.uniqe <- unique(taxa.all)
        position <- match(taxa.uniqe, taxa.all, nomatch = 0)
        position.addition <- cbind(1, t(1 + position))
        KI <- KI[position, position.addition]
    }
    if (!is.null(KI)) {
        taxa.all <- KI[, 1]
        position <- order(taxa.all)
        position.addition <- cbind(1, t(1 + position))
        KI <- KI[position, position.addition]
    }
    print("Remove duplicates for Z (column wise)...")
    if (!is.null(Z)) {
        taxa.all <- as.matrix(Z[1, ])
        taxa.uniqe <- intersect(taxa.all, taxa.all)
        position <- match(taxa.uniqe, taxa.all, nomatch = 0)
        Z <- Z[, position]
    }
    print("Maching Z with Kinship colwise...")
    if (!is.null(KI)) {
        taxa.all <- KI[, 1]
        taxa.kinship <- unique(taxa.all)
    }
    if (!is.null(Z) & !is.null(KI)) {
        taxa.Z <- as.matrix(Z[1, ])
        taxa.Z_K_common <- ifelse(is.null(KI), taxa.Z, intersect(taxa.kinship, 
            taxa.Z))
        Z <- cbind(Z[, 1], Z[, match(taxa.Z_K_common, taxa.Z, 
            nomatch = 0)])
        print("Maching Z without origin...")
        Z1 = Z[-1, -1]
        Z2 = data.frame(Z1)
        Z3 = as.matrix(Z2)
        Z4 = as.numeric(Z3)
        Z5 = matrix(data = Z4, nrow = nrow(Z1), ncol = ncol(Z1))
        RS = rowSums(Z5) > 0
        Z <- Z[c(TRUE, RS), ]
        print("Maching GT and CV...")
        if (length(Z) <= 1) 
            stop("GAPIT says: there is no place to match IDs!")
    }
    taxa = intersect(Y[, 1], Y[, 1])
    if (!is.null(Z)) 
        taxa = intersect(Z[-1, 1], taxa)
    if (!is.null(GT)) 
        taxa = intersect(taxa, taxa.kept)
    if (!is.null(CV)) 
        taxa = intersect(taxa, CV[, 1])
    if (!is.null(GK)) 
        taxa = intersect(taxa, GK[, 1])
    if (length(taxa) <= 1) 
        stop("GAPIT says: There is no individual ID matched to covariate. Please check!")
    if (!is.null(Z)) {
        t = c(TRUE, Z[-1, 1] \%in\% taxa)
        if (length(t) <= 2) 
            stop("GAPIT says: There is no individual ID matched among data. Please check!")
        Z <- Z[t, ]
        print("QC final process...")
        Z1 = Z[-1, -1]
        Z2 = data.frame(Z1)
        Z3 = as.matrix(Z2)
        Z4 = as.numeric(Z3)
        Z5 = matrix(data = Z4, nrow = nrow(Z1), ncol = ncol(Z1))
        CS = colSums(Z5) > 0
        Z <- Z[, c(TRUE, CS)]
    }
    Y <- Y[Y[, 1] \%in\% taxa, ]
    if (!is.null(CV)) 
        CV = CV[CV[, 1] \%in\% taxa, ]
    if (!is.null(GK)) 
        GK = GK[GK[, 1] \%in\% taxa, ]
    if (!is.null(GT)) 
        taxa.kept = data.frame(taxa.kept[taxa.kept \%in\% taxa])
    print("size of taxa.kept")
    print(dim(taxa.kept))
    Y = Y[order(Y[, 1]), ]
    CV = CV[order(CV[, 1]), ]
    if (!is.null(GK)) 
        GK = GK[order(GK[, 1]), ]
    if (!is.null(Z)) 
        Z = Z[c(1, 1 + order(Z[-1, 1])), ]
    position = match(taxa.kept[, 1], GT[, 1], nomatch = 0)
    order.taxa.kept = order(taxa.kept[, 1])
    GTindex = position[order.taxa.kept]
    flag = nrow(Y) == nrow(Z) - 1 & nrow(Y) == nrow(GT) & nrow(Y) == 
        nrow(CV)
    print("GAPIT.QC accomplished successfully!")
    return(list(Y = Y, KI = KI, GT = GT, CV = CV, Z = Z, GK = GK, 
        GTindex = GTindex, flag = flag))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
