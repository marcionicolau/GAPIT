\name{GAPIT.Log}
\alias{GAPIT.Log}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Log(Y = Y, KI = KI, Z = Z, CV = CV, SNP.P3D = SNP.P3D, group.from = group.from, group.to = group.to, group.by = group.by, kinship.cluster = kinship.cluster, kinship.group = kinship.group, ngrid = ngrid, llin = llin, ulim = ulim, esp = esp, name.of.trait = name.of.trait)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{KI}{
%%     ~~Describe \code{KI} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{CV}{
%%     ~~Describe \code{CV} here~~
}
  \item{SNP.P3D}{
%%     ~~Describe \code{SNP.P3D} here~~
}
  \item{group.from}{
%%     ~~Describe \code{group.from} here~~
}
  \item{group.to}{
%%     ~~Describe \code{group.to} here~~
}
  \item{group.by}{
%%     ~~Describe \code{group.by} here~~
}
  \item{kinship.cluster}{
%%     ~~Describe \code{kinship.cluster} here~~
}
  \item{kinship.group}{
%%     ~~Describe \code{kinship.group} here~~
}
  \item{ngrid}{
%%     ~~Describe \code{ngrid} here~~
}
  \item{llin}{
%%     ~~Describe \code{llin} here~~
}
  \item{ulim}{
%%     ~~Describe \code{ulim} here~~
}
  \item{esp}{
%%     ~~Describe \code{esp} here~~
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
function (Y = Y, KI = KI, Z = Z, CV = CV, SNP.P3D = SNP.P3D, 
    group.from = group.from, group.to = group.to, group.by = group.by, 
    kinship.cluster = kinship.cluster, kinship.group = kinship.group, 
    ngrid = ngrid, llin = llin, ulim = ulim, esp = esp, name.of.trait = name.of.trait) 
{
    facto <- list(NULL)
    value <- list(NULL)
    facto[[1]] = "Trait"
    value[[1]] = paste(dim(Y))
    facto[[2]] = "group.by "
    value[[2]] = group.by
    facto[[3]] = "Trait name "
    value[[3]] = name.of.trait
    facto[[4]] = "Kinship"
    value[[4]] = dim(KI)
    facto[[5]] = "Z Matrix"
    value[[5]] = dim(Z)
    facto[[6]] = "Covariate"
    value[[6]] = dim(CV)
    facto[[7]] = "SNP.P3D"
    value[[7]] = SNP.P3D
    facto[[8]] = "Clustering algorithms"
    value[[8]] = kinship.cluster
    facto[[9]] = "Group kinship"
    value[[9]] = kinship.group
    facto[[10]] = "group.from "
    value[[10]] = group.from
    facto[[11]] = "group.to "
    value[[11]] = group.to
    theLog = as.matrix(cbind(facto, value))
    colnames(theLog) = c("Model", "Value")
    file = paste("GAPIT.", name.of.trait, ".Log.csv", sep = "")
    write.table(theLog, file, quote = FALSE, sep = ",", row.names = FALSE, 
        col.names = TRUE)
    return(theLog)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
