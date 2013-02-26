\name{GAPIT.PCA}
\alias{GAPIT.PCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.PCA(X, taxa, PC.number = min(ncol(X), nrow(X)), file.output = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{taxa}{
%%     ~~Describe \code{taxa} here~~
}
  \item{PC.number}{
%%     ~~Describe \code{PC.number} here~~
}
  \item{file.output}{
%%     ~~Describe \code{file.output} here~~
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
function (X, taxa, PC.number = min(ncol(X), nrow(X)), file.output = TRUE) 
{
    print("Callingg prcomp...")
    PCA.X <- prcomp(X)
    print("Creating PCA graphs...")
    if (file.output & PC.number > 1) {
        pdf("GAPIT.PCA.eigenValue.pdf", width = 12, height = 12)
        par(mar = c(5, 5, 5, 5))
        screeplot(PCA.X, type = "lines")
        dev.off()
        pdf("GAPIT.PCA.pdf", width = 12, height = 12)
        par(mar = c(5, 5, 5, 5))
        maxPlot = min(as.numeric(PC.number[1]), 3)
        for (i in 1:(maxPlot - 1)) {
            for (j in (i + 1):(maxPlot)) {
                plot(PCA.X$x[, i], PCA.X$x[, j], xlab = paste("PC", 
                  i, sep = ""), ylab = paste("PC", j, sep = ""), 
                  pch = 1, col = "red", cex = 1, cex.lab = 1.5, 
                  cex.axis = 1.2, lwd = 2, las = 1)
            }
        }
        dev.off()
    }
    print("Joining taxa...")
    PCs <- cbind(taxa, as.data.frame(PCA.X$x))
    print("Exporting PCs...")
    if (file.output) 
        write.table(PCs, "GAPIT.PCA.csv", quote = FALSE, sep = ",", 
            row.names = FALSE, col.names = TRUE)
    return(list(PCs = PCs, EV = PCA.X$sdev^2, nPCs = NULL))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
