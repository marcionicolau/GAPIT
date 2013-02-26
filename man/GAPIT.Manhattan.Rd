\name{GAPIT.Manhattan}
\alias{GAPIT.Manhattan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Manhattan(GI.MP = NULL, name.of.trait = "Trait", plot.type = "Genomewise", DPP = 50000, cutOff = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GI.MP}{
%%     ~~Describe \code{GI.MP} here~~
}
  \item{name.of.trait}{
%%     ~~Describe \code{name.of.trait} here~~
}
  \item{plot.type}{
%%     ~~Describe \code{plot.type} here~~
}
  \item{DPP}{
%%     ~~Describe \code{DPP} here~~
}
  \item{cutOff}{
%%     ~~Describe \code{cutOff} here~~
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
function (GI.MP = NULL, name.of.trait = "Trait", plot.type = "Genomewise", 
    DPP = 50000, cutOff = 0.01) 
{
    print("Manhattan ploting...")
    if (is.null(GI.MP)) 
        return
    GI.MP = matrix(as.numeric(as.matrix(GI.MP)), nrow(GI.MP), 
        ncol(GI.MP))
    GI.MP <- GI.MP[!is.na(GI.MP[, 1]), ]
    GI.MP <- GI.MP[!is.na(GI.MP[, 2]), ]
    GI.MP <- GI.MP[GI.MP[, 3] > 0, ]
    numMarker = nrow(GI.MP)
    bonferroniCutOff = -log10(cutOff/numMarker)
    GI.MP[, 3] <- -log10(GI.MP[, 3])
    y.lim <- ceiling(max(GI.MP[, 3]))
    chm.to.analyze <- unique(GI.MP[, 1])
    chm.to.analyze = chm.to.analyze[order(chm.to.analyze)]
    numCHR = length(chm.to.analyze)
    if (plot.type == "Chromosomewise") {
        print("Manhattan ploting Chromosomewise")
        pdf(paste("GAPIT.", name.of.trait, ".Manhattan-Plot.Chromosomewise.pdf", 
            sep = ""), width = 10)
        par(mar = c(5, 5, 4, 3), lab = c(8, 5, 7))
        for (i in 1:numCHR) {
            subset = GI.MP[GI.MP[, 1] == chm.to.analyze[i], ]
            y.lim <- ceiling(max(subset[, 3]))
            if (length(subset) > 3) {
                x <- as.numeric(subset[, 2])/10^(6)
                y <- as.numeric(subset[, 3])
            }
            else {
                x <- as.numeric(subset[2])/10^(6)
                y <- as.numeric(subset[3])
            }
            order = order(y, decreasing = TRUE)
            y = y[order]
            x = x[order]
            index = GAPIT.Pruning(y, DPP = round(DPP/numCHR))
            x = x[index]
            y = y[index]
            plot(y ~ x, type = "p", ylim = c(0, y.lim), xlim = c(min(x), 
                max(x)), col = "navy", xlab = expression(Base ~ 
                Pairs ~ (x10^-6)), ylab = "-Log Base 10 p-value", 
                main = paste("Chromosome", chm.to.analyze[i], 
                  sep = " "))
            abline(h = bonferroniCutOff, col = "forestgreen")
        }
        dev.off()
        print("manhattan plot on chromosome finished")
    }
    if (plot.type == "Genomewise") {
        print("Manhattan ploting Genomewise")
        GI.MP <- GI.MP[order(GI.MP[, 2]), ]
        GI.MP <- GI.MP[order(GI.MP[, 1]), ]
        color.vector <- rep(c("orangered", "navyblue"), numCHR)
        ticks = NULL
        lastbase = 0
        for (i in chm.to.analyze) {
            index = (GI.MP[, 1] == i)
            ticks <- c(ticks, lastbase + mean(GI.MP[index, 2]))
            GI.MP[index, 2] = GI.MP[index, 2] + lastbase
            lastbase = max(GI.MP[index, 2])
        }
        x0 <- as.numeric(GI.MP[, 2])
        y0 <- as.numeric(GI.MP[, 3])
        z0 <- as.numeric(GI.MP[, 1])
        position = order(y0, decreasing = TRUE)
        index0 = GAPIT.Pruning(y0[position], DPP = DPP)
        index = position[index0]
        x = x0[index]
        y = y0[index]
        z = z0[index]
        pdf(paste("GAPIT.", name.of.trait, ".Manhattan-Plot.Genomewise.pdf", 
            sep = ""), width = 10)
        par(mar = c(5, 5, 5, 1))
        plot(y ~ x, xlab = expression(Chromosome), ylab = expression(-log[10](italic(p))), 
            cex.lab = 2, col = ifelse(z\%\%2 == 0, "orangered", 
                "navy"), axes = FALSE, type = "p", pch = 20, 
            main = paste(name.of.trait, sep = " "))
        abline(h = bonferroniCutOff, col = "forestgreen")
        axis(1, at = ticks, cex.axis = 1.5, labels = chm.to.analyze, 
            tick = F)
        axis(2, at = 1:y.lim, cex.axis = 1.5, labels = 1:y.lim, 
            tick = F)
        box()
        dev.off()
        print("Manhattan done Genomewise")
    }
    print("GAPIT.Manhattan accomplished successfully!")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
