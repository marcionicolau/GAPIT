\name{GAPIT.QQ}
\alias{GAPIT.QQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.QQ(P.values, plot.type = "log_P_values", name.of.trait = "Trait", DPP = 50000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P.values}{
%%     ~~Describe \code{P.values} here~~
}
  \item{plot.type}{
%%     ~~Describe \code{plot.type} here~~
}
  \item{name.of.trait}{
%%     ~~Describe \code{name.of.trait} here~~
}
  \item{DPP}{
%%     ~~Describe \code{DPP} here~~
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
function (P.values, plot.type = "log_P_values", name.of.trait = "Trait", 
    DPP = 50000) 
{
    print("Sorting p values")
    print(paste("Number of P values: ", length(P.values)))
    if (length(P.values[P.values > 0]) < 1) 
        return(NULL)
    DPP = round(DPP/4)
    P.values <- P.values[order(P.values)]
    print("Setting p_value_quantiles...")
    p_value_quantiles <- (1:length(P.values))/(length(P.values) + 
        1)
    if (plot.type == "log_P_values") {
        log.P.values <- -log10(P.values)
        log.Quantiles <- -log10(p_value_quantiles)
        index = GAPIT.Pruning(log.P.values, DPP = DPP)
        log.P.values = log.P.values[index]
        log.Quantiles = log.Quantiles[index]
        pdf(paste("GAPIT.", name.of.trait, ".QQ-Plot.pdf", sep = ""))
        par(mar = c(5, 5, 5, 5))
        qqplot(log.Quantiles, log.P.values, xlim = c(0, max(log.Quantiles)), 
            ylim = c(0, max(log.P.values)), cex.axis = 1.5, cex.lab = 2, 
            lty = 1, lwd = 1, col = "Blue", xlab = expression(Expected ~ 
                ~-log[10](italic(p))), ylab = expression(Observed ~ 
                ~-log[10](italic(p))), main = paste(name.of.trait, 
                sep = " "))
        abline(a = 0, b = 1, col = "red")
        dev.off()
    }
    if (plot.type == "P_values") {
        pdf(paste("QQ-Plot_", name.of.trait, ".pdf", sep = ""))
        par(mar = c(5, 5, 5, 5))
        qqplot(p_value_quantiles, P.values, xlim = c(0, 1), ylim = c(0, 
            1), type = "l", xlab = "Uniform[0,1] Theoretical Quantiles", 
            lty = 1, lwd = 1, ylab = "Quantiles of P-values from GWAS", 
            col = "Blue", main = paste(name.of.trait, sep = " "))
        abline(a = 0, b = 1, col = "red")
        dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
