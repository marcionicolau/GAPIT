\name{GAPIT.Compression.Visualization}
\alias{GAPIT.Compression.Visualization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Compression.Visualization(Compression = Compression, name.of.trait = name.of.trait)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Compression}{
%%     ~~Describe \code{Compression} here~~
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
function (Compression = Compression, name.of.trait = name.of.trait) 
{
    if (length(Compression) <= 6) 
        Compression = t(as.matrix(Compression[which(Compression[, 
            4] != "NULL" | Compression[, 4] != "NaN"), ]))
    if (length(Compression) == 6) 
        Compression = matrix(Compression, 1, 6)
    if (length(Compression) > 6) 
        Compression = Compression[which(Compression[, 4] != "NULL" | 
            Compression[, 4] != "NaN"), ]
    if (length(Compression) < 1) 
        return()
    print("Pie chart")
    LL = as.numeric(Compression[, 4])
    Compression.best = Compression[1, ]
    variance = as.numeric(Compression.best[5:6])
    colors <- c("grey50", "grey70")
    labels0 <- round(variance/sum(variance) * 100, 1)
    labels <- paste(labels0, "\%", sep = "")
    LL.best0 = as.numeric(Compression.best[4])
    LL.best = floor(LL.best0 * 100)/100
    theOptimum = paste(c(Compression.best[c(1:3)], LL.best))
    pdf(paste("GAPIT.", name.of.trait, ".Optimum.pdf", sep = ""), 
        width = 14)
    par(mfrow = c(1, 1), mar = c(1, 1, 5, 5), lab = c(5, 5, 7))
    pie(variance, col = colors, labels = labels, angle = 45)
    legend(1, 0.5, c("Genetic variance", "Residual variance"), 
        cex = 1.5, fill = colors)
    text(1.5, 0, "The optimum compression", col = "red")
    for (i in 1:4) {
        text(1.5, -0.1 * i, theOptimum[i], col = "red")
    }
    dev.off()
    Compression = Compression[order(as.numeric(Compression[, 
        3])), ]
    if (length(Compression) == 6) 
        return()
    if (length(unique(Compression[, 3])) > 1) {
        color.vector.basic <- c("red", "blue", "black", "blueviolet", 
            "indianred", "cadetblue", "orange")
        color.vector.addition <- setdiff(c(colors()[grep("red", 
            colors())], colors()[grep("blue", colors())]), color.vector.basic)
        color.vector.addition.mixed <- sample(color.vector.addition, 
            max(0, ((length(unique(Compression[, 1])) * length(unique(Compression[, 
                2]))) - length(color.vector.basic))))
        color.vector <- c(color.vector.basic, color.vector.addition.mixed)
        line.vector <- rep(1:(length(unique(Compression[, 1])) * 
            length(unique(Compression[, 2]))))
        pdf(paste("GAPIT.", name.of.trait, ".Compression.multiple.group.", 
            ".pdf", sep = ""), width = 14)
        par(mfrow = c(2, 3), mar = c(5, 5, 1, 1), lab = c(5, 
            5, 7))
        k <- 1
        for (i in 1:length(unique(Compression[, 1]))) {
            for (j in 1:length(unique(Compression[, 2]))) {
                if ((i == 1) & (j == 1)) {
                  Compression.subset <- Compression[which((Compression[, 
                    1] == as.character(unique(Compression[, 1])[i])) & 
                    (Compression[, 2] == as.character(unique(Compression[, 
                      2])[j]))), ]
                  x <- as.numeric(Compression.subset[, 3])
                  y <- as.numeric(Compression.subset[, 4])
                  plot(y ~ x, type = "l", pch = 30, lty = line.vector[i], 
                    ylim = c(min(as.numeric(Compression[, 4])), 
                      max(as.numeric(Compression[, 4]))), xlim = c(min(as.numeric(Compression[, 
                      3])), max(as.numeric(Compression[, 3]))), 
                    col = color.vector[j], xlab = "Number of Groups", 
                    ylab = "-2Log Likelihood", )
                  label = paste(c(as.character(unique(Compression[, 
                    1]))[k], " ", as.character(unique(Compression[, 
                    2]))[j]), collapse = "")
                }
                if ((i != 1) | (j != 1)) {
                  k <- k + 1
                  Compression.subset <- Compression[which((Compression[, 
                    1] == as.character(unique(Compression[, 1])[i])) & 
                    (Compression[, 2] == as.character(unique(Compression[, 
                      2])[j]))), ]
                  x <- as.numeric(Compression.subset[, 3])
                  y <- as.numeric(Compression.subset[, 4])
                  lines(y ~ x, type = "l", pch = 30, lty = line.vector[i], 
                    col = color.vector[j])
                  label = c(label, paste(c(as.character(unique(Compression[, 
                    1]))[i], " ", as.character(unique(Compression[, 
                    2]))[j]), collapse = ""))
                }
            }
        }
        k <- 1
        for (i in 1:length(unique(Compression[, 1]))) {
            for (j in 1:length(unique(Compression[, 2]))) {
                if ((i == 1) & (j == 1)) {
                  Compression.subset <- Compression[which((Compression[, 
                    1] == as.character(unique(Compression[, 1])[i])) & 
                    (Compression[, 2] == as.character(unique(Compression[, 
                      2])[j]))), ]
                  x <- as.numeric(Compression.subset[, 3])
                  y <- as.numeric(Compression.subset[, 5])
                  plot(y ~ x, type = "l", pch = 17, lty = line.vector[i], 
                    ylim = c(min(as.numeric(Compression[, 5])), 
                      max(as.numeric(Compression[, 5]))), xlim = c(min(as.numeric(Compression[, 
                      3])), max(as.numeric(Compression[, 3]))), 
                    col = color.vector[j], xlab = "Number of Groups", 
                    ylab = "Genetic Variance", )
                }
                if ((i != 1) | (j != 1)) {
                  k <- k + 1
                  Compression.subset <- Compression[which((Compression[, 
                    1] == as.character(unique(Compression[, 1])[i])) & 
                    (Compression[, 2] == as.character(unique(Compression[, 
                      2])[j]))), ]
                  x <- as.numeric(Compression.subset[, 3])
                  y <- as.numeric(Compression.subset[, 5])
                  lines(y ~ x, type = "l", pch = 17, lty = line.vector[i], 
                    col = color.vector[j])
                }
            }
        }
        k <- 1
        for (i in 1:length(unique(Compression[, 1]))) {
            for (j in 1:length(unique(Compression[, 2]))) {
                if ((i == 1) & (j == 1)) {
                  Compression.subset <- Compression[which((Compression[, 
                    1] == as.character(unique(Compression[, 1])[i])) & 
                    (Compression[, 2] == as.character(unique(Compression[, 
                      2])[j]))), ]
                  x <- as.numeric(Compression.subset[, 3])
                  y <- as.numeric(Compression.subset[, 6])
                  plot(y ~ x, type = "l", pch = 17, ylim = c(min(as.numeric(Compression[, 
                    6])), max(as.numeric(Compression[, 6]))), 
                    xlim = c(min(as.numeric(Compression[, 3])), 
                      max(as.numeric(Compression[, 3]))), col = color.vector[j], 
                    xlab = "Number of Groups", ylab = "Residual Variance", 
                    )
                }
                if ((i != 1) | (j != 1)) {
                  k <- k + 1
                  Compression.subset <- Compression[which((Compression[, 
                    1] == as.character(unique(Compression[, 1])[i])) & 
                    (Compression[, 2] == as.character(unique(Compression[, 
                      2])[j]))), ]
                  x <- as.numeric(Compression.subset[, 3])
                  y <- as.numeric(Compression.subset[, 6])
                  lines(y ~ x, type = "l", pch = 17, lty = line.vector[i], 
                    col = color.vector[j])
                }
            }
        }
        heritablilty.vector <- as.numeric(Compression[, 5])/(as.numeric(Compression[, 
            5]) + as.numeric(Compression[, 6]))
        totalVariance.vector <- as.numeric(as.numeric(Compression[, 
            5]) + as.numeric(Compression[, 6]))
        Compression.h2 <- cbind(Compression, heritablilty.vector, 
            totalVariance.vector)
        k <- 1
        for (i in 1:length(unique(Compression.h2[, 1]))) {
            for (j in 1:length(unique(Compression.h2[, 2]))) {
                if ((i == 1) & (j == 1)) {
                  Compression.subset <- Compression.h2[which((Compression.h2[, 
                    1] == as.character(unique(Compression.h2[, 
                    1])[i])) & (Compression.h2[, 2] == as.character(unique(Compression.h2[, 
                    2])[j]))), ]
                  x <- as.numeric(Compression.subset[, 3])
                  y <- as.numeric(Compression.subset[, 8])
                  plot(y ~ x, type = "l", pch = 17, lty = line.vector[k], 
                    ylim = c(min(as.numeric(Compression.h2[, 
                      8])), max(as.numeric(Compression.h2[, 8]))), 
                    xlim = c(min(as.numeric(Compression.h2[, 
                      3])), max(as.numeric(Compression.h2[, 3]))), 
                    col = color.vector[1], xlab = "Number of Groups", 
                    ylab = "Total Variance", )
                }
                if ((i != 1) | (j != 1)) {
                  k <- k + 1
                  Compression.subset <- Compression.h2[which((Compression.h2[, 
                    1] == as.character(unique(Compression.h2[, 
                    1])[i])) & (Compression.h2[, 2] == as.character(unique(Compression.h2[, 
                    2])[j]))), ]
                  x <- as.numeric(Compression.subset[, 3])
                  y <- as.numeric(Compression.subset[, 8])
                  lines(y ~ x, type = "l", pch = 17, lty = line.vector[i], 
                    col = color.vector[j])
                }
            }
        }
        k <- 1
        for (i in 1:length(unique(Compression[, 1]))) {
            for (j in 1:length(unique(Compression[, 2]))) {
                if ((i == 1) & (j == 1)) {
                  Compression.subset <- Compression.h2[which((Compression.h2[, 
                    1] == as.character(unique(Compression.h2[, 
                    1])[i])) & (Compression.h2[, 2] == as.character(unique(Compression.h2[, 
                    2])[j]))), ]
                  x <- as.numeric(Compression.subset[, 3])
                  y <- as.numeric(Compression.subset[, 7])
                  plot(y ~ x, type = "l", pch = 17, lty = line.vector[k], 
                    ylim = c(min(as.numeric(Compression.h2[, 
                      7])), max(as.numeric(Compression.h2[, 7]))), 
                    xlim = c(min(as.numeric(Compression.h2[, 
                      3])), max(as.numeric(Compression.h2[, 3]))), 
                    col = color.vector[1], xlab = "Number of Groups", 
                    ylab = "Heritability", )
                }
                if ((i != 1) | (j != 1)) {
                  k <- k + 1
                  Compression.subset <- Compression.h2[which((Compression.h2[, 
                    1] == as.character(unique(Compression.h2[, 
                    1])[i])) & (Compression.h2[, 2] == as.character(unique(Compression.h2[, 
                    2])[j]))), ]
                  x <- as.numeric(Compression.subset[, 3])
                  y <- as.numeric(Compression.subset[, 7])
                  lines(y ~ x, type = "l", lty = line.vector[i], 
                    pch = 17, col = color.vector[j])
                }
            }
        }
        legend.col = 1 + floor(length(unique(Compression[, 1])) * 
            length(unique(Compression[, 2]))/20)
        line.style = rep(1:length(unique(Compression[, 1])), 
            each = length(unique(Compression[, 2])))
        line.color = rep(1:length(unique(Compression[, 2])), 
            length(unique(Compression[, 1])))
        plot(0 ~ 0, axes = FALSE, type = "l", ylab = "", xlab = "", 
            frame.plot = FALSE)
        legend("topleft", label, col = color.vector[line.color], 
            lty = line.style, ncol = legend.col, horiz = FALSE)
        dev.off()
    }
    if (length(unique(Compression[, 3])) == 1 & length(unique(Compression[, 
        1])) * length(unique(Compression[, 2])) > 1) {
        pdf(paste("GAPIT.Compression.single.group.", name.of.trait, 
            ".pdf", sep = ""), width = 14)
        par(mfrow = c(2, 2), mar = c(5, 5, 1, 1), lab = c(5, 
            5, 7))
        nkt = length(unique(Compression[, 1]))
        nca = length(unique(Compression[, 2]))
        kvr = rep(c(1:nkt), nca)
        kvc0 = rep(c(1:nca), nkt)
        kvc = as.numeric(t(matrix(kvc0, nca, nkt)))
        kt.name = Compression[1:nkt, 1]
        ca.index = ((1:nca) - 1) * nkt + 1
        ca.name = Compression[ca.index, 2]
        KG <- t(tapply(as.numeric(Compression[, 4]), list(kvr, 
            kvc), mean))
        colnames(KG) = kt.name
        barplot(as.matrix(KG), ylab = "-2 Log Likelihood", beside = TRUE, 
            col = rainbow(length(unique(Compression[, 2]))))
        KG <- t(tapply(as.numeric(Compression[, 5]), list(kvr, 
            kvc), mean))
        colnames(KG) = kt.name
        barplot(as.matrix(KG), ylab = "Genetic varaince", beside = TRUE, 
            col = rainbow(length(unique(Compression[, 2]))))
        KG <- t(tapply(as.numeric(Compression[, 6]), list(kvr, 
            kvc), mean))
        colnames(KG) = kt.name
        barplot(as.matrix(KG), ylab = "Residual varaince", beside = TRUE, 
            col = rainbow(length(unique(Compression[, 2]))))
        KG <- t(tapply(as.numeric(Compression[, 5])/(as.numeric(Compression[, 
            5]) + as.numeric(Compression[, 6])), list(kvr, kvc), 
            mean))
        colnames(KG) = kt.name
        barplot(as.matrix(KG), ylab = "Heritability", beside = TRUE, 
            col = rainbow(length(unique(Compression[, 2]))), 
            ylim = c(0, 1))
        legend("topleft", paste(t(ca.name)), cex = 0.8, bty = "n", 
            fill = rainbow(length(unique(Compression[, 2]))), 
            horiz = TRUE)
        dev.off()
    }
    print("GAPIT.Compression.Visualization accomplished successfully!")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
