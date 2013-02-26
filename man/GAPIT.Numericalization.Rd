\name{GAPIT.Numericalization}
\alias{GAPIT.Numericalization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Numericalization(x, bit = 2, effect = "Add", impute = "None", Create.indicator = FALSE, Major.allele.zero = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{bit}{
%%     ~~Describe \code{bit} here~~
}
  \item{effect}{
%%     ~~Describe \code{effect} here~~
}
  \item{impute}{
%%     ~~Describe \code{impute} here~~
}
  \item{Create.indicator}{
%%     ~~Describe \code{Create.indicator} here~~
}
  \item{Major.allele.zero}{
%%     ~~Describe \code{Major.allele.zero} here~~
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
function (x, bit = 2, effect = "Add", impute = "None", Create.indicator = FALSE, 
    Major.allele.zero = FALSE) 
{
    if (bit == 1) {
        x[x == "X"] = "N"
        x[x == "-"] = "N"
        x[x == "+"] = "N"
        x[x == "/"] = "N"
        x[x == "K"] = "Z"
    }
    if (bit == 2) {
        x[x == "XX"] = "N"
        x[x == "--"] = "N"
        x[x == "++"] = "N"
        x[x == "//"] = "N"
        x[x == "NN"] = "N"
    }
    n = length(x)
    lev = levels(as.factor(x))
    lev = setdiff(lev, "N")
    len = length(lev)
    count = 1:len
    for (i in 1:len) {
        count[i] = length(x[(x == lev[i])])
    }
    if (Major.allele.zero) {
        if (len > 1 & len <= 3) {
            if (bit == 1) {
                count.temp = cbind(count, seq(1:len))
                if (len == 3) 
                  count.temp = count.temp[-3, ]
                count.temp <- count.temp[order(count.temp[, 1], 
                  decreasing = TRUE), ]
                if (len == 3) 
                  order = c(count.temp[, 2], 3)
                else order = count.temp[, 2]
            }
            if (bit == 2) {
                count.temp = cbind(count, seq(1:len))
                if (len == 3) 
                  count.temp = count.temp[-2, ]
                count.temp <- count.temp[order(count.temp[, 1], 
                  decreasing = TRUE), ]
                if (len == 3) 
                  order = c(count.temp[1, 2], 2, count.temp[2, 
                    2])
                else order = count.temp[, 2]
            }
            count = count[order]
            lev = lev[order]
        }
    }
    if (bit == 1 & len == 3) {
        temp = count[2]
        count[2] = count[3]
        count[3] = temp
    }
    position = order(count)
    if (len <= 1 | len > 3) 
        x = 0
    if (len == 2) 
        x = ifelse(x == "N", NA, ifelse(x == lev[1], 0, 2))
    if (bit == 1) {
        if (len == 3) 
            x = ifelse(x == "N", NA, ifelse(x == lev[1], 0, ifelse(x == 
                lev[3], 1, 2)))
    }
    else {
        if (len == 3) 
            x = ifelse(x == "N", NA, ifelse(x == lev[1], 0, ifelse(x == 
                lev[3], 2, 1)))
    }
    if (impute == "Middle") {
        x[is.na(x)] = 1
    }
    if (len == 3) {
        if (impute == "Minor") {
            x[is.na(x)] = position[1] - 1
        }
        if (impute == "Major") {
            x[is.na(x)] = position[len] - 1
        }
    }
    else {
        if (impute == "Minor") {
            x[is.na(x)] = 2 * (position[1] - 1)
        }
        if (impute == "Major") {
            x[is.na(x)] = 2 * (position[len] - 1)
        }
    }
    if (effect == "Dom") 
        x = ifelse(x == 1, 1, 0)
    if (effect == "Left") 
        x[x == 1] = 0
    if (effect == "Right") 
        x[x == 1] = 2
    return(matrix(x, n, 1))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
