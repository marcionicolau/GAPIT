\name{GAPIT.Memory}
\alias{GAPIT.Memory}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Memory(Memory = NULL, Infor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Memory}{
%%     ~~Describe \code{Memory} here~~
}
  \item{Infor}{
%%     ~~Describe \code{Infor} here~~
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
function (Memory = NULL, Infor) 
{
    gc()
    size <- memory.size()
    if (is.null(Memory)) {
        Increased = 0
        Memory = cbind(Infor, size, Increased)
    }
    else {
        Increased = 0
        Memory.current = cbind(Infor, size, Increased)
        Memory = rbind(Memory, Memory.current)
        Memory[nrow(Memory), 3] = as.numeric(as.matrix(Memory[nrow(Memory), 
            2])) - as.numeric(as.matrix(Memory[nrow(Memory) - 
            1, 2]))
    }
    return(Memory)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
