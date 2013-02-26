\name{GAPIT.Memory.Object}
\alias{GAPIT.Memory.Object}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Memory.Object(name.of.trait = "Trait")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (name.of.trait = "Trait") 
{
    print(paste("R is using", memory.size(), "MB out of limit", 
        memory.limit(), "MB"))
    object.sizes <- function() {
        return(rev(sort(sapply(ls(envir = .GlobalEnv), function(object.name) object.size(get(object.name))))))
    }
    memory = object.sizes()
    file = paste("GAPIT.", name.of.trait, ".Memory.Object.csv", 
        sep = "")
    write.table(memory, file, quote = FALSE, sep = ",", row.names = TRUE, 
        col.names = TRUE)
    pdf(paste("GAPIT.", name.of.trait, ".Memory.Object.pdf", 
        sep = ""))
    barplot(object.sizes(), main = "Memory usage by object", 
        ylab = "Bytes", xlab = "Variable name", col = heat.colors(length(object.sizes())))
    dotchart(object.sizes(), main = "Memory usage by object", 
        xlab = "Bytes")
    pie(object.sizes(), main = "Memory usage by object")
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
