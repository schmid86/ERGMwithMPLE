\name{summary.MPLEergm}
\alias{summary.MPLEergm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing ERGM Model Fits
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.MPLEergm(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{A MPLEergm model object
%%     ~~Describe \code{x} here~~
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
\author{Christian Schmid (songhyo86@gmail.com)
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
# load Faux Mesa High friendship data
data(faux.mesa.high)

# estimate model
m1 <- MPLEergm(faux.mesa.high ~ edges + nodematch("Sex")+ gwesp(0.25, fixed=TRUE),
R=500, method="bootstrap")
summary(m1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
