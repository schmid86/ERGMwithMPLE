\name{simulate.MPLEergm}
\alias{simulate.MPLEergm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate new networks from MPLEergm objects
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simulate.MPLEergm(object, nsim=1, seed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{A MPLEergm model object
%%     ~~Describe \code{x} here~~
}
\item{nsim}{Number of networks to be simulated}
\item{seed}{Seed to be set}
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

# simulate 10 new networks from m1
simulate(m1, nsim=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
