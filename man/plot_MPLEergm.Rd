\name{plot.MPLEergm}
\alias{plot.MPLEergm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot density and trace plots for MPLEergm objects
%%  ~~function to do ... ~~
}
\description{Trace plots and density plots of model parameters for each simulated network. The values have been subtracted by the observed values. The solid red line indicates the variable's value in the observed network. Ideally, this line lies in the center of the plots. The solid line also indicates the expected value obtained from a distribution described by the MLE and/or the Method of Moments estimator.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.MPLEergm(object, vars.per.page=3, cex.axis=1, cex.lab=1, cex.main=1, lwd=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{A MPLEergm object to be diagnosed
%%     ~~Describe \code{x} here~~
}
\item{vars.per.page}{Number of rows (one variable per row) per plotting page. Default is \code{vars.per.page=3}.}
\item{cex.axis}{Specify the size of the tick label numbers/text with a numeric value of length 1.}
\item{cex.lab}{Specify the size of the axis label text with a numeric value of length 1.}
\item{cex.main}{Specify the size of the title text with a numeric value of length 1.}
\item{lwd}{line width relative to the default (default=1)}
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
plot(m1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
