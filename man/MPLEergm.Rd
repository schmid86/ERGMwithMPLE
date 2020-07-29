\name{MPLEergm}
\alias{MPLEergm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MPLEergm: MPLE methods for Exponential Random Graph Models
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MPLEergm(formula, R=500, method=c("bootstrap", "godambe"), control.mple = NULL,
control.sim = NULL, control.boot=NULL, checkforseparateddata=TRUE,
nowonit=100, cpu.cores=1, constraints =~.)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Formula for the ERGM. Model construction works like in the ergm package
with the same model terms etc. (for a list of terms, see help("\link[ergm]{ergm-terms}")).
%%     ~~Describe \code{x} here~~
}
\item{R}{The number of bootstrap replications or simulated networks for the variability matrix approximation. More replications lead to more accurate results, but slow down the estimation.}
\item{method}{If "godambe", then the Godambe covariance adjustment is calculated. If "bootstrap", then parametric bootstrap intervals are being calculated.}
\item{control.mple}{ergm controls for \code{ergm} function in the ergm-package. Controls the \code{ergm} function that estimates the MPLE of the observed network. See \link[ergm]{control.ergm} for details.}
\item{control.sim}{simulation control for \code{simulate.ergm} function. See \link[ergm]{control.simulate.ergm} for details.}
\item{control.boot}{ergm controls for \code{ergm} function in the ergm-package. Controls the \code{ergm}-function that estimates the MPLE of the bootstrapped networks. See \link[ergm]{control.ergm} for details.}
\item{checkforseparateddata}{If "TRUE", then the data is checked for separation using the
R-package "safeBinaryRegression" (Konis, 2007). If the data is separated, the MPLE does
not exist and the function will stop and indicate which variable(s) cause the separation.}
\item{nowonit}{Returns: 'Now on Iteration \code{nowonit}' after every \code{nowonit}'s iteration. Keeps track of the progress of the estimation. Default is \code{nowonit=100}.}
\item{cpu.cores}{Number of CPU cores used. Default is \code{cpu.cores=1}. If the number of CPU cores used is greater than 1, the code will automatically parallelize the computation using the \bold{doParallel} R-package.}
\item{constraints}{Constraints of the ERGM. See \link[ergm]{ergm} for details.}



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
\references{C. S. Schmid and B. A. Desmarais. Exponential random graph models with big networks: Maximum pseudolikelihood estimation and the parametric bootstrap. 2017 IEEE International Conference on Big Data (Big Data), pages 116–121, Dec 2017

%% ~put references to the literature/web site here ~
}
\author{Christian S. Schmid (songhyo86@gmail.com) \cr
David R. Hunter (dhunter@stat.psu.edu) \cr
Bruce A. Desmarais (http://brucedesmarais.com/)

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

# estimating model uncertainty using parametric bootstrap
m1 <- MPLEergm(faux.mesa.high ~ edges + nodematch("Sex")+ gwesp(0.25, fixed=TRUE),
R=500, method="bootstrap")
summary(m1)

# estimating model uncertainty using the Godambe adjustment
m2 <- MPLEergm(faux.mesa.high ~ edges + nodematch("Sex")+ gwesp(0.25, fixed=TRUE),
R=500, method="godambe")
summary(m2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
