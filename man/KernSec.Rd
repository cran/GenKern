\name{KernSec}

\alias{KernSec}

\title{Section from a generalised bi-variate binned kernel density estimate}

\description{Computes section of bi-variate generalised binned kernel density estimate}

\usage{KernSec(x, y, yvalue, gridsize=100, correlation=0, xbandwidth,
	ybandwidth, range.x)}

\arguments{
\item{x}{vector of x values}
\item{y}{vector of y values}
\item{yvalue}{value for y at which f(x,y) is to be evaluated}
\item{gridsize}{integer of bins in the x estimate;  default=100}
\item{correlation}{value corresponding to the x,y correlation; default=cor(x,y)}
\item{xbandwidth}{value for the x window width; default=dpik(x)}
\item{ybandwidth}{value for the y window width; default=dpik(y)}
}

\value{returns two vectors:
\item{xvals}{vector of bin centres}
\item{yden}{vector of density estimates for f(x,y)|y}
}

\author{David Lucy <\email{d.j.lucy@bradford.ac.uk}>
\cr Robert Aykroyd <\email{robert@amsta.leeds.ac.uk}>\url{http://www.amsta.leeds.ac.uk/~robert/}}

\references{Robertson, I. Lucy, D, Baxter, L. Pollard, A.M. Aykroyd, R.G. Carter, A.H.C. Switsur, V.R. and Waterhouse, J.S.(1999) A kernel based Bayesian approach to climatic reconstruction.\emph{Holocene} \bold{9}(4): 495-500}

\note{fast code suitable for use in calibration}

\seealso{\code{\link{KernSur}} \code{\link{per}} \code{\link[KernSmooth]{dpik}}}

\section{Acknowledgements}{Written in collaboration with A.M.Pollard <\email{a.m.pollard@bradford.ac.uk}> with the financial support of the Natural Environment Research Council (NERC) grant GR3/11395}

\examples{
x <- c(2,4,6,8) 			# make up some x-y data
y <- x     
z <- KernSec(x, y, 5, gridsize=100, correlation=0, xbandwidth=0.5,
			     ybandwidth=3, range.x=c(0,10)) 

plot(z$xvals, z$yden, type="l") 	#should plot a wiggly line
}

\keyword{density}
