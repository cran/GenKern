\name{KernSec}

\alias{KernSec}

\title{Univariate kernel density estimate}

\description{Computes univariate kernel density estimate using Gaussian kernels which can also use non-equally spaced ordinates and adaptive bandwidths}

\usage{KernSec(x, gridsize=100, xbandwidth, range.x, na.rm=FALSE)}

\arguments{
\item{x}{vector of \code{x} values}
\item{xgridsize}{integer for number of ordinates at which to calculate the smoothed estimate: default=100}
\item{xbandwidth}{value of \code{x} window width, or vector of local window widths: default=\code{dpik(x)}}
\item{range.x}{total range of the estimate in the \code{x} dimension, or a vector giving the \code{x} ordinates: default=\code{range +- 1.5 * mean bandwidth}}
\item{na.rm}{NA behaviour: \code{TRUE} drops cases with NA's, \code{FALSE} stops function with a warning if NA's are detected: default=\code{FALSE}}
}

\value{returns two vectors:
\item{xvals}{vector of ordinates}
\item{yden}{vector of density estimates corresponding to each \code{x} ordinate}
}

\author{David Lucy <\email{d.j.lucy@bradford.ac.uk}> 
\cr Robert Aykroyd <\email{robert@amsta.leeds.ac.uk}>\url{http://www.amsta.leeds.ac.uk/~robert/}}

\references{Robertson, I. Lucy, D. Baxter, L. Pollard, A.M. Aykroyd, R.G. Carter, A.H.C. Switsur, V.R. and Waterhouse, J.S.(1999) A kernel based Bayesian approach to climatic reconstruction. \emph{Holocene} \bold{9}(4): 495-500}

\note{\note{Slow code suitable for visualisation and display of p.d.f where highly generalised k.p.d.fs are needed - \code{\link[KernSmooth]{bkde}} is faster when uniformly grided, single bandwidth, k.p.d.fs are required, although in the single dimensional case you won't notice the difference.

This function doesn't use bins as such, it calculates the density at a set of points. These points can be thought of as 'bin centres' but in reality they're not.

From version 1.00 onwards a number of improvements have been made: NA's are now handled semi-convincingly by dropping if required. A multi-element vector of bandwidths associated with each case can be sent, so it is possible to accept the default, give a fixed bandwidth, or a bandwidth associated with each case.

It should be noted that if a multi-element vector is sent for bandwidth, it must be of the same length as the data vector. Furthermore, multi-element vectors which approximate the bin centres, can be sent rather than the extreme limits of the range; which means that the points at which the density is to be calculated need not be uniformly spaced.

If the default \code{xbandwidth} is to be used there \bold{must} be at least five unique values for in the \code{x} vector. If not the function will return an error. If you don't have five unique values in the vector then send a value, or vector for \code{xbandwidth}

The number of ordinates defaults to the length of \code{range.x} if \code{range.x} is a vector of ordinates, otherwise it is \code{xgridsize}, or 100 if that isn't specified.

Finally, the various modes of sending parameters can be mixed, ie: the extremes of the range can be sent to define the range for \code{x}, but a multi-element vector could be sent to define the ordinates in the \code{y} dimension, or, a vector could be sent to describe the bandwidth for each case in \code{x}.
}

\seealso{\code{\link{KernSur}} \code{\link{per}}  \code{\link[base]{density}} \code{\link[base]{hist}} \code{\link[KernSmooth]{bkde}} \code{\link[KernSmooth]{bkde2D}} \code{\link[KernSmooth]{dpik}}}

\section{Acknowledgements}{Written in collaboration with A.M.Pollard <\email{a.m.pollard@bradford.ac.uk}> with the financial support of the Natural Environment Research Council (NERC) grant GR3/11395}

\examples{
x <- c(2,4,6,8) 			# make up some x data
z <- KernSec(x, xbandwidth=2, range.x=c(0,10))
plot(z$xvals, z$yden, type="l") 
# use a defined vector for the ordinates and bandwidths
ords <- seq(from=0, to=10, length=100)
bands <- x/15
z <- KernSec(x, xbandwidth=bands, range.x=ords)
plot(z$xvals, z$yden, type="l") 	# should plot a wriggly line
}

\keyword{distribution}
