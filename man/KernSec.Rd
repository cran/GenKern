\name{KernSec}

\alias{KernSec}

\title{Univariate kernel density estimate}

\description{Computes univariate kernel density estimate using Gaussian kernels which can also use non-equally spaced ordinates and adaptive bandwidths and local bandwidths}

\usage{KernSec(x, xgridsize=100, xbandwidth, range.x)}

\arguments{
\item{x}{vector of \code{x} values}
\item{xgridsize}{integer for number of ordinates at which to calculate the smoothed estimate: default=100}
\item{xbandwidth}{value of \code{x} window width, or vector of local window widths, one for each \code{x}, or one for each \code{range.x}, or a vector of length \code{xgridsize}: default=\code{dpik(x)}}
\item{range.x}{total range of the estimate in the \code{x} dimension, or a vector giving the \code{x} ordinates: default=\code{range +- 1.5 * mean bandwidth}}
}

\value{returns two vectors:
\item{xords}{vector of ordinates}
\item{yden}{vector of density estimates corresponding to each \code{x} ordinate}
}

\author{David Lucy <\email{d.lucy@lancaster.ac.uk}> \url{http://www.maths.lancs.ac.uk/~lucy/}
\cr Robert Aykroyd <\email{r.g.aykroyd@leeds.ac.uk}>\url{http://www.amsta.leeds.ac.uk/~robert/}}

\references{Lucy, D. Aykroyd, R.G. & Pollard, A.M.(2002) Non-parametric calibration for age estimation . \emph{Applied Statistics} \bold{51}(2): 183-196}

\note{Slow code suitable for visualisation and display of p.d.f where highly generalised k.p.d.fs are needed - \code{\link[KernSmooth]{bkde}} is faster when uniformly grided, single bandwidth, k.p.d.fs are required, although in the univariate case you won't notice the difference.

This function doesn't use bins as such, it calculates the density at a set of points. These points can be thought of as 'bin centres' but in reality they're not.

For version 1.10 on local kernel density estimates can now be sent, so that a vector of bandwidths can be send which is the same length as that of the observations. This will give a density which is has a unique bandwidth for each observation. Or a vector of bandwidths can be sent which is the same length as that of the number of bins. This will give a unique bandwidth for each ordinate, and is described in Wand & Jones (1995) \emph{Kernal Smoothing}. It is for the user to supply this vector of bandwidths, possibly with some form of \emph{pilot estimation}.

It should be noted that multi-element vectors which approximate the bin centres, can be sent rather than the extreme limits of the range; which means that the points at which the density is to be calculated need not be uniformly spaced.

If the default \code{xbandwidth} is to be used there \bold{must} be at least five unique values for in the \code{x} vector. If not the function will return an error. If you don't have five unique values in the vector then send a value, or vector for \code{xbandwidth}

The number of ordinates defaults to the length of \code{range.x} if \code{range.x} is a vector of ordinates, otherwise it is \code{xgridsize}, or 100 if that isn't specified.

The option \code{na.rm} is no longer supported. The function will automatically remove NAs where appropriate and possible, and will return a warning.

Finally, the various modes of sending parameters can be mixed, ie: the extremes of the range can be sent to define the range for \code{x}, but a multi-element vector could be sent to define the ordinates in the \code{y} dimension, or, a vector could be sent to describe the bandwidth for each case in \code{x}.
}


\seealso{\code{\link{KernSur}} \code{\link{per}}  \code{\link[stats]{density}} \code{\link[graphics]{hist}} \code{\link[KernSmooth]{bkde}} \code{\link[KernSmooth]{bkde2D}} \code{\link[KernSmooth]{dpik}}}

\section{Acknowledgements}{Written in collaboration with A.M.Pollard <\email{mark.pollard@rlaha.ox.ac.uk}> with the financial support of the Natural Environment Research Council (NERC) grant GR3/11395}

\examples{
x <- c(2,4,6,8,10) 

z <- KernSec(x)				# simplest invocation
plot(z$xords, z$yden, type="l")

z <- KernSec(x, xbandwidth=2, range.x=c(0,8))
plot(z$xords, z$yden, type="l")

# local bandwidths
ords <- seq(from=0, to=14, length=100)
bands <- x/15
z <- KernSec(x, xbandwidth=bands, range.x=ords)
plot(z$xords, z$yden, type="l")         # should plot a wriggly line

bands <- seq(from=1, to=4, length=100)	# improvise a pilot estimate
z <- KernSec(x, xbandwidth=bands, range.x=ords)
plot(z$xords, z$yden, type="l") 
}

\keyword{distribution}
\keyword{smooth}
