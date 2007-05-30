\name{KernSur}

\alias{KernSur}

\title{Bivariate kernel density estimation}

\description{Compute bivariate kernel density estimate using five parameter Gaussian kernels which can also use non equally spaced and adaptive bandwidths}

\usage{KernSur(x, y, xgridsize=100, ygridsize=100, correlation, xbandwidth,
ybandwidth, range.x, range.y, na.rm=FALSE)}

\arguments{
\item{x}{vector of \code{x} values}
\item{y}{vector of \code{y} values}
\item{xgridsize}{integer for number of ordinates at which to calculate the smoothed estimate: default=100}
\item{ygridsize}{integer for number of ordinates at which to calculate the smoothed estimate: default=100}
\item{correlation}{\code{x,y} correlation, or vector of local correlations: default=\code{cor(x,y)}}
\item{xbandwidth}{value of \code{x} window width, or vector of local window widths: default=\code{dpik(x)}}
\item{ybandwidth}{value of \code{y} window width, or vector of local window widths: default=dpik(y)}
\item{range.x}{total range of the estimate in the \code{x} dimension, or a vector giving the \code{x} ordinates: default=\code{range +- 1.5 * mean bandwidth}}
\item{range.y}{total range of the estimate in the \code{y} dimension, or a vector giving the \code{y} ordinates: default=\code{range +- 1.5 * mean bandwidth}}
\item{na.rm}{NA behaviour: \code{TRUE} drops cases with NA's, \code{FALSE} stops function with a warning if NA's are detected: default=\code{FALSE}}
}
}
\value{returns two vectors and a matrix:
\item{xords}{vector of ordinates at which the density has been estimated in the \code{x} dimension}
\item{yords}{vector of ordinates at which the density has been estimated in the \code{y} dimension}
\item{zden}{matrix of density for \code{f(x,y)} with dimensions \code{xgridsize}, \code{ygridsize}}
}

\author{David Lucy <\email{d.lucy@lancaster.ac.uk}> \url{http://www.maths.lancs.ac.uk/~lucy/}
\cr Robert Aykroyd <\email{r.g.aykroyd@leeds.ac.uk}>\url{http://www.amsta.leeds.ac.uk/~robert/}}

\references{Lucy, D. Aykroyd, R.G. & Pollard, A.M.(2002) Non-parametric calibration for age estimation. \emph{Applied Statistics} \bold{51}(2): 183-196}

\note{Slow code suitable for visualisation and display of correlated p.d.f, where highly generalised k.p.d.fs are needed - \code{\link[KernSmooth]{bkde2D}} is much faster when uncorrelated, uniformly grided, single bandwidth, k.p.d.fs are required.

This function doesn't use bins as such, it calculates the density at a set of points in each dimension. These points can be thought of as 'bin centres' but in reality they're not.

From version 1.00 onwards a number of improvements have been made: NA's are now handled semi-convincingly by dropping if required. A multi-element vector of bandwidths associated with each case can be sent for either dimension, so it is possible to accept the default, give a fixed bandwidth, or a bandwidth associated with each case. A multi-element vector of correlations can be sent, rather than a single correlation.

It should be noted that if a vector is sent for correlation, or either bandwidth, they must be of the same length as the data vectors. Furthermore, vectors which approximate the bin centres, can be sent rather than the extreme limits in the range; which means that the points at which the density is to be calculated need not be uniformly spaced.

Unlike \code{\link{KernSec}} this function does not yet support local bandwidths.

If the default \code{bandwidth} is to be used there \bold{must} be at least five unique values for in the \code{x} and \code{y} vectors. If not the function will return an error. If you don't have five unique values in the vector then send a value, or vector for \code{bandwidth}

The number of ordinates defaults to the length of \code{range.x} if \code{range.x} is a vector of ordinates, otherwise it is \code{xgridsize}, or 100 if that isn't specified.

Finally, the various modes of sending parameters can be mixed, ie: the extremes of the range can be sent to define the range for \code{x}, but a multi-element vector could be sent to define the ordinates in the \code{y} dimension, or, a vector could be sent to describe the bandwidth for each case in the \code{x} direction, and a single-element vector defines all bandwidths in the \code{y}.

Version 1.1-0 has a bugfix in that it now outputs the magnetude of the density function at the specified bi-variate points, not an approximation to the volumes.}

\seealso{\code{\link{KernSec}} \code{\link{per}}  \code{\link[base]{density}} \code{\link[base]{hist}} \code{\link[KernSmooth]{bkde}} \code{\link[KernSmooth]{bkde2D}} \code{\link[KernSmooth]{dpik}}}

\section{Acknowledgements}{Written in collaboration with A.M.Pollard <\email{mark.pollard@rlaha.ox.ac.uk}> with the financial support of the Natural Environment Research Council (NERC) grant GR3/11395}

\examples{
x <- c(2,4,6,8,10)                         # make up some x-y data
y <- x 

# calculate and plot a surface with zero correlation based on above data
op <- KernSur(x,y, xgridsize=50, ygridsize=50, correlation=0, 
        xbandwidth=1, ybandwidth=1, range.x=c(0,13), range.y=c(0,13))
image(op$xords, op$yords, op$zden, col=terrain.colors(100), axes=TRUE)
contour(op$xords, op$yords, op$zden, add=TRUE)
box()

# re-calculate and re-plot the above using a 0.8 correlation
op <- KernSur(x,y, xgridsize=50, ygridsize=50, correlation=0.8,
         xbandwidth=1, ybandwidth=1, range.x=c(0,13), range.y=c(0,13))
image(op$xords, op$yords, op$zden, col=terrain.colors(100), axes=TRUE)
contour(op$xords, op$yords, op$zden, add=TRUE)
box()

# calculate and plot a surface of the above data with an ascending
# correlation and bandwidths and a vector of equally spaced ordinates
bands <- c(1,1.1,1.2,1.3,1.0)
cors <- c(0,-0.2,-0.4,-0.6, -0.7)
rnge.x <- seq(from=0, to=13, length=100)

op <- KernSur(x,y, xgridsize=50, ygridsize=50, correlation=cors, 
        xbandwidth=bands, ybandwidth=bands, range.x=rnge.x, range.y=c(0,13))
image(op$xords, op$yords, op$zden, col=terrain.colors(100), axes=TRUE)
contour(op$xords, op$yords, op$zden, add=TRUE)
box()

}
\keyword{distribution}
\keyword{smooth}
