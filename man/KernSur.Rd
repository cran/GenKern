\name{KernSur}

\alias{KernSur}

\title{Generalised bi-variate binned kernel density estimate}

\description{Compute generalised bivariate binned kernel density estimate}

\usage{KernSur(x, y, xgridsize=100, ygridsize=100, correlation=0, xbandwidth,
ybandwidth, range.x, range.y)}

\arguments{
\item{x}{vector of x values}
\item{y}{vector of y values}
\item{xgridsize}{integer for number of bins into which x is to be divided: default=100}
\item{ygridsize}{integer for number of bins into which y is to be divided: default=100}
\item{correlation}{x,y correlation: default=cor(x,y)}
\item{xbandwidth}{value of x window width: default=dpik(x)}
\item{ybandwidth}{value of y window width: default=dpik(y)}
\item{range.x}{total range of the estimate in the x dimension: default=total x range +- 1.5 *hx}
\item{range.y}{total range of the estimate in the x dimension: default=total y range +- 1.5 *hy}
}
}
\value{returns two vectors and a matrix:
\item{xvals}{vector of bin center locations for x}
\item{yvals}{vector of bin centre locations for y}
\item{yden}{matrix of density for f(x,y) with dimensions xgridsize, ygridsize}
}

\author{David Lucy <\email{d.j.lucy@bradford.ac.uk}>
\cr Robert Aykroyd <\email{robert@amsta.leeds.ac.uk}>\url{http://www.amsta.leeds.ac.uk/~robert/}}

\references{Robertson, I. Lucy, D, Baxter, L. Pollard, A.M. Aykroyd, R.G. Carter, A.H.C. Switsur, V.R. and Waterhouse, J.S.(1999) A kernel based Bayesian approach to climatic reconstruction.\emph{Holocene} \bold{9}(4): 495-500}

\note{slow code suitable for visualisation and display of correlated p.d.f - \code{\link[KernSmooth]{bkde2D}} is much faster when uncorrelated k.p.d.fs are required}

\seealso{\code{\link{KernSec}} \code{\link{per}} \code{\link[KernSmooth]{dpik}}}

\section{Acknowledgements}{Written in collaboration with A.M.Pollard <\email{a.m.pollard@bradford.ac.uk}> with the financial support of the Natural Environment Research Council (NERC) grant GR3/11395}

\examples{
x <- c(2,4,6,8) 			# make up some x-y data
y <- x 

# calculate and plot a surface with zero correlation based on above data
op <- KernSur(x,y, xgridsize=50, ygridsize=50, correlation=0, 
	xbandwidth=1, ybandwidth=1, range.x=c(0,10), range.y=c(0,10))
image(op$xvals, op$yvals, op$yden, col=terrain.colors(100), axes=TRUE)
contour(op$xvals, op$yvals, op$yden, add=TRUE)
box()

# re-calculate and re-plot the above using a 0.8 correlation
op <- KernSur(x,y, xgridsize=50, ygridsize=50, correlation=0.8,
	 xbandwidth=1, ybandwidth=1, range.x=c(0,10), range.y=c(0,10))
image(op$xvals, op$yvals, op$yden, col=terrain.colors(100), axes=TRUE)
contour(op$xvals, op$yvals, op$yden, add=TRUE)
box()

}
\keyword{density}
