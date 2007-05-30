\name{per}

\alias{per}

\title{Locate value for ith percentage point in a binned distribution}

\description{Calculates the value for the ith point in a binned distribution}

\usage{per(den, vals, point, na.rm=FALSE, neg.rm=FALSE)}

\arguments{
\item{den}{vector of frequency or density values}
\item{vals}{vector of values corresponding to the centres of the bins in \code{den}, or the bin break points}
\item{point}{percentage point of the distribution ie: 0.50 is median}
\item{na.rm}{behaviour for NA's in the vector of density values: \code{FALSE} (default) per() will fail with warning if NA's are detected, \code{TRUE} per() will assume that these values are really zeros}
\item{neg.rm}{per() will also fail if any member of the density vector is negative (which can happen occasionally from density functions based on FFT), set this to \code{TRUE} to treat these values as zeros}
}

\value{returns a value:
\item{x}{value of \code{vals} corresponding to the \code{point} position}
}

\author{David Lucy <\email{d.lucy@lancaster.ac.uk}> \url{http://www.maths.lancs.ac.uk/~lucy/}
\cr Robert Aykroyd <\email{r.g.aykroyd@leeds.ac.uk}>\url{http://www.amsta.leeds.ac.uk/~robert/}}

\note{Not restricted to uniform bin widths but due to linear interpolation gets less accurate as bin widths deviate from uniformity. The vectors must be in 
ascending order of bin centres bin break points. The density can be a frequency in that it doesn't have to sum to unity.

Out of character for the rest of the GenKern package this function does assume proper bins rather than ordinates, although if a density estimate has been generated using \code{\link{KernSec}} then the ordinate vector can be used as a first order approximation to bin centres.}

\seealso{\code{\link{KernSur}} \code{\link{per}}  \code{\link[base]{density}} \code{\link[base]{hist}} \code{\link[KernSmooth]{bkde}} \code{\link[KernSmooth]{bkde2D}} \code{\link[KernSmooth]{dpik}}}

\section{Acknowledgements}{Written in collaboration with A.M.Pollard <\email{mark.pollard@rlaha.ox.ac.uk}> with the financial support of the Natural Environment Research Council (NERC) grant GR3/11395}

\examples{
# make up some x-y data
x <- seq(1,100)
y <- dnorm(x, mean=40, sd=10)
plot(x,y)                   
# mark the median, 0.1 and 0.9 positions with vertical lines
abline(v=per(y,x,0.5))
abline(v=per(y,x,0.9))
abline(v=per(y,x,0.1)) 
# for a bimodal distribution which doesn't sum to one
x <- c(1:5) 
y <- c(2,3,4,3,4)
per(y,x,0.5) # should return 3.25
# change the previous example to bin extremes
x <- c(1:6)
per(y,x,0.5) # should return 3.75
}

\keyword{nonparametric}
