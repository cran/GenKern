\name{per}

\alias{per}

\title{Locate value for ith percentage point in a binned probability distribution}

\description{Calculates the value for the ith point in a binned probability distribution}

\usage{per(den, vals, point)}

\arguments{
\item{den}{vector of density values}
\item{vals}{vector of values corresponding to the centres of the bins in den}
\item{point}{percentage point of the distribution ie: 0.50 is median}
}

\value{returns a value:
\item{x}{value of \code{vals} corresponding to the \code{points} position}
}

\author{David Lucy <\email{d.j.lucy@bradford.ac.uk}>
\cr Robert Aykroyd <\email{robert@amsta.leeds.ac.uk}>\url{http://www.amsta.leeds.ac.uk/~robert/}}

\references{Robertson, I. Lucy, D, Baxter, L. Pollard, A.M. Aykroyd, R.G. Carter, A.H.C. Switsur, V.R. and Waterhouse, J.S.(1999) A kernel based Bayesian approach to climatic reconstruction.\emph{Holocene} \bold{9}(4): 495-500}

\note{not limited to uniform bin widths but due to linear interpolation gets less accurate}

\seealso{\code{\link{KernSur}} \code{\link{KernSec}}}

\section{Acknowledgements}{Written in collaboration with A.M.Pollard <\email{a.m.pollard@bradford.ac.uk}> with the financial support of the Natural Environment Research Council (NERC) grant GR3/11395}

\examples{
# make up some x-y data
x <- seq(1,100)
y <- dnorm(x, mean=40, sd=10)
plot(x,y)                   
#mark the median, 0.1 and 0.9 positions with vertical lines
abline(v=per(y,x,0.5))
abline(v=per(y,x,0.9))
abline(v=per(y,x,0.1)) 
}

\keyword{density}
