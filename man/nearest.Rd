\name{nearest}

\alias{nearest}

\title{Index of a vector nearest in value to a supplied value}

\description{Returns the index of a vector which contains the value closest to an arbitary value}

\usage{nearest(x, xval, outside=FALSE, na.rm=FALSE)}

\arguments{
\item{x}{vector of values}
\item{xval}{value to find the nearest value in \code{x} to}
\item{outside}{if not set to \code{TRUE} the function returns an error if \code{xval} is outside the range of \code{x} - default \code{FALSE}}
\item{na.rm}{NA behaviour: \code{TRUE} drops cases with NA's, \code{FALSE} stops function with a warning if NA's are detected: default=\code{FALSE}}
}

\value{returns an integer:
\item{index}{the index of \code{x} with the value nearest to \code{xval}}
}

\author{David Lucy <\email{dlucy@maths.ed.ac.uk}>\url{http://www.maths.ed.ac.uk/~dlucy/}
\cr Robert Aykroyd <\email{robert@amsta.leeds.ac.uk}>\url{http://www.amsta.leeds.ac.uk/~robert/}}

\note{The vector doesn't have to be in any particular order - this routine will just give the index of the nearest number. The only inconsistancy is that if the value of \code{xval} are not strictly within the range of the vector the function will return an error. To prevent this call with the \code{outside=TRUE} flag enabled. If there are many values which match the 'nearest' value then the function will return a vector of their indicies.}

\section{Acknowledgements}{Written in collaboration with A.M.Pollard <\email{a.m.pollard@bradford.ac.uk}> with the financial support of the Natural Environment Research Council (NERC) grant GR3/11395}

\examples{
# make up a vector
x <- c(1,2,2,2,2,2,3,4,5,6,7,8,9,10)
# conventional useage - xval within range should return 9
nearest(x, 4.7)
# xval - outside the range of x should return 14
nearest(x, 12.7, outside=TRUE)
# many 'nearest' values in x - should return - 2 3 4 5 6 
nearest(x, 1.7)
# make x[3] an NA
x[3] <- NA
# returns - 2 4 5 6 - by enabling na.rm
nearest(x, 1.7, na.rm=TRUE)
}

\keyword{arith}
