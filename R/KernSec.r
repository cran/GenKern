KernSec <- function(
		    x,
		    xgridsize=100,
		    xbandwidth,
		    range.x,
		    na.rm=FALSE
		   )
{
# multiplier for the xbandwidth which gives the tails at either end for the
# final correlated kernel density estimate only operates as default
deadspace <- 1.5       
cases <- length(x)

# flags to increase efficiency in parameter selection
flag1 <- 1
flag2 <- 1

# xbandwidth selection
	# default bandwidth
        if(missing(xbandwidth))
                {
		flag1 <- 0
		xbandwidth <- bandwidthselect(x, bandwidths=FALSE)
		} 
	# user supplied bandwidth 
	if(flag1){xbandwidth <- bandwidthselect(x, xbandwidth)}


# Do the NA handling
# put it all together into a data frame or na.omit doesn't work
z <- data.frame(x, xbandwidth)
	# if NAs not allowed fail the function
	if(na.rm == FALSE){na.fail(z)}
	# get rid of NA cases
	if(na.rm == TRUE){z <- na.omit(z)}
# reassign the vectors with NAs removed
x <- z$x; xbandwidth <- z$xbandwidth


# range selection
	# default range of xvalues
        if(missing(range.x))
                {
		flag2 <- 0
		xvals <- rangeselect(x, rnge=FALSE, xgridsize, xbandwidth, deadspace)
		} 
	# user supplied ranges
	if(flag2){xvals <- rangeselect(x, range.x, xgridsize, xbandwidth, deadspace)}


# generate a vector of length xords with zeros in it to contain
# the density estimate
xordslen <- length(xvals)
est <- rep(0, xordslen)


# invoke the .c module
out <- .C(
	 "GenKernSec",
	 as.double(x),
	 as.integer(length(x)),
	 as.double(xvals),
	 as.double(xbandwidth),
	 as.double(est),
	 as.integer(xordslen)
	 )

# assign the return values
yden <- out[[5]]

return(xvals, yden)
}

