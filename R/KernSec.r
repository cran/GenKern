KernSec <- function(
		 x,
		 y,
		 yvalue,
		 gridsize=100,
		 correlation,
		 xbandwidth,
		 ybandwidth,
		 range.x
		 )
{

# multiplier for the xbandwidth which gives the tails at either end for the
# final correlated kernel density estimate only operates as default
deadspace <- 1.5       


# set up the default values for a few required variables

        if(missing(correlation))                                                                                if(missing(correlation))
                {correlation <- cor(x,y)}

	if(missing(xbandwidth))                                                                                        if(missing(xbandwidth))
                {xbandwidth <- dpik(x)}

        if(missing(ybandwidth))
                {ybandwidth <- dpik(y)}

        if(missing(range.x))
                {range.x <- c( (min(x) - (deadspace * xbandwidth)),(max(x) + (deadspace * xbandwidth)) )}
            
# if not default extract the range
min.x <- range.x[1]
max.x <- range.x[2]


binsinx <- round(gridsize)
xperbin <- (max.x - min.x)/binsinx

# define the 
corker <- rep(0, binsinx)


# generate a vector of x values which identify the x value of each bin in corker
# where we're bound to get picked up is what bit of the bin the value refers to
# the minx is the left of bin 1 - maxx is the right of the final bin - therefore
# we have some vernier type scaling of the central values for each bin and need
# to adjust each accordingly
xvals <- seq(min.x, max.x, length=binsinx)
mangle <- seq(0.5 , -0.5 , length=binsinx)
xvals <- xvals + (xperbin * mangle)


# invoke the .c module
out <- .C(
	 "CorKernSec",
	 as.double(corker),
	 as.integer(binsinx),
	 as.double(x),
	 as.double(xvals),
	 as.integer(length(x)),
	 as.double(y),
	 as.double(yvalue),
	 as.double(xbandwidth),
	 as.double(ybandwidth),
	 as.double(correlation)
	 )

# assign the return values
yden <- out[[1]]


#print(xbandwidth)
#print(ybandwidth)


return(xvals, yden)
}

