KernSur <- function(
		 x, 
		 y, 
		 xgridsize=100, 
		 ygridsize=100,	 
		 correlation, 
		 xbandwidth, 
		 ybandwidth, 
		 range.x,
		 range.y 
		 ) 

{

# multiplier for the xbandwidth which gives the tails at either end for the 
# final correlated kernel density estimate only operates as default 
deadspace <- 1.5        
 
# set up the default values for a few required variables 
 
        if(missing(correlation))
                {correlation <- cor(x,y)} 
 
        if(missing(xbandwidth))
                {xbandwidth <- dpik(x)} 
 
        if(missing(ybandwidth)) 
                {ybandwidth <- dpik(y)} 
 
        if(missing(range.x)) 
                {range.x <- c( (min(x) - (deadspace * xbandwidth)),(max(x) + (deadspace * xbandwidth)) )} 

        if(missing(range.y)) 
                {range.y <- c( (min(y) - (deadspace * ybandwidth)),(max(y) + (deadspace * ybandwidth)) )} 

	# trap any correlations of 1 as they will cause division by zero errors later
	if(correlation > 0.999)
		{correlation <- 0.999}

	if(correlation < -0.999)
		{correlation <- -0.999}


# if not default extract the ranges
min.x <- range.x[1] 
max.x <- range.x[2] 

min.y <- range.y[1] 
max.y <- range.y[2] 

binsinx <- round(xgridsize) 
xperbin <- (max.x - min.x)/binsinx 

binsiny <- round(ygridsize) 
yperbin <- (max.y - min.y)/binsiny 

#define the  kernelsurface array
corker <- rep(0, (binsinx * binsiny)) 
 
# generate vectors of x and y values which identify the x's and y values of each
# bin in corker where we're bound to get picked up is what bit of the bin the value
# refers to the minx is the left of bin 1 - maxx is the right of the final bin - therefore 
# we have some vernier type scaling of the central values for each bin and need 
# to adjust each accordingly 
xvals <- seq(min.x, max.x, length=binsinx) 
mangle <- seq(0.5 , -0.5 , length=binsinx) 
xvals <- xvals + (xperbin * mangle) 

yvals <- seq(min.y, max.y, length=binsiny) 
mangle <- seq(0.5 , -0.5 , length=binsiny) 
yvals <- yvals + (yperbin * mangle)


# invoke the .c module 
out <- .C( 
	"CorKernSur",
	 as.double(corker), 
	 as.integer(binsinx), 
	 as.integer(binsiny),
	 as.double(x), 
	 as.double(y), 
	 as.double(xvals), 
	 as.double(yvals), 
	 as.double(xbandwidth), 
	 as.double(ybandwidth), 
	 as.double(correlation), 
	 as.integer(length(x))
	 ) 
 
# assign the return values 
yden <- out[[1]] 

dim(yden) <- c(xgridsize, ygridsize) 
 
return(xvals, yvals, yden) 
}

