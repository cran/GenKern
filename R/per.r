# per.r find the value of the ith percentage point of a
# probability distribution

per <- function(den, vals, point)
	{
	# prepare the input vectors
	den <- as.double(den)
	vals <- as.double(vals)
	point <- as.double(point)
	lenden <- as.integer(length(den))
	ans <- as.double(0)

	# invoke the c module
	out <- .C("per", den, lenden, vals, point, ans)

	# assign output
	pong <- out[[5]]
	return(pong)
	}

