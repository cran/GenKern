/*
getlims.c - part of David Lucy and Robert Aykroyd's CorKern package for the R
statistical language

Copyright (C) David Lucy and Robert Aykroyd see LICENCE

David Lucy: d.j.lucy@bradford.ac.uk
Robert Aykroyd: robert@amsta.leeds.ac.uk
http://www.amsta.leeds.ac.uk/~robert/


This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software Foundation;
either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program;
if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA
*/


/*
getlims - outputs a vector of bin limits based on a vector of bin centres

For Linux and Solaris systems equiped with gcc
compile with (not within R): gcc -lm getlims.c -c

link to R by (again not from within R): R SHLIB -o libGenKerngetlims.so per.o

invoke from R or S with:
out <- .C(
	 "getlims",
	 as.double(vals),
	 as.double(limits),
	 as.integer(bins),
	 )

before doing so all array variables MUST be doubles from R - single integers can
be passed as int - also a dyn.load("libenKerngetlims.so") never goes amiss in the calling
function.


variable list and corresponding useage:

vals	double *vals	*(vals+ctr)	vector of bin centres
limits	double *limits	*(limits + ctr)	vector of bin limits
bins	int *bins	*bins		number of bins

temporary:

int ctr		ctr 				counter


algorthm:
per.c requires a vector of bin limits rather than centres so this makes the 
vector of centres coming from functions such as KernSur and KernSec, BKDE and density,
into a best guess for bin limits. Obviously there will be n+1 (where n is the number of
bins) bin limits - it is up to the user of this code to provide this vector. This routine
simply finds the half way point between bin centres, the first and last elements being
calculated on the basis of the bins being symetric about the centre. This works for
equally spaced bins, but gets worse as the bins deviate from being of uniform size, so can
only be treated as a first order aproximation for non-uniform bins.
*/
/* The bin values MUST be in ascending order */

getlims(double *vals, double *limits, int *bins)
{
int ctr;

/* assign a value to the lower edge of the first bin */
*limits = *vals - ((*(vals + 1) - *vals)/2);

/* do all the other bins */
	for(ctr=1; ctr< *bins; ctr++)
		{
		*(limits + ctr) = *(vals + ctr) - ((*(vals + ctr) - *(vals + (ctr -1))) /2);
		}

/* assign a value to the higher edge of the topmost bin */
*(limits + *bins) = *(vals + *bins - 1) + ((*(vals + *bins -1) - *(vals + *bins - 2)) /2);

return(*limits);
}
