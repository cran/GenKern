/*
GenKernSec.c - part of David Lucy and Robert Aykroyd's CorKern package for the R
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
GenKernSec - outputs a uni-variate kernel density function of x. Possible to use this
instead of bkde() or density() with the bonuses that one may use non uniformly spaced
points and local bandwidths. Unlike the corresponding GenKernSur(), GenKernSec() executes
rapidly enough for no discernable difference to eb noticed between it and bkde()


For Linux and Solaris systems equiped with gcc
compile with (not within R): gcc -lm CorKernSec.c -c

link to R by (again not from within R): R SHLIB -o libGenKernSec.so GenKernSec.o

invoke from R or S with:

out <- .C(
	 "GenKernSec",
	 as.double(x),
	 as.integer(length(x)),
	 as.double(xvals),
	 as.double(xbandwidth),
	 as.double(est),
	 as.integer(binsinx)
	 )

before doing so all array variables MUST be doubles from R - single integers can
be passed as int.


variable list and corresponding useage:

in the joint p.d.f equation:

x	double *x	*(x+ctr)	vector of raw x scores
xvals	double *xvals	*(xvals + ctr1)	vector of x's corresponding to bins
hx	double *hx	*hx		vector of xbandwidths

control flow:

estimate	double *estimate	*(estimate + ctr1)	p.d.f at each *xval
binsinx		int *binsinx		*binsinx		no bins in estimate
lenx		int *lenx		*lenx			no of cases


temporary:

int ctr ctr1	ctr ctr1			counters
float A B	A B				(reusable) bits of f(x)
float PI	PI				constant pi


algorthm:

As simple as it gets - just toggles through all cases calculating the summed pdf at
each point specified in *xvals

math.h is required for the sqrt() and pow() functions

The output doesn't sum to unity at the moment
*/

#include <math.h>

void GenKernSec(
       double *x,
       int *lenx,
       double *xvals, 
       double *hx,
       double *estimate, 
       int *lenestimate
       )
{

int ctr=0, ctr1=0;
float PI=3.141593, A=0.0, B=0.0;


for(ctr = 0; ctr < *lenx; ctr++)
	{

	A = pow(*(hx + ctr), 2) * 2;

	for(ctr1 = 0; ctr1 < *lenestimate; ctr1++)
		{

		B = 1 / (sqrt(PI * A));
		
		*(estimate + ctr1) += B * exp( -(pow( (*(x + ctr) - *(xvals + ctr1)),2)/A));

		}
	}
return;
}


