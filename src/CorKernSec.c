/*
CorKernSec.c - part of David Lucy and Robert Aykroyd's CorKern package for the R
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
CorKernSec - forth attempt to produce a reasonable correlated kernel routine
this one gets passed the values for x and y then calculates the z values from
the bi-variate density function rather than messing about with producing
kernels which are passed to the function

CorKernSec outputs a section from the bi-variate density function of x and y this
section being propto P(x,y) given a specific y - use to calculate f(x|y) for Bayesian
estimates of continuous functions. Possible to use instead of bkde() or density() as
in the univariate case the time diffrence in execution is so small as to be unoticable


For Linux and Solaris systems equiped with gcc
compile with (not within R): gcc -lm CorKern.c -c

link to R by (again not from within R): R SHLIB -o libCorKern.so CorKern.o

invoke from R or S with:
out <- .C(
	 "CorKernSec",
	 as.double(corker),
	 as.integer(binsinx),
	 as.double(x),
	 as.double(xvals),
	 as.integer(length(xvals)),
	 as.double(y),
	 as.double(yvalue),
	 as.double(xbandwidth),
	 as.double(ybandwidth),
	 as.double(correlation)
	 )

before doing so all array variables MUST be doubles from R - single integers can
be passed as int - also a dyn.load("libKerAdd.so") never goes amiss in the calling
function.


variable list and corresponding useage:

in the joint p.d.f equation:

mux	double *x	*(x+ctr)	vector of raw x scores
x	double *xvals	*(xvals + ctr1)	vector of x's corresponding to bins
y	double *y	*(y + ctr)	vector of raw y scores
muy	double *yvalue	*yvalue		y value for which we're finding f(x|y)
rho	double *rho	*rho		x-y correlation coefficient
hx	double *hx	*hx		xbandwidth
hy	double *hy	*hy		ybandwidth

control flow:

estimate	double *estimate	*(estimate + ctr1)	p.d.f at y = yvalue
binsinx		int *binsinx		*binsinx		no bins in estimate
lenx		int *lenx		*lenx			no of cases


temporary:

int ctr ctr1	ctr ctr1			counters
float A B	A B				(reusable) bits of f(x|y)
float PI	PI				constant pi
float rhosq	rhosq				correlation ^ 2 (reuseable)


algorthm:

we're trying to calculate a slice of a bi-variate correlated kernel surface using
Gaussian kernels - all we do is get an array of x values, an array of which value
each bin corresponds to and a y value from which to calculate corresponding densities
which is done directly

the p.d.f. is somewhat long so to make it managable I cut it into chunks and assigned
each to a variable - there are some repeated units so was also used for optimisation

math.h is required for the sqrt() and pow() functions

The results from here will be propto the conditional, but are in fact the joint - to
get the true conditional each binned value really needs to be /n - alternatively you 
could just normailse the estimate
*/

#include <math.h>

CorKernSec(
       double *estimate, 
       int *lenestimate,
       double *x,
       double *xvals, 
       int *lenx,
       double *y,
       double *yvalue,
       double *hx,
       double *hy,
       double *rho
       )
{

int ctr=0, ctr1=0;
float rhosq=0.0, PI=3.141593, A=0.0, B=0.0;

rhosq = pow(*rho,2);

for(ctr = 0; ctr < *lenx; ctr++)
	{

	for(ctr1 = 0; ctr1 < *lenestimate; ctr1++)
		{

		A = ( *(xvals + ctr1) - *(x + ctr) ) / *hx;
		B = ( *yvalue - *(y + ctr)) / *hy;

		*(estimate + ctr1) += 1 / (*hy * *hx * (sqrt(2 * PI * (1 - rhosq))))
				      	*
				       	exp
				       	( 
				       
				       	(-1 / (2 * (1 - rhosq))) 
				       	*
				       	(pow(A, 2) - 2 * *rho * A * B + pow(B, 2))
				       
				       	);


		}
	}
return;
}


