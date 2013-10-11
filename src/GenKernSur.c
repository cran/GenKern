/*
GenKernSur.c - part of David Lucy and Robert Aykroyd's GenKern package for the R
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
GenKernSur - fifth attempt to produce a reasonable correlated kernel routine
this one gets passed the values for x and y then calculates the z values from
the bi-variate density function rather than messing about with producing
kernels which are passed to the function

GebKernSur calculates the joint probability density function P(x,y) for all
values of x and y within the specified x and y ranges - use for visual examination
of these pdfs rather than any specific calculation - a useful test of the routine is to
sink the correlation to zero then compare the binned values to bkde2D() values from 
the KernSmooth() package - if ALL the passed parameters are the same this routine gives
very nearly the same values in each bin (any difference is of an order of magnetude
where it can be attributed to rounding errors). Don't use unless correlated kernels
are required as bkde2D is much quicker. The advantage with this version is that there are
no restrictions about equal sized bins, you can define bin centres anywhere. No restrictions
on variable window widths in either dmension - you can even send the directions of the 
kernels as a variable.


invoke from R or S with:
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
	 as.double(correlationsq),
	 as.integer(length(x))
	 ) 
 
before doing so all array variables MUST be doubles from R - single integers can
be passed as int - also a dyn.load("GenKerSur.so") never goes amiss.


variable list and corresponding useage:

in the joint p.d.f equation:
mux	double *x	*(x + ctr)	vector of raw x scores
x	double *x	*(x + ctr2)	vector of raw x scores
y	double *y	*(y + ctr2)	vector of raw y scores
xvals	double *xvals	*(xvals + ctr1)	vector of x values corresponding to the x bins
yvals	double *yvals	*(yvals + ctr)  vector of x values corresponding to the x bins
rho	double *rho	*(rho + ctr2)	vector of x-y correlation coefficients
rhosq	double *rhosq	*(rhosq + ctr2)	vector of squared x-y correlation coefficients
hx	double *hx	*(hx + ctr2)	vector of xbandwidths
hy	double *hy	*(hy + ctr2)	vector of ybandwidths

control flow:
estimate	double *estimate	*(estimate + ctr1)	p.d.f at y = yvalue
lenxest		int *lenxest		*lenxest		no of bins in x dimension
lenyest		int *lenyest		*lenyest		no of bins in y dimension
lenx		int *lenx		*lenx			no of xy pairs

temporary:
int ctr ctr1	ctr ctr1			counters
int ctr2	ctr2				counter
float A B	A B				(reusable) bits of f(x|y)
float PI	PI				constant pi
float rhosq	rhosq				correlation ^ 2 (reuseable)


algorthm:

we're trying to calculate a bi-variate correlated kernel surface using
Gaussian kernels - all we do is get an array of x values, an array of which value
each bin corresponds to and a y value from which to calculate corresponding densities
which is done directly

the p.d.f. is somewhat long so to make it managable I cut it into chunks and assigned
each to a variable - there are some repeated units so was also used for optimisation

math.h is required for the sqrt() and pow() functions

The results from here are standardised so the volume is unity if you don't want that
then you'll have to multiply the output by some factor external to this function
*/

#include<math.h>

void GenKernSur(
          double *estimate,
	  int *lenxest,
	  int *lenyest,
	  double *x,
	  double *y,
	  double *xvals, 
	  double *yvals,
	  double *hx,
	  double *hy,
	  double *rho,
	  double *rhosq,
	  int *lenx
	  )
{

int ctr=0, ctr1=0, ctr2=0;
float PI=3.141593, A=0.0, B=0.0, sum=0.0;

 
/* ctr2 for all pairs in x and y */
for(ctr2 = 0; ctr2 < *lenx; ctr2++)
	{

	/* ctr for all values of y */
	for(ctr = 0; ctr < *lenyest; ctr++)
		{
		/* ctr1 for all values of x */
		for(ctr1=0; ctr1 < *lenxest; ctr1++)
			{
			A = ( *(xvals + ctr1) - *(x + ctr2) ) / *(hx + ctr2); 
			B = ( *(yvals + ctr) - *(y + ctr2)) / *(hy + ctr2); 

			/* calculate the value in the bin */
			*(estimate + ((ctr * *lenxest) + ctr1)) += 1 / (*(hy + ctr2)* *(hx + ctr2) * (sqrt(2 * PI * (1 - *(rhosq + ctr2))))) 
				      	* 
				       	exp 
				       	(  				        
				       	(-1 / (2 * (1 - *(rhosq + ctr2))))  
				       	* 
				       	(pow(A, 2) - 2 * *(rho + ctr2) * A * B + pow(B, 2))

					) ;
			}
		}
	}

/* do not do this as we want the density function so we don't get
confused when we try to use non-uniform bins
/* make the volume equal unity */

/*
for(ctr=0; ctr<(*lenxest * *lenyest); ctr++)
	{
	sum += *(estimate + ctr);
	}

for(ctr=0; ctr<(*lenxest * *lenyest); ctr++)
	{
	*(estimate + ctr) /= sum;
	}
*/

return;
}
