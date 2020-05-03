/* mathZernike() -- recursive calculation of the Zernike polynomials. 
   D.Wells, NRAO-CV.
   1995-07-28: first version in C
   1995-08-04: changes
   1997-10-15: ZINDEX() changed to two-digit order numbers (n,m<=99) 
*/

/* -=-=-=-=-=-=-=-=- Begin Copyright Notice -=-=-=-=-=-=-=-=-=-=-=-=-=-

   Copyright (C) 1995 Associated Universities, Inc. Washington DC, USA.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

   Correspondence concerning GBT software should be addressed as follows:
   GBT Operations, National Radio Astronomy Observatory, P. O. Box 2,
   Green Bank, WV 24944-0002 USA

   -=-=-=-=-=-=-=-=-=- End Copyright Notice -=-=-=-=-=-=-=-=-=-=-=-=-=- */

/* The Zernike polynomials are the complete set of orthogonal polynomials
   on the unit circle.  See Section 9.2.1 and Appendix VII of Born & Wolf
   (1959).  This function is a transliteration from Fortran-66 to ANSI-C of
   a subroutine which was originally coded circa 1974 by Larry E. Goad (now
   at Itek Optical System) while he was at Kitt Peak National
   Observatory. Argument nord[] gives the maximum radial order n>=0 for
   each angular order m, up to an angular order for which
   nord[]=-1. Table-XXI on p.464 of B&W can be reproduced with
   nord[]={8,8,8,8,8,8,8,8,-1} and the primary (Seidel) abberation terms in
   Table XXII on p.469 can be computed with nord[]={4,3,2,-1}. If mathZernike()
   is called with u[]==NULL, it returns only the #coeffs as the function
   value, so that arrays can be malloc-ed. The ui[] array returned by
   mathZernike() contains integers which should be printed as five decimal
   digits using format "%05d"; the digits are the radial order, angular
   order and sin/cos indicator (cos=0). */

#include "ray.h"
#define ZINDEX(n,m,i) (((2 * n + m) * 100 + m) * 10 + i)

int mathZernike(                       /* returns #coeffs in u[] & ui[]    */
		double r,              /* radius in range 0-->1.0          */
		double cost,           /* cosine of theta                  */
		double sint,           /* sine   of theta                  */
		int    nord[],         /* order limits                     */
		double u[],            /* Zernike coefficients return here */
		int    ui[])           /* Indicies of coeffs return here   */
{
    double cssn[2], x, km0, k2m0, k2m2, m2, p0, p1, p, q;
    int nu, m, nmax, kmax, ncssn, k, i;

    x = 2.0 * r * r - 1.0;
    cssn[0] = 1.0;
    cssn[1] = 0.0;
    m = 0;
    nu = 0;
    if (u != NULL) {
	u[nu] = 1.0;
	ui[nu] = ZINDEX(0,0,0);
    }

    while ((nmax = nord[m]) >= 0) {        /* loop on m (angular order): */
	if (nmax > 0) {
	    kmax = (nmax - m) / 2;
	    ncssn = (m == 0) ? 1 : 2;
	    k = 0;                	        /* perform the k=0 case: */
	    p0 = 1.0;
	    if (m > 0)
		for (i = 0; i < ncssn; i++) {
		    nu++;
		    if (u != NULL) {
			u[nu] = p0 * cssn[i];
			ui[nu] = ZINDEX(k,m,i);
		    }
		}
	    if (k < kmax) {
		k = 1;                 		/* perform the k=1 case: */
		p1 = 0.5 * ((m + 2.0) * x - m);
		for (i = 0; i < ncssn; i++) {
		    nu++;
		    if (u != NULL) {
			u[nu] = p1 * cssn[i];
			ui[nu] = ZINDEX(k,m,i);
		    }
		}
		if (k < kmax) {
		    km0 = k + m;    /* initialize radial recursion loop: */
		    k2m0 = 2 * k + m;
		    k2m2 = 2 * k + m + 2;
		    m2 = m * m;
		    do {     /* now perform radial recursion for k-2..n: */
			p = (-k * km0 * k2m2 * p0) +
			    (0.5 * (k2m0 + 1) * (k2m0 * k2m2 * x - m2) * p1);
			p0 = p1;
			km0 += 1.0;
			k++;
			p1 = p / (k * km0 * k2m0);
			k2m0 = k2m2;
			k2m2 += 2.0;
			for (i = 0; i < ncssn; i++) {
			    nu++;
			    if (u != NULL) {
				u[nu] = p1 * cssn[i];
				ui[nu] = ZINDEX(k,m,i);
			    }
			}
		    } while (k < kmax);
		}
	    }
	}
	m++;	                             /* recursion for m*theta: */
	q       = (cssn[0] * cost - cssn[1] * sint) * r;
	cssn[1] = (cssn[1] * cost + cssn[0] * sint) * r;
	cssn[0] = q;
    }
    return(++nu);
}
