/* mathSyminv2() --- compute inverse of a symmetric matrix */

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

/* mathSyminv2() obtains the inverse of a symmetric matrix a[][] of order n.
   It is a translation from Fortran-II to ANSI-C of a subroutine named
   SMINV which D.Wells has used since 1963. SMINV was/is a translation
   from Algol-60 to Fortran-II of CACM Algorithm 150, named "syminv2".
   That name plus the comments, some variable names, and some of the
   detailed (C-like) logic from the original Algol-60 have been restored
   in this ANSI-C version of the algorithm. Also, the first dimension
   has been added to the arguments in order to invert partially-filled
   matricies. D.Wells, NRAO-CV, 1995-08-01. */

#include "ray.h"
#define A(I,J) a[I*dim2+J]

int mathSyminv2(                   /* returns TRUE if matrix is singular    */
		double a[],        /* input a[n][n], overwritten by inverse */
		int n,             /* order of 2-D matrix                   */
		int dim2)          /* size of fastest changing axis         */
{
    double       *p, *q, bigajj;
    int          *r, i, j, k, fail=FALSE;

    if ((q = (double *) listAlloc(n, sizeof(double))) == NULL) {
	printf("mathSyminv2: could not malloc q[%d]!?!\n", n);
	exit(EXIT_FAILURE);
    }
    if ((p = (double *) listAlloc(n, sizeof(double))) == NULL) {
	printf("mathSyminv2: could not malloc p[%d]!?!\n", n);
	exit(EXIT_FAILURE);
    }
    if ((r = (int *) listAlloc(n, sizeof(int))) == NULL) {
	printf("mathSyminv2: could not malloc r[%d]!?!\n", n);
	exit(EXIT_FAILURE);
    }
    for (i = 0; i < n; i++) r[i] = TRUE;
    for (i = 0; i < n; i++) {                         /* grand loop */
	for (j = 0, bigajj = 0.0; j < n; j++)         /* search for pivot */
	    if (r[j] && (fabs(A(j,j)) > bigajj)) {
		bigajj = fabs(A(j,j));
		k = j;
	    }
	if (bigajj == 0.0) {
	    printf("mathSyminv2: matrix a[%d][%d] is singular!?!\n", n, n);
	    fail = TRUE;
	    break;
	}
	r[k] = FALSE;               /* preparation of elimination step i */
	q[k] = 1.0 / A(k,k);
	p[k] = 1.0;
	A(k,k) = 0.0;
	for (j = 0; j < k; j++) {
	    p[j] = A(j,k);
	    q[j] = (r[j] ? -A(j,k) : +A(j,k)) * q[k];
	    A(j,k) = 0.0;
	}
	for (j = k+1; j < n; j++) {
	    p[j] = r[j] ? +A(k,j) : -A(k,j);
	    q[j] = -A(k,j) * q[k];
	    A(k,j) = 0.0;
	} 
	for (j = 0; j < n; j++)                   /* elimination proper */
	    for (k = j; k < n; k++)
		A(j,k) += p[j] * q[k];
    }
    if (!fail)                              /* copy to lower triangular */
	for (k = 0; k < (n - 1); k++)
	    for (j = (k + 1); j < n; j++)
		A(j,k) = A(k,j);
    listFree(p, "m.S: p");
    listFree(q, "m.S: q");
    listFree(r, "m.S: r");
    return(fail);
}
