/* rayTrace() --- trace a list of lists of Rays through a list of Surfaces
   D.Wells, NRAO-CV
   1993-10/1994-01: initial development in C
   1995-06/-07:     changes
   1997-10-15: mathVectorMatrix.h has been changed
   1997-11-20: change more calling sequences to new mathVectorMatrix style.
   1998-06-11: 'C3'-->'VC', 'V3'-->'VS', 'MV'-->'MS'
   1998-10-23: delete 'l' code on 'g'and 'f' formats
   1998-10-25: recover a redundant memory segment.
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

/* Function rayTrace() traces the set of rays RayBundleSet through the
   system of surfaces System. Intersections of rays with aspheric
   surfaces are determined to the specified tolerance.  The function
   returns a RayBundleSet containing the final positions, direction
   cosines and path_lengths of the rays. It also returns the list of
   Segments of the traced rays. */

#include "ray.h"

/* The following definition is for backward compatibiiity: */
/* Compute rotation matrix AA[3][3] from Euler(?) angles. See "VNR Concise
   Encylopedia of Mathematics'., p.535: E[0] = phi = first Euler angle
   which rotates about the X-axis, E[1]= psi = second Euler angle,
   about Y-axis, E[2] = chi = third angle about Z-axis. Positive
   angles rotate in right- hand-rule sense: if thumb points to +X,
   fingers curl toward positive phi. NOTE: the Euler(?) angle convention
   used here, as given by VNR, differs from that given by Goldstein
   and used in EULER_2_MATRIX(); eventually this difference of
   conventions must be explained properly. */
#define VNR_2_MATRIX(AA, E) ANGLES_2_MATRIX(MC, AA, X, VC(E, 0), Y, VC(E, 1), Z, VC(E, 2), NOP, NOP, NOP, NOP);

#define TOLMAX 6

struct Node *rayTrace(				 /* returns traced RayBundleSet */
		      struct Node *RayBundleSet, /* the lists of input rays */
		      struct Node *System,	 /* the list of surfaces    */
		      double tolerance,		 /* for aspheric intercepts */
		      struct Node *Segments)	 /* init this list! appends */
{
	double e, M_1x, M_1_2, xi_1, xi_1_p, L, r_1_2;
	double s_2, x_bar_0, G_0;
	double O_2, g_bar_1;

	double temp, N[3], delta, delta_length[3], index;
	double euler[3], a[3][3];
	const double t = 0.0; /* $t$ was used by Feder, is not used now (see below) */
	int j, k;
	int n_bundles, n_rays, n_NaN_bundles, n_NaN_rays, ray_reject, bundle_reject;
	int n_NaN_a, n_NaN_b, n_NaN_c;
	struct Surface *S;
	struct Ray *R, *RayIn, Rp, Rsave;
	struct Segment *Seg, *Seg2;
	struct Node *set, *bundle, *s, *FinalBundle, *FinalBundleSet;
	char finalname[NAMEMAX], *segname, traced[] = "";
	MVM_PRIVATE_VARIABLES;

	set = RayBundleSet;
	listNodeCheck(set, "r.T: RBSin");
	j = strlen(set->item);
	if (j < NAMEMAX)
		strcpy(finalname, set->item);
	else {
		strncpy(finalname, set->item, (NAMEMAX - 1));
		finalname[NAMEMAX - 1] = '\0';
	}
	j = strlen(finalname);
	k = strlen(traced);
	if ((j + k) < NAMEMAX)
		strcat(finalname, traced);
	FinalBundleSet = listInitialize(finalname);

	if (strlen(Segments->item) == 0) {
		segname = (char *)listAlloc(strlen(System->item) + strlen(finalname) + 5, sizeof(char));
		strcpy(segname, System->item);
		strcat(segname, "; ");
		strcat(segname, finalname);
		listFree(Segments->item, ""); /* 1998-10-25: recover this segment! */
		Segments->item = segname;
	}

	/* loop over the set of bundles: */
	n_bundles = n_NaN_bundles = n_rays = n_NaN_rays = n_NaN_a = n_NaN_b = n_NaN_c = 0;
	while (set->next->next != set->next) {
		set = set->next;
		listNodeCheck(set, "r.T: set");

		/* trace next bundle of rays: */
		bundle = (struct Node *)set->item;
		j = strlen(bundle->item);
		if (j < NAMEMAX)
			strcpy(finalname, bundle->item);
		else {
			strncpy(finalname, bundle->item, (NAMEMAX - 1));
			finalname[NAMEMAX - 1] = '\0';
		}
		j = strlen(finalname);
		k = strlen(traced);
		if ((j + k) < NAMEMAX)
			strcat(finalname, traced);
		FinalBundle = listInitialize(finalname);
		listAppend(FinalBundle, FinalBundleSet);
		n_bundles++;
		bundle_reject = FALSE;

		/* Loop over the input rays: */
		while (bundle->next->next != bundle->next) {
			bundle = bundle->next;
			listNodeCheck(bundle, "r.T: nextray");

			/* trace next ray: */
			RayIn = (struct Ray *)bundle->item;
			R = (struct Ray *)listAlloc(1, sizeof(struct Ray));
			for (j = 0; j < 3; j++) {
				R->T[j] = RayIn->T[j];
				R->Q[j] = RayIn->Q[j];
			}
			R->Length = RayIn->Length;
			R->Intensity = RayIn->Intensity;
			R->ColorCode = RayIn->ColorCode;
			n_rays++;
			ray_reject = FALSE;

			/* Looping over the surfaces: */
			listNodeCheck(System, "r.T: System");
			s = System;
			index = 1.0; /* assume we start with air (~vacuum) */
			while (s->next->next != s->next) {
				s = s->next;
				listNodeCheck(s, "r.T: s");
				S = (struct Surface *)s->item;

				/* Save state of ray R before tracing surface: */
				for (j = 0; j < 3; j++) {
					Rsave.T[j] = R->T[j];
					Rsave.Q[j] = R->Q[j];
				}
				Rsave.Length = R->Length;
				Rsave.Intensity = R->Intensity;
				Rsave.ColorCode = R->ColorCode;

				/* Transform ray R to coordinate system of surface S: */
				euler[0] = S->E[0]; /* first Euler angle rotates about x-axis */
				euler[1] = S->E[1]; /* second Euler angle, about y-axis */
				euler[2] = S->E[2]; /* third angle, about z-axis */
				VNR_2_MATRIX(a, euler);
				for (j = 0; j < 3; j++) { /* subtract vertex posn of surface */
					Rp.T[j] = (R->T[j] - S->S[j]);
					Rp.Q[j] = R->Q[j];
				}
				MATRIX_VECTOR_MULT(VC, R->T, MC, a, VC, Rp.T);
				MATRIX_VECTOR_MULT(VC, R->Q, MC, a, VC, Rp.Q);

				/* Get length of ray between surfaces: NOTE: variable 't' was
				   used by Feder as the vertex separation between previous
				   surface and this surface. In the new scheme, in which rays
				   are transformed to the coordinate system of the surface
				   before tracing, 't' is _zero_. It is still present in this
				   code to enable comparison with the Feder paper; the
				   optimizing compiler will eliminate it from the
				   expressions. */
				/* R->Q[0] is X, R->Q[1] is Y and R->Q[2] is Z */
				/* (X,Y,Z) is the vector along the ray to the surface */
				/* R->T[0] is x, R->T[1] is y and R->T[3] is z */
				/* (x,y,z) is the vector form of the vertex of the surface */
				/* Feder paper equation (1) */
				e = (t * R->Q[0]) - VECTOR_DOT(VC, R->T, VC, R->Q);
				/* Feder paper equation (2) */
				M_1x = R->T[0] + e * R->Q[0] - t;
				/* Feder paper equation (3) */
				M_1_2 = VECTOR_DOT(VC, R->T, VC, R->T) - (e * e) + (t * t) - (2.0 * t * R->T[0]);
				r_1_2 = 1. / (S->c_1 * S->c_1);
				if (M_1_2 > r_1_2) {
					M_1_2 = r_1_2; /* SPECIAL RULE! 96-01-22 */
				}
				/* Feder paper equation (4) */
				xi_1 = sqrt((R->Q[0] * R->Q[0]) - S->c_1 * (S->c_1 * M_1_2 - 2.0 * M_1x));
				if (isnan(xi_1)) { /* NaN! reject this ray! */
					ray_reject = bundle_reject = TRUE;
					n_NaN_rays++;
					n_NaN_a++;
					break;
				}
				/* Feder paper equation (5) */
				L = e + (S->c_1 * M_1_2 - 2.0 * M_1x) / (R->Q[0] + xi_1);

				/* Get intercept with new (spherical) surface: */
				for (j = 0; j < 3; j++)
					delta_length[j] = -R->T[j];
				VECTOR_EXTEND(VC, R->T, VC, R->Q, L);
				R->T[0] -= t;
				/* Now R->T has x1, y1, z1 */

				/* The ray has been traced to the osculating sphere with
				   curvature c1. Now we will iterate to get the intercept with
				   the nearby aspheric surface. Suppose the (rotationally
				   symmetric) aspheric is given by $$x = f(y,z)$$, and is a
				   function of $y^2 + z^2$ only. For a spherical surface, one
				   has $$x = r - (r^2 - s^2)^{1\over2}$$, where $s^2 = y^2 +
				   z^2$. For a general surface one may add deformation terms
				   to this expression and obtain $$x = c s^2 / (1 + (1 - c^2
				   s^2)^{1\over2})) + (A_2 s^2 + A_4 s^4 + ...) = f$$.  The
				   equation is expressed in this form in order to avoid
				   indeterminacy as c approaches zero, and in order to
				   represent surfaces that are nearly spherical. Near-spheres
				   cannot be handled well by a power series alone, especially
				   in the neighborhood of $s = 1 / c$.

				   In this implementation we include a term for the numerical
				   eccentricity so that we can trace any pure conic section
				   without using the $A_i$ terms. */

				j = 0;
				do {
					/* Get square of radius of intercept: */
					/* Feder equation s^2 = y^2 + z^2, section E */
					s_2 = R->T[1] * R->T[1] + R->T[2] * R->T[2];

					/* Get the point on aspheric which is at the same radius as
					   the intercept of the ray. Then compute a tangent plane to
					   the aspheric at this point and find where it intersects
					   the ray.  This point will lie very close to the aspheric
					   surface.  The first step is to compute the x-coordinate
					   on the aspheric surface using $\overline{x}_0 = f(y_0,
					   z_0)$. */
					/* (1 - c^2*s^2)^(1/2) - part of equation (12) */
					temp = sqrt(1.0 - S->c_1 * S->c_1 * s_2 * (1.0 - S->eps * S->eps));
					if (isnan(temp)) {
						ray_reject = bundle_reject = TRUE;
						n_NaN_rays++;
						n_NaN_b++;
						break;
					}
					/* Feder equation (12) */
					/* But using c*s^2/[1 + (1 - c^2*s^2)^(1/2)] + aspheric A_2*s^2 + A_4*s^4 + ... */
					x_bar_0 = (S->c_1 * s_2) / (1.0 + temp) + (S->a_2 + S->a_4 * s_2) * s_2;
					delta = fabs(R->T[0] - x_bar_0);

					/* Get the direction numbers for the normal to the
					   aspheric: */
					/* Feder equation (13), l */
					N[0] = temp;
					temp = S->c_1 + N[0] * (2.0 * S->a_2 + 4.0 * S->a_4 * s_2);
					/* Feder equation (14), m */
					N[1] = -R->T[1] * temp;
					/* Feder equation (15), n */
					N[2] = -R->T[2] * temp;

					/* Get the distance from aspheric point to ray intercept */
					G_0 = N[0] * (x_bar_0 - R->T[0]) / VECTOR_DOT(VC, R->Q, VC, N);

					/* and compute new estimate of intercept point: */
					VECTOR_EXTEND(VC, R->T, VC, R->Q, G_0);
				} while ((delta > tolerance) && (++j < TOLMAX));
				if (ray_reject)
					break;
				if (j >= TOLMAX) {
					printf("rayTrace: delta=%g, reached %d iterations!?!\n", delta, j);
					exit(EXIT_FAILURE);
				}

				for (j = 0; j < 3; j++)
					delta_length[j] += R->T[j];
				delta_length[0] += t;
				index /= S->mu_1; /* effective_distance = dl*n */
				R->Length += (VECTOR_LENGTH(VC, delta_length) * fabs(index));

				/* Now to do the refraction. First get the cosines of
				   incidence and refraction: */

				xi_1 = VECTOR_DOT(VC, R->Q, VC, N);
				O_2 = VECTOR_DOT(VC, N, VC, N);
				xi_1_p = sqrt(O_2 * (1.0 - S->mu_1 * S->mu_1) + S->mu_1 * S->mu_1 * xi_1 * xi_1);
				if (isnan(xi_1_p)) {
					ray_reject = bundle_reject = TRUE;
					n_NaN_rays++;
					n_NaN_c++;
					break;
				}
				g_bar_1 = (xi_1_p - S->mu_1 * xi_1) / O_2;
				VECTOR_EXTEND(VC, R->Q, VC, R->Q, (S->mu_1 - 1.0));
				VECTOR_EXTEND(VC, R->Q, VC, N, g_bar_1);

				/* Transform ray R back to original coordinate system: */
				for (j = 0; j < 3; j++) { /* rotate back: a^{-1} * R --> Rp */
					Rp.T[j] = Rp.Q[j] = 0.0;
					for (k = 0; k < 3; k++) {
						Rp.T[j] += (R->T[k] * a[j][k]);
						Rp.Q[j] += (R->Q[k] * a[j][k]);
					}
				}
				for (j = 0; j < 3; j++) { /* add old origin to ray position */
					R->T[j] = (Rp.T[j] + S->S[j]);
					R->Q[j] = Rp.Q[j];
				} /* R is now back in global coordinate system */

				/* Append this segment to the segment list: */
				Seg = (struct Segment *)listAlloc(1, sizeof(struct Segment));
				for (j = 0; j < 3; j++) {
					Seg->T1[j] = Rsave.T[j];
					Seg->T2[j] = R->T[j];
					Seg->ColorCode = R->ColorCode;
				}
				listAppend(Seg, Segments);
			}
			if (fabs(index) - 1.0 > 0.0001)
				printf("# rayTrace: index = %8.6f, != +/-1.0!?!\n", index);
			if (ray_reject)
				continue;
			listAppend(R, FinalBundle);
			/* Extend the last segment of this ray by 10%: */
			Seg2 = (struct Segment *)listAlloc(1, sizeof(struct Segment));
			for (j = 0; j < 3; j++) {
				Seg2->T1[j] = Seg->T2[j];
				Seg2->T2[j] = 0.10 * (Seg->T2[j] - Seg->T1[j]) + Seg2->T1[j];
				Seg2->ColorCode = Seg->ColorCode;
			}
			listAppend(Seg2, Segments);
		}
		if (bundle_reject)
			n_NaN_bundles++;
	}
	if (n_NaN_bundles != 0) {
		printf("# rayTrace: %d NaN_rays (of %d) seen in %d of %d bundles! (%d,%d,%d)\n", n_NaN_rays, n_rays,
		       n_NaN_bundles, n_bundles, n_NaN_a, n_NaN_b, n_NaN_c);
	}
	return (FinalBundleSet);
}
