/* rayGetFoci() --- function to locate the foci of bundles of converging rays
   D.Wells, NRAO-CV, Nov93/July95,95-12-12 */

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

#include "ray.h"

struct Node *rayGetFoci(			   /* returns list of Focus */
			struct Node *RayBundleSet, /* list of lists of rays */
			double tolerance,	   /* position accuracy     */
			struct Node *Segments)	   /* appends to Segments   */
{
	const int debug = 0;
	const double radius_factor = 50.0; /* rms multiplier for wavefront radius */
	struct Node *set, *bundle, *foci, *p;
	struct Ray *R;
	struct Focus *focus;
	struct Segment *seg;
	int i, j, k, k0, k1, k2, m;
	double r[4], left[4], right, weight, c[4][4], o_minus_c, x[4], Isum, size, sphere, radius;
	union ad1d2 {
		double d1[16];
		double d2[4][4];
	} a;

	set = RayBundleSet;
	if (set == NULL) {
		printf("rayGetFoci: RayBundleSet is NULL!?!\n");
		exit(EXIT_FAILURE);
	}
	listNodeCheck(set, "r.G.F: set");
	foci = listInitialize(set->item);

	while (set->next->next != set->next) { /* loop over the set of bundles: */
		set = set->next;
		bundle = (struct Node *)set->item;
		listNodeCheck(bundle, "r.G.F: bundle");
		if ((focus = (struct Focus *)listAlloc(1, sizeof(struct Focus))) == NULL) {
			printf("rayGetFoci: failed to listAlloc focus struct!?!\n");
			exit(EXIT_FAILURE);
		}
		strcpy(focus->name, bundle->item);

		/* Accumulate normal eqns for this bundle: */
		for (i = 0; i < 4; i++) { /* first zero the eqns */
			r[i] = 0.0;
			for (j = 0; j < 4; j++)
				a.d2[i][j] = 0.0;
		}
		focus->n = 0;
		p = bundle;
		while (p->next->next != p->next) { /* loop over rays in bundle */
			p = p->next;
			R = p->item;
			(focus->n)++;
			for (k = 0; k < 3; k++) { /* loop over xyz */
				for (i = 0; i < 4; i++)
					left[i] = 0.0;
				left[k] = 1.0; /* define this normal eqn */
				left[3] = -R->Q[k];
				right = (R->T[k] - R->Q[k] * R->Length);
				weight = R->Intensity;
				for (i = 0; i < 4; i++) { /* accumulate normal eqn */
					for (j = i; j < 4; j++)
						a.d2[i][j] += (left[i] * left[j]) * weight;
					r[i] += (right * left[i]) * weight;
				}
			}
		}
		for (i = 0; i < 4; i++) /* copy to lower triangular */
			for (j = i; j < 4; j++)
				a.d2[j][i] = a.d2[i][j];
		if (debug)
			printf("# rayGetFoci correlation matrix c[4][4]:\n");
		for (i = 0; i < 4; i++) { /* compute correlation matrix */
			if (debug)
				printf("#        ");
			for (j = 0; j < 4; j++) {
				c[i][j] = a.d2[i][j] / sqrt(a.d2[i][i] * a.d2[j][j]);
				if (debug)
					printf(" %6.3f", c[i][j]);
			}
			if (debug)
				printf("\n");
		}
		if (mathSyminv2(a.d1, 4, 4)) /* get inverse of normal eqns */
			printf("rayGetFoci: mathSyminv2(a,4) returned TRUE\n");
		for (i = 0; i < 4; i++) /* get the answers */
			for (j = 0, x[i] = 0.0; j < 4; j++)
				x[i] += a.d2[i][j] * r[j]; /* x = r * a^-1 */
		for (i = 0; i < 3; i++)
			focus->xyz[i] = x[i]; /* focal point position */
		focus->lc = x[3];	      /* path length to f.p. */
		p = bundle;		      /* now get xyz rms: */
		for (i = 0; i < 3; i++)
			focus->xyzs[i] = 0.0;
		Isum = 0.0;
		p = bundle;
		while (p->next->next != p->next) { /* loop over rays in bundle */
			p = p->next;
			R = p->item;
			for (i = 0; i < 3; i++) {
				o_minus_c = (R->T[i] - R->Q[i] * R->Length) - focus->xyz[i] + (R->Q[i] * focus->lc);
				weight = R->Intensity;
				focus->xyzs[i] += o_minus_c * o_minus_c * weight;
				Isum += weight;
			}
		}
		for (i = 0; i < 3; i++)
			focus->xyzs[i] = sqrt(focus->xyzs[i] / Isum);
		/* Get wavefront(length) rms: */
		/* We will choose a sphere radius, will move back along the
		   rays by that amount, will compute the radii from the
		   computed focal point to those points on the rays, and will
		   produce the rms of difference between those radii and the
		   chosen sphere radius.  The sphere radius will be large
		   compared to the rms size of the focal point. */
		for (i = 0, size = 0.0; i < 3; i++)
			size += focus->xyzs[i] * focus->xyzs[i];
		size = sqrt(size / 3.0);
		if (size < tolerance)
			size = tolerance;
		sphere = radius_factor * size;
		focus->ls = 0.0;
		p = bundle;
		while (p->next->next != p->next) { /* loop over rays in bundle */
			p = p->next;
			R = p->item;
			for (i = 0, radius = 0.0; i < 3; i++) { /* get point on ray */
				r[i] = R->T[i] + ((focus->lc - sphere) - R->Length) * R->Q[i];
				o_minus_c = (r[i] - focus->xyz[i]);
				radius += o_minus_c * o_minus_c;
			}
			o_minus_c = sqrt(radius) - sphere; /* the wavefront error */
			weight = R->Intensity;
			focus->ls += o_minus_c * o_minus_c * weight;
		}
		focus->ls = sqrt(focus->ls / Isum); /* rms wavefront error */
		if (debug)
			printf("# rayGetFoci: sphere = %8.4f, focus->ls=%8.3g\n", sphere, focus->ls);

		listAppend(focus, foci);

		/* Now plot error box for this focal point: */
		for (k = 0; k < 3; k++) {
			k0 = (k + 0) % 3;
			k1 = (k + 1) % 3;
			k2 = (k + 2) % 3;
			m = 3;
			for (i = -m; i <= +m; i += 2 * m) {
				for (j = -m; j <= +m; j += 2 * m) {
					if ((seg = (struct Segment *)listAlloc(1, sizeof(struct Segment))) == NULL) {
						printf("rayGetFoci: listAlloc(struct Segment)=>NULL!?!\n");
						exit(EXIT_FAILURE);
					}
					seg->T1[k0] = focus->xyz[k0] - focus->xyzs[k0] * m;
					seg->T2[k0] = focus->xyz[k0] + focus->xyzs[k0] * m;
					seg->T1[k1] = seg->T2[k1] = focus->xyz[k1] + focus->xyzs[k1] * i;
					seg->T1[k2] = seg->T2[k2] = focus->xyz[k2] + focus->xyzs[k2] * j;
					seg->ColorCode = 2001;
					listAppend(seg, Segments);
					if ((seg = (struct Segment *)listAlloc(1, sizeof(struct Segment))) == NULL) {
						printf("rayGetFoci: listAlloc(struct Segment)=>NULL!?!\n");
						exit(EXIT_FAILURE);
					}
					seg->T1[k0] = seg->T2[k0] = focus->xyz[k0] + focus->xyzs[k0] * i;
					seg->T1[k1] = focus->xyz[k1] - focus->xyzs[k1] * m;
					seg->T2[k1] = focus->xyz[k1] + focus->xyzs[k1] * m;
					seg->T1[k2] = seg->T2[k2] = focus->xyz[k2] + focus->xyzs[k2] * j;
					seg->ColorCode = 2001;
					listAppend(seg, Segments);
				}
			}
		}
	}
	return (foci);
}
