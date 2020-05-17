/* rayGetPlanes() --- function to characterize nearly-plane wavefronts */

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

/*
   D.Wells, NRAO-CV
   1995-07-28: initial
   1998-10-25: fixed major bug: freeing u[] & ui[] ==> Segmentation Fault!
*/

#define SEIDEL 1
#include "ray.h"
#define A(I, J) a[I * AMNIMAX + J]

struct Node *rayGetPlanes(			     /* return list of planes */
			  struct Node *RayBundleSet, /* list of lists of rays */
			  double yz_axis[],	     /* Zernike analysis axis */
			  double yz_radius,	     /* Zernike radius        */
			  char print_mode[])	     /* "verbose"|"silent"    */
{
	struct Node *set, *bundle, *planes, *p;
	struct Ray *R;
	struct Plane *plane;
	double temp, radius, cost, sint, zlength, Isum;
	double u[AMNIMAX], a[AMNIMAX * AMNIMAX], r[AMNIMAX], c[AMNIMAX][AMNIMAX];
	int i, j, nu, ui[AMNIMAX], first_bundle;
	int nord[] = {4, 3, 2, -1}; /* max_orders for Seidel abberations */
	const char *aberration, nomatch[] = "<no_match>";

	set = RayBundleSet;
	if (set == NULL) {
		printf("rayGetPlanes: RayBundleSet is NULL!?!\n");
		exit(EXIT_FAILURE);
	}
	listNodeCheck(set, "r.G.P: set_in");
	planes = listInitialize(set->item);
	first_bundle = TRUE;

	nu = mathZernike(0.5, 1.0, 0.0, nord, NULL, NULL);
	if (nu > AMNIMAX) {
		printf("rayGetPlanes: nu=%d, > AMNIMAX=%d, TOO BIG!?!\n", nu, AMNIMAX);
		exit(EXIT_FAILURE);
	}

	while (set->next->next != set->next) { /* loop over the set of bundles: */
		set = set->next;
		listNodeCheck(set, "r.G.P: bundle");
		bundle = (struct Node *)set->item;
		if ((plane = (struct Plane *)listAlloc(1, sizeof(struct Plane))) == NULL) {
			printf("rayGetPlanes: failed to listAlloc plane struct!?!\n");
			exit(EXIT_FAILURE);
		}
		strcpy(plane->name, bundle->item);

		/* fit Zernike polynomials to pathlengths for this bundle: */
		for (i = 0; i < nu; i++) { /* zero normal equations */
			r[i] = 0.0;
			for (j = 0; j < nu; j++)
				A(i, j) = 0.0;
		}
		plane->n = 0;
		p = bundle;
		while (p->next->next != p->next) { /* loop over rays in bundle */
			p = p->next;
			listNodeCheck(p, "r.G.P: ray1");
			R = p->item;
			(plane->n)++;
			for (i = 0, radius = 0.0; i < 2; i++) {
				temp = (R->T[i + 1] - yz_axis[i]);
				radius += (temp * temp);
			}
			radius = sqrt(radius); /* get polar coordinates */
			cost = (R->T[1] - yz_axis[0]) / radius;
			sint = (R->T[2] - yz_axis[1]) / radius;
			nu = mathZernike((radius / yz_radius), /* get Zernike coefficients */
					 cost, sint, nord, u, ui);
			for (i = 0; i < nu; i++) { /* accumulate normal equations */
				for (j = i; j < nu; j++)
					A(i, j) += (u[i] * u[j]) * R->Intensity;
				r[i] += (R->Length * u[i]) * R->Intensity;
			}
		}
		for (i = 0; i < nu; i++) /* copy to lower triangular */
			for (j = i; j < nu; j++)
				A(j, i) = A(i, j);
		for (i = 0; i < nu; i++) /* compute correlation matrix */
			for (j = 0; j < nu; j++)
				c[i][j] = A(i, j) / sqrt(A(i, i) * A(j, j));
		if (mathSyminv2(a, nu, AMNIMAX)) /* get inverse of normal equations */
			printf("rayGetPlanes: mathSyminv2(a,%d) returned TRUE\n", nu);
		for (i = 0; i < nu; i++) {
			for (j = 0, plane->Amni[i] = 0.0; j < nu; j++) /* Amni = r * a^-1 */
				plane->Amni[i] += A(i, j) * r[j];
			plane->mni[i] = ui[i];
		}
		plane->nu = nu;
		plane->sig = Isum = 0.0; /* now to get RMS of pathlength: */
		p = bundle;
		while (p->next->next != p->next) { /* loop over rays in bundle */
			p = p->next;
			listNodeCheck(p, "r.G.P: ray2");
			R = p->item;
			for (i = 0, radius = 0.0; i < 2; i++) {
				temp = (R->T[i + 1] - yz_axis[i]);
				radius += (temp * temp);
			}
			radius = sqrt(radius); /* get polar coordinates */
			cost = (R->T[1] - yz_axis[0]) / radius;
			sint = (R->T[2] - yz_axis[1]) / radius;
			nu = mathZernike((radius / yz_radius), /* get Zernike coefficients*/
					 cost, sint, nord, u, ui);
			for (i = 0, zlength = 0.0; i < nu; i++)
				zlength += plane->Amni[i] * u[i];
			temp = (R->Length - zlength);
			plane->sig += (temp * temp) * R->Intensity;
			Isum += R->Intensity;
		}
		plane->sig = sqrt(plane->sig / Isum);
		for (i = 0; i < nu; i++) /* compute errors of Amni[] */
			plane->Asig[i] = plane->sig * sqrt(A(i, i));

		if (first_bundle && /* print solution for first bundle only */
		    (strcmp(print_mode, "verbose") == 0)) {
			printf("# First bundle solution <%s>:\n", plane->name);
			printf("# Coeffs: i  mni     name         A_mni[i]\n");
			for (i = 0; i < nu; i++) {
				for (j = 0, aberration = nomatch; j < Seidel_n; j++) {
					if (Seidel_list[j].mni == ui[i]) {
						aberration = Seidel_list[j].name;
						break;
					}
				}
				printf("#        %2d %04d = %-13s=%9.4f +/- %6.4f,\n", i, ui[i], aberration,
				       plane->Amni[i], plane->Asig[i]);
			}
			printf("#                    RMS of fit =%9.4f\n", plane->sig);
			printf("# Correl: i  mni"); /* print correlation matrix */
			for (j = 0; j < nu; j++)
				printf("  %04d", ui[j]);
			printf("\n");
			for (i = 0; i < nu; i++) {
				printf("#        %2d %04d", i, ui[i]);
				for (j = 0; j < nu; j++) {
					if (fabs(c[i][j]) > 0.01)
						printf("%6.2f", c[i][j]);
					else
						printf("      ");
				}
				printf("\n");
			}
			first_bundle = FALSE;
		}

		plane->yz_radius = yz_radius;
		listAppend(plane, planes);
	}
	return (planes);
}
