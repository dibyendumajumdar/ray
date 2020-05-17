/* rayPrtPlanes() --- print a list of Plane results
   D.Wells, NRAO-CV
   1993-11: initial version
   1995-07: mods
   1996-01-05: more mods
   1997-11-21: change Seidel tilt term indicies 110&111-->1010&1011
   1998-10-23: delete unused variable 'char *indicies[AMNIMAX]'
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

#define SEIDEL 1
#include "ray.h"

void rayPrtPlanes(struct Node *planes_list, /* list of Plane structs      */
		  int d)		    /* digits after decimal point */
{
	struct Node *bundle;
	struct Plane *p;
	char head[FORMMAX], body[FORMMAX], temp[FORMMAX], form[FORMMAX], temp2[FORMMAX], *match;
	const char *aberrations[AMNIMAX];
	double small, val;
	int i, j, k, w;
	struct {
		const char *original;
		const char *replace;
	} zeros[] = {{"-0.", " -."}, {" 0.", "  ."}};

	if (planes_list == NULL) {
		printf("rayPrtPlanes: planes_list is NULL!?!\n");
		exit(EXIT_FAILURE);
	}
	bundle = planes_list;
	printf("\n\t\t-=-< %s >-=-\n", (char *)bundle->item);
	small = 0.5 * pow(10.0, -(double)d);
	i = 0;
	/* loop over the set of bundles: */
	while (bundle->next->next != bundle->next) {
		bundle = bundle->next;
		p = (struct Plane *)bundle->item;
		i++;
		if (i == 1) {
			w = d + 5;
			sprintf(head, "%%3s%%%ds%%3s%%%ds", NAMEMAX, w);
			sprintf(body, "%%3d%%%ds%%3d%%%d.%dlf", NAMEMAX, w, d);
			w = d + 3;
			sprintf(temp, "A_%04d", p->mni[0]);
			printf(head, " ", " ", " ", temp);
			sprintf(form, "%%%ds", w);
			for (j = 1; j < p->nu; j++) {
				sprintf(temp, "A_%03d", p->mni[j]);
				printf(form, temp);
			}
			printf("\n");
			for (j = 0; j < p->nu; j++) {
				for (k = 0, aberrations[j] = "-na-"; k < Seidel_n; k++) {
					if (Seidel_list[k].mni == p->mni[j]) {
						aberrations[j] = Seidel_list[k].name;
						break;
					}
				}
			}
			printf(head, "i", "bundle_name    ", "n", aberrations[0]);
			for (j = 1; j < p->nu; j++) {
				strcpy(temp, " ");
				strncat(temp, aberrations[j], w - 1);
				temp[w] = '\0';
				printf(form, temp);
			}
			printf("\n");
			printf(head, "--", "--------------------", "--", "-------");
			for (j = 1; j < p->nu; j++)
				printf(form, "-----");
			printf("\n");
			sprintf(form, "%%%d.%dlf", w, d);
		}
		printf(body, i, p->name, p->n, p->Amni[0]);
		snprintf(temp, sizeof temp, "%%%ds", w);
		for (j = 1; j < p->nu; j++) {
			val = p->Amni[j];
			/* convert the Seidel wave tilt terms to milliradians: */
			if ((p->mni[j] == 1010) || (p->mni[j] == 1011))
				val *= (1000.0 / p->yz_radius);
			if (fabs(val) > small) {
				snprintf(temp2, sizeof temp2, form, val);
				for (k = 0; k < 2; k++) { /* delete redundant zeroes */
					if ((match = strstr(temp2, zeros[k].original)) != NULL) {
						strncpy(match, zeros[k].replace, 3);
						break;
					}
				}
				printf(temp, temp2);
			} else
				printf(temp, " ");
		}
		printf("\n");
		/* print RMS wavefront error if it is significant: */
		if (p->sig > small) {
			sprintf(temp2, "%%%d.%dlf", w, d);
			sprintf(temp, temp2, p->sig);
			printf(head, " ", "NOTE! RMS residuals", "=", temp);
			printf("\n");
		}
	}
	if (i == 0) {
		printf("rayPrtPlanes: list <%s> is empty.\n\n", (char *)bundle->item);
	} else {
		printf("\n");
	}
}
