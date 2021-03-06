/* rayAddSurface.c -- append a new surface description to list of surfaces */

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

struct Node *rayAddSurface(				   /* returns ptr to new node     */
			   struct Node *list,		   /* list of surfaces            */
			   char surfname[],		   /* descriptive string          */
			   double curvature,		   /* 1/r                         */
			   double k,			   /* eccentricity of conic       */
			   double A_2,			   /* deformation terms           */
			   double A_4, double index_ratio, /* N/N_1, -1 means mirror      */
			   double S[],			   /* XYZ of vertex               */
			   double E[],			   /* Euler angles of vertex tilt */
			   double A_6, double A_8, double A_10, double A_12, double A_14
)
{
	int l;
	struct Surface *s;

	s = (struct Surface *)listAlloc(1, sizeof(struct Surface));
	l = strlen(surfname);
	if (l < (NAMEMAX - 1))
		strcpy(s->name, surfname);
	else {
		strncpy(s->name, surfname, (NAMEMAX - 2));
		s->name[NAMEMAX - 1] = '\0';
	}
	s->c_1 = curvature;
	s->k = k;
	s->a_2 = A_2;
	s->a_4 = A_4;
	s->a_6 = A_6;
	s->a_8 = A_8;
	s->a_10 = A_10;
	s->a_12 = A_12;
	s->a_14 = A_14;
	s->mu_1 = index_ratio;
	for (l = 0; l < 3; l++) {
		s->S[l] = S[l];
		s->E[l] = E[l];
	}
	return (listAppend(s, list));
}
