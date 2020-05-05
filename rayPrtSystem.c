/* rayPrtSystem() --- produce tables of parameters of a list of Surfaces
   D.Wells, NRAO-CV, Jan94/July95 */

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

void rayPrtSystem(struct Node *list, /* list of surfaces           */
		  int d)	     /* digits after decimal point */
{
	struct Node *t;
	struct Surface *s;
	int i;
	char head_format[FORMMAX], body_format[FORMMAX], temp[FORMMAX];

	if (list == NULL) {
		printf("\n\t-=-< System_list is NULL! >-=-\n\n");
		return;
	}
	t = list;
	printf("\n\t-=-< %s >-=-\n", (char *)t->item);

	i = 0; /* list names, curvatures, .. */
	while (t->next->next != t->next) {
		t = t->next;
		s = (struct Surface *)t->item;
		i++;
		if (i == 1) {
			strcpy(head_format, "%2s%14s");
			strcpy(body_format, "%2d%14s");
			strcat(head_format, " %12s%7s%8s%7s%9s\n");
			strcat(body_format, " %12.6lg%8.5lg%7.2lg%7.2lg%9.5lg\n");
			printf("Surface Properties:\n");
			printf(head_format, "i", "name ", "curv", "eps", "A_2", "A_4", "mu");
			printf(head_format, "--", "--------", "-----", "----", "----", "----", "----");
		}
		printf(body_format, i, s->name, s->c_1, s->eps, s->a_2, s->a_4, s->mu_1);
	}

	t = list; /* list vertex positions & tilts */
	i = 0;
	while (t->next->next != t->next) {
		t = t->next;
		s = (struct Surface *)t->item;
		i++;
		if (i == 1) {
			strcpy(head_format, "%2s%14s ");
			strcpy(body_format, "%2d%14s ");
			strcat(head_format, " %9s%9s%9s");
			sprintf(temp, " %%9.%dlf%%9.%dlf%%9.%dlf", d, d, d);
			strcat(body_format, temp);
			strcat(head_format, " %10s%10s%10s\n");
			strcat(body_format, " %10.5lf%10.5lf%10.5lf\n");
			printf("Vertex Positions and Tilts:\n");
			printf(head_format, "i", "name", "S[0]", "S[1]", "S[2]", "E[0]", "E[1]", "E[2]");
			printf(head_format, "--", "---------", "----", "----", "----", "-----", "-----", "-----");
		}
		printf(body_format, i, s->name, s->S[0], s->S[1], s->S[2], s->E[0], s->E[1], s->E[2]);
	}

	if (i == 0) {
		printf("listPrintSystem: <%s> is empty.\n\n", (char *)t->item);
	} else {
		printf("\n");
	}
}
