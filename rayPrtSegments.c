/* rayPrtSegments() --- print a listing of coordinates of traced ray segments
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

void rayPrtSegments(struct Node *list, /* list of ray segments       */
		    int d)	       /* digits after decimal point */
{
	struct Node *t;
	struct Segment *s;
	int i;
	const int w = 10;
	char head_format[FORMMAX], body_format[FORMMAX];

	if ((t = list) == NULL) {
		printf("rayPrtSegments: list==NULL -->ABORT!\n");
		exit(EXIT_FAILURE);
	}
	printf("\n\t-=-< %s >-=-\n\n", (char *)t->item);
	i = 0;
	while (t->next->next != t->next) {
		t = t->next;
		s = (struct Segment *)t->item;
		i++;
		if (i == 1) {
			sprintf(head_format, "%%4s %%%ds%%%ds%%%ds %%%ds%%%ds%%%ds %%3s\n", w, w, w, w, w, w);
			sprintf(body_format, "%%4d %%%d.%dlf%%%d.%dlf%%%d.%dlf %%%d.%dlf%%%d.%dlf%%%d.%dlf %%3d\n", w,
				d, w, d, w, d, w, d, w, d, w, d);
			printf(head_format, "i", "T1[0] ", "T1[1] ", "T1[2] ", "T2[0] ", "T2[1] ", "T2[2] ", "Code");
			printf(head_format, "--", "-------", "-------", "-------", "-------", "-------", "-------",
			       "---");
		}
		printf(body_format, i, s->T1[0], s->T1[1], s->T1[2], s->T2[0], s->T2[1], s->T2[2], s->ColorCode);
	}
	if (i == 0) {
		printf("rayPrtSegments: list <%s> is empty.\n\n", (char *)t->item);
	} else {
		printf("\n");
	}
}
