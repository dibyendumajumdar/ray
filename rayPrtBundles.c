/* rayPrtBundles() --- print a listing of lists of Ray structs
   D.Wells, NRAO-CV, Nov93/July95 */

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

void rayPrtBundles(struct Node *set,           /* list of lists of rays      */
		   int d)                      /* digits after decimal point */
{
  struct Node *s, *bundle, *t;
  struct Ray *r;
  int i, b;
  const int w=10, dc=8;
  char head_format[FORMMAX], body_format[FORMMAX];

  if (set == NULL) {
    printf("rayPrtBundles: set==NULL -->ABORT!\n");
    exit(EXIT_FAILURE);
  }
  printf("\n\t-=-< %s >-=-\n", (char *) set->item);
  s = set;
  b = 0;
  while (s->next->next != s->next) {
    s = s->next;
    if ((bundle = (struct Node *) s->item) == NULL) {
      printf("rayPrtBundles: bundle==NULL -->ABORT!\n");
      exit(EXIT_FAILURE);
    }
    t = bundle;
    b++;
    printf("\n\t-=-< %s >-=-\n\n", (char *) t->item);
    i = 0;
    while (t->next->next != t->next) {
      t = t->next;
      r = (struct Ray *) t->item;
      i++;
      if (i == 1) {
	sprintf(head_format,
	       "%%2s%%%ds%%%ds%%%ds %%%ds%%%ds%%%ds%%%ds%%7s%%4s\n",
	       w, w, w, dc, dc, dc, w);
	/* printf("rayPrtBundles: head=<%s>\n", head_format); */
	sprintf(body_format,
	       "%%2d%%%d.%dlf%%%d.%dlf%%%d.%dlf %%%d.5lf%%%d.5lf%%%d.5lf%%%d.%dlf%%7.4lf%%4d\n",
		w, d, w, d, w, d, dc, dc, dc, w, d);
	/* printf("rayPrtBundles: body=<%s>\n", body_format); */
	printf(head_format,
	       "i", "T[0] ", "T[1] ", "T[2] ",
	       "Q[0] ", "Q[1] ", "Q[2] ", "Length ", "Inten", " Color");
	printf(head_format,
	       "--", "-------", "-------", "-------",
	       "-------", "-------", "-------", "-------", "------", "---");
      }
      printf(body_format, i, r->T[0], r->T[1], r->T[2],
	     r->Q[0], r->Q[1], r->Q[2], r->Length, r->Intensity, r->ColorCode);
    }
    if (i == 0) {
      printf("rayPrtBundles: bundle <%s> contains no rays.\n\n",
	     (char *) t->item);
    } else {
      printf("\n");
    }
  }
  if (b == 0) printf("rayPrtBundles: set <%s> contains no bundles of rays.\n\n",
		     (char *)set->item);
}
