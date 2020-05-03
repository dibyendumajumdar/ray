/* rayPrtFoci() --- print a list of foci results
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

void rayPrtFoci (struct Node *foci_list,       /* list of Focus structs      */
		 int d)                        /* digits after decimal point */
{
    struct Node *bundle;
    struct Focus *f;
    char head[FORMMAX], body[FORMMAX], temp[FORMMAX];
    int i, w;
    
    if (foci_list == NULL) {
	printf("rayPrtFoci: foci_list is NULL!?!\n");
	exit(EXIT_FAILURE);
    }
    bundle = foci_list;
    printf("\n\t\t-=-< %s >-=-\n", (char *) bundle->item);
    i = 0;
    /* loop over the set of bundles: */
    while (bundle->next->next != bundle->next) {
	bundle = bundle->next;
	f = (struct Focus *) bundle->item;
	i++;
	if (i == 1) {
	    w = d + 5;
	    sprintf(head, "%%3s%%%ds%%4s%%%ds%%%ds%%%ds%%%ds",
		    NAMEMAX, w, w, w, w);
	    sprintf(body, "%%3d%%%ds%%4d%%%d.%dlf%%%d.%dlf%%%d.%dlf%%%d.%dlf",
		    NAMEMAX, w, d, w, d, w, d, w, d);
	    w = d + 3; 
	    sprintf(temp, "%%%ds%%%ds%%%ds%%%ds\n",
		    w, w, w, w);
	    strcat(head, temp);
	    sprintf(temp, "%%%d.%dlf%%%d.%dlf%%%d.%dlf%%%d.%dlf\n",
		    w, d, w, d, w, d, w, d);
	    strcat(body, temp);
	    printf(head, "i", "bundle_name", "n",
		   "xc", "yc", "zc", "lc",
		   "xs", "ys", "zs", "ls");
	    printf(head, "--", "----------------", "---",
		   "------", "------", "------", "-------",
		   "-----", "-----", "-----", "-----"); 
	}
	printf (body, i, f->name, f->n,
		f->xyz[0], f->xyz[1], f->xyz[2], f->lc,
		f->xyzs[0], f->xyzs[1], f->xyzs[2], f->ls);
    }
    if (i == 0) {
	printf("rayPrtFoci: list <%s> is empty.\n\n", (char *) bundle->item);
    } else {
	printf("\n");
    }
}
