/* rayPltSystem() -- append Segments to show elements of System */

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

#define SETSEG(s,x1,y1,z1,x2,y2,z2,c) \
    s->T1[0] = x1; s->T1[1] = y1; s->T1[2] = z1; \
    s->T2[0] = x2; s->T2[1] = y2; s->T2[2] = z2; \
    s->ColorCode = c;

void rayPltSystem(                       /* returns lines in segments list */
		  struct Node *system,   /* list of optical elements       */
		  struct Node *segments) /* list of line segments          */
{
    struct Node *s;
    struct Surface *Surf;
    struct Segment *Seg;


   s = system;
    while (s->next->next != s->next) { /* loop over Surfaces */
	s = s->next;
	Surf = (struct Surface *)s->item;
	
    }

    /* Draw X,Y,Z axes: */
    Seg = (struct Segment *)listAlloc(1, sizeof(struct Segment));
    SETSEG(Seg,-65.,   0.,  0.,+8.,  0.,  0.,2001); /* X axis */
    listAppend(Seg, segments);
    Seg = (struct Segment *)listAlloc(1, sizeof(struct Segment));
    SETSEG(Seg,  0.,-105.,  0., 0.,+10.,  0.,2001); /* Y axis */
    listAppend(Seg, segments);
    Seg = (struct Segment *)listAlloc(1, sizeof(struct Segment));
    SETSEG(Seg,  0.,   0.,-55., 0.,  0.,+55.,2001); /* Z axis */
    listAppend(Seg, segments);
  
}
