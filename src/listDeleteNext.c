/* listDeleteNext() -- delete next node of list and associated malloc-ed data */

/* -=-=-=-=-=-=-=-=- Begin Copyright Notice -=-=-=-=-=-=-=-=-=-=-=-=-=-

   Copyright (C) 1998 Associated Universities, Inc. Washington DC, USA.

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

int listDeleteNext(                            /* return TRUE if node deleted */
		   struct Node *t)             /* ptr to node in a list */
{
    struct Node *n;

    if (t->next->next == t->next)
	return(FALSE);                         /* FALSE if no node to delete */
    n = t->next; 
    listNodeCheck(n, "l.D.N: n");
    listFree(n->item, "l.D.N: data item struct");
    t->next = n->next;                         /* cut Node out of chain      */
    listNodeDelete(n, "l.D.N: next");
    return(TRUE);                              /* TRUE if node deleted       */
}
