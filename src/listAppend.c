/* listAppend() -- add a new item to tail of a list 
   D.Wells, NRAO-CV
   1998-03-06: change tabs to blanks
*/
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

struct Node *listAppend (                      /* returns ptr to new node */
                         void *newitem,        /* ptr to new item struct  */
                         struct Node *list)    /* ptr to existing list    */
{
  struct Node *t;

  t = list;
  listNodeCheck(t, "l.A: list");
  if (t == NULL) {
      printf("listAppend: pointer to list is NULL!?!\n");
      exit(EXIT_FAILURE);
  }
  while (t->next->next != t->next) {
    listNodeCheck(t, "l.A: t_node");
    t = t->next;
  }
  return(listInsertAfter (newitem, t));
}
