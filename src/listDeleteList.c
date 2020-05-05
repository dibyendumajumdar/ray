/* listDeleteList() -- delete a list of nodes and malloc-ed items */

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

void listDeleteList (struct Node *list)        /* list of items */
{
  if (list == NULL) {
      printf("listDeleteList: pointer to list is NULL!?!\n");
      exit(EXIT_FAILURE);
  }
  listNodeCheck(list, "l.D.L: list");
  while (listDeleteNext(list)) {};      /* returns TRUE if delete successful */
  listNodeDelete(list->next, "l.D.L: tail");   /* delete the tail Node       */
  listFree(list->item, "l.D.L: head string");  /* delete the head string     */
  listNodeDelete(list, "l.D.L: head");         /* delete the head Node       */
}
